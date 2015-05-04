# vim: fdm=marker
'''
author:     Fabio Zanini/Richard Neher
date:       25/04/2015
content:    Data access module HIV patients.
'''
# Modules
import numpy as np
import pandas as pd
from samples import *
from Bio import SeqIO


# Classes
class Patient(pd.Series):
    '''HIV patient'''

    def __init__(self, *args, **kwargs):
        '''Initialize a patient with all his samples'''
        include_cell = kwargs.pop('include_cell', False)
        super(Patient, self).__init__(*args, **kwargs)
        self.samples = sorted(load_samples_sequenced(patients=[self.name], include_cell=include_cell), 
                              key = lambda x:x.date)
        self._dates = [x.date for x in self.samples]
        self._cd4 = [x['CD4+ count'] for x in self.samples]
        self._viral_load = [x['viral load'] for x in self.samples]

        # We take 400 ul of serum
        # We typically have 6 reactions with that total volume (plus the F4 dilution
        # series, but each of those uses only 0.1x template which is very little)
        self._n_templates_viral_load = np.array([x*0.4/6.1 for x in self._viral_load], dtype = float)
        #self._n_templates_dilutions = np.ma.masked_invalid([x.get_n_template_dilutions() for x in self.samples])
        self._times = []
        self.reference = self.load_reference()
        self.annotation = {x.qualifiers['note'][-1]:x for x in self.reference.features}

    @classmethod
    def load(cls, pname):
        from hivwholeseq.sequencing.filenames import table_filename
        patients = pd.read_excel(table_filename, 'Patients', index_col=1)
        patients.index = pd.Index(map(str, patients.index))
        if pname in patients.index:
            return cls(patients.loc[pname])
        else:
            return cls(patients.loc[patients.code == pname].iloc[0])

    @property
    def _constructor(self):
        return Patient

    @property
    def transmission_date(self):
        '''The most likely time of transmission'''
        return self['last negative date'] + \
                (self['first positive date'] - self['last negative date']) / 2


    @property
    def dates(self):
        return self._dates
    @property
    def viral_load(self):
        return self._viral_load
    @property
    def cd4(self):
        return self._cd4

    @property
    def dsi(self):
        return self.times(unit='days')
    @property
    def msi(self):
        return self.times(unit='month')
    @property
    def ysi(self):
        return self.times(unit='year')
    
    
    def times(self, unit='days'):
        '''Get the times from transmission'''
        delta_days = [d.toordinal() - self.transmission_date.toordinal() for d in self.dates]
        if unit.startswith('day'):
            return np.array(delta_days, dtype=float)
        elif unit.startswith('year'):
            return np.array(delta_days, dtype=float)/365.24
        elif unit.startswith('month'):
            return np.array(delta_days, dtype=float)/365.24*12
        else:
            raise ValueError("bad time unit")

    @property
    def n_templates_dilutions(self):
        '''Get the time course of the number of templates to PCR, limiting depth'''
        return self._n_templates_dilutions

    @property
    def n_templates_viral_load(self):
        '''Get the number of templates, estimated from the viral load'''
        return self._n_templates_viral_load


    def get_n_templates_roi(self, roi):
        '''Get number of templates, roi specific from overlap frequencies'''
        # FIXME: this loads stuff from get_roi, which then goes back to patient
        fragments = self.get_fragments_covered(roi)
        n = [min(sample[fr+'q'] for fr in fragments)
             for _, sample in self.samples.iterrows()]
        n = np.ma.masked_invalid(n)
        return n

    @property
    def initial_sample(self):
        '''The initial sample used as a mapping reference'''
        return self.samples[0]


    def load_reference(self):
        from hivwholeseq.patients.filenames import get_initial_reference_filename
        return SeqIO.read(get_initial_reference_filename(self.name, "genomewide", format='gb'), 'gb')


    def _annotation_to_fragment_indices(self, anno):
        coordinates = {}
        ind = [int(x) for x in self.annotation[anno]]
        coordinates['start'] = min(ind)
        coordinates['length'] = len(ind)
        fragments = ['F'+str(i) for i in xrange(1,7)]
        if anno not in fragments:
            for frag in fragments:
                frag_ind = set([int(x) for x in self.annotation[frag]])
                anno_indices_on_fragment = sorted(frag_ind.intersection(ind))
                if len(anno_indices_on_fragment):
                    anno_indices_self = np.arange(coordinates['length'])[np.in1d(ind, anno_indices_on_fragment)]
                    coordinates[frag] = (anno_indices_self, 
                                       np.array(anno_indices_on_fragment)- int(self.annotation[frag].location.start))
        else:
            coordinates[anno] = (np.arange(coordinates['length']), np.arange(coordinates['length']))
        return coordinates

    def get_coverage_trajectories(self, region, **kwargs):
        '''Get coverage as a function of time'''
        coordinates = self._annotation_to_fragment_indices(region)
        cov = np.ma.array([tmp_sample.get_coverage(coordinates, **kwargs) for tmp_sample in self.samples])
        return cov

    def get_allele_count_trajectories(self, region, safe=False, **kwargs):
        '''Get the allele count trajectories from files
        
        Args:
          region (str): region to study, a fragment or a genomic feature (e.g. V3)
          **kwargs: passed down to the function (VERBOSE, etc.).

        Note: the genomewide counts are currently saved to file.
        '''
        coordinates = self._annotation_to_fragment_indices(region)
        act = np.ma.array([tmp_sample.get_allele_counts(coordinates, **kwargs) for tmp_sample in self.samples])
        # set very low frequencies to zero, these are likely sequencing errors
        return act

    def get_allele_frequency_trajectories(self, region, safe=False,error_rate = 2e-3,  **kwargs):
        '''Get the allele count trajectories from files
        
        Args:
          region (str): region to study, a fragment or a genomic feature (e.g. V3)
          **kwargs: passed down to the function (VERBOSE, etc.).

        Note: the genomewide counts are currently saved to file.
        '''
        coordinates = self._annotation_to_fragment_indices(region)
        aft = np.ma.array([tmp_sample.get_allele_frequencies(coordinates, **kwargs) for tmp_sample in self.samples])
        # set very low frequencies to zero, these are likely sequencing errors
        aft[aft<error_rate]=0
        return aft

    @staticmethod
    def get_initial_consensus_noinsertions(aft, VERBOSE=0, return_ind=False):
        '''Make initial consensus from allele frequencies, keep coordinates and masked
        
        Args:
          aft (np.ma.ndarray): 3d masked array with the allele frequency trajectories

        Returns:
          np.ndarray: initial consensus, augmented with later time points at masked
          positions, with Ns if never covered
        '''
        from ..utils.sequence import alpha

        af0 = aft[0]
        # Fill the masked positions with N...
        cons_ind = af0.argmax(axis=0)
        cons_ind[af0[0].mask] = 5
    
        # ...then look in later time points
        if aft.shape[0] == 1:
            if return_ind:
                return cons_ind
            else:
                return alpha[cons_ind]

        for af_later in aft[1:]:
            cons_ind_later = af_later.argmax(axis=0)
            cons_ind_later[af_later[0].mask] = 5
            ind_Ns = (cons_ind == 5) & (cons_ind_later != 5)
            cons_ind[ind_Ns] = cons_ind_later[ind_Ns]
        if return_ind:
            return cons_ind
        else:
            return alpha[cons_ind]

    def get_map_coordinates_reference(self, roi, refname='HXB2', in_patient = True):
        from hivwholeseq.patients.filenames import get_coordinate_map_filename
        genomewide_map = np.loadtxt(get_coordinate_map_filename(self.name, 'genomewide', refname=refname), dtype = int)
        if roi in self.annotation:
            roi_pos = np.array([x for x in self.annotation[roi]], dtype = int)
            ind = np.in1d(genomewide_map[:,1], roi_pos)
            roi_indices = np.in1d(roi_pos, genomewide_map[:,1]).nonzero()[0]
            return np.vstack((genomewide_map[ind].T, [roi_indices])).T
        elif roi == "genomewide":
            return np.vstack((genomewide_map.T, [genomewide_map[:,1]])).T            
        else:
            try:
                start, stop = map(int, roi)
                start_ind = np.searchsorted(genomewide_map[:,in_patient], start)
                stop_ind = np.searchsorted(genomewide_map[:,in_patient], stop)
                return np.vstack((genomewide_map[start_ind:stop_ind].T, [genomewide_map[start_ind:stop_ind,in_patient]-start])).T
            except:
                raise ValueError("ROI not understood")

    def positions_to_features(self):
        '''
        map of positions to features, including the number of proteins, RNA, etc this pos is part of
        '''
        self.pos_to_feature = [{'protein':0, 'RNA':0, 'LTR':0, 'codons':[]} 
                                for pos in xrange(len(self.reference))]
        for fname, feature in self.annotation.iteritems():
            for ii, pos in enumerate(feature):
                if feature.type=='protein':
                    self.pos_to_feature['protein']+=1
                    self.codons.append((fname, ii//3, pos%3))
                elif feature.type=='RNA':
                    self.pos_to_feature['RNA']+=1
                elif 'LTR' in fname:
                    self.pos_to_feature['LTR']+=1

if __name__=="__main__":
    from matplotlib import pyplot as plt
    plt.ion()
    p = Patient.load('20097')
    p = Patient.load('p3')

