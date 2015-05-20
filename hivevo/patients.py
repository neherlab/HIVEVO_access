# vim: fdm=marker
'''
author:     Fabio Zanini/Richard Neher
date:       25/04/2015
content:    Data access module HIV patients.
'''
# Modules
import numpy as np
import pandas as pd
from Bio import SeqIO
from .samples import *
from .af_tools import *
from .sequence import alpha, alphaa



# Classes
class Patient(pd.Series):
    '''
    Class providing access to longitudinal sequencing data of HIV-1 populations
    in participants of the HIVEVO study. The class contains time-ordered samples
    and access methods to single nucleotide variants, pair frequencies, and genomic
    features of the the HIV poputions
    '''

    def __init__(self, *args, **kwargs):
        '''Initialize a patient with all his samples'''
        include_cell = kwargs.pop('include_cell', False)
        super(Patient, self).__init__(*args, **kwargs)
        self.samples = sorted(load_samples_sequenced(patients=[self.name]), 
                              key = lambda x:x.date)
        self._dates = [x.date for x in self.samples]
        self._cd4 = [x['CD4+ count'] for x in self.samples]
        self._viral_load = [x['viral load'] for x in self.samples]

        # We take 400 ul of serum
        # We typically have 6 reactions with that total volume (plus the F4 dilution
        # series, but each of those uses only 0.1x template which is very little)
        self._n_templates_viral_load = np.array([x*0.4/6.1 for x in self._viral_load], dtype = float)
        self._n_templates_dilutions = np.ma.masked_invalid([x.get_n_templates_dilutions() for x in self.samples])
        self._times = []
        self.reference = self.load_reference()
        # translate genbank encoded sequence features into a dictionary
        self.annotation = {x.qualifiers['note'][-1]:x for x in self.reference.features}
        self._initial_consensus_noinsertions()

    @classmethod
    def load(cls, pname):
        from .filenames import table_filename
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
        return self['last negative date'] + (self['first positive date'] - self['last negative date']) / 2

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

    @property
    def initial_sample(self):
        '''The initial sample used as a mapping reference'''
        return self.samples[0]

    def load_reference(self):
        from .filenames import get_initial_reference_filename
        return SeqIO.read(get_initial_reference_filename(self.name, "genomewide", format='gb'), 'gb')

    def _region_to_indices(self,region):
        '''returns a list of positions corresponding to a genomic region'''
        if region=='genomewide':
            return np.arange(len(self.reference))
        elif region in self.annotation:
            return np.array([int(x) for x in self.annotation[region]], dtype=int)
        else:
            raise ValueError('no annotation with name '+region)

    def _annotation_to_fragment_indices(self, anno):
        '''
        returns coordinates of a region specified in the annotation
        in terms of the fragments F1 to F5. This is needed to extract 
        region specific allele counts, frequencies etc.
        returns a dict containing 'length', 'start' (of the region in the genome)
        and for each fragment involved 'F1': (indices in the region of interest, indices on the fragment)
        '''
        coordinates = {}
        region_indices = self._region_to_indices(anno)
        coordinates['start'] = min(region_indices)
        coordinates['length'] = len(region_indices)
        fragments = ['F'+str(i) for i in xrange(1,7)]
        if anno not in fragments:
            for frag in fragments: # loop over fragments and extract the indices of the region on this fragment
                frag_ind = set(self._region_to_indices(frag))       # indices of the fragment
                region_indices_on_fragment = sorted(frag_ind.intersection(region_indices))  # intersection of region and fragment positions
                if len(region_indices_on_fragment): # attach indices in region and on fragment
                    anno_indices_self = np.arange(coordinates['length'])[np.in1d(region_indices, region_indices_on_fragment)]
                    coordinates[frag] = (anno_indices_self, 
                                       np.array(region_indices_on_fragment)- int(self.annotation[frag].location.start))
        else: # if requested region is a fragment, return only this fragment
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
        act = np.ma.array([tmp_sample.get_allele_counts(coordinates, **kwargs)
                           for tmp_sample in self.samples], hard_mask=True, shrink=False)
        return act

    def get_allele_frequency_trajectories(self, region, safe=False,error_rate = 2e-3,  **kwargs):
        '''Get the allele count trajectories from files
        
        Args:
          region (str): region to study, a fragment or a genomic feature (e.g. V3)
          **kwargs: passed down to the function (VERBOSE, etc.).

        Note: the genomewide counts are currently saved to file.
        '''
        coordinates = self._annotation_to_fragment_indices(region)
        aft = np.ma.array([tmp_sample.get_allele_frequencies(coordinates, **kwargs) 
                          for tmp_sample in self.samples], hard_mask=True, shrink=False)
        # set very low frequencies to zero, these are likely sequencing errors
        aft[aft<error_rate]=0
        return aft

    def _initial_consensus_noinsertions(self, VERBOSE=0, return_ind=False):
        '''Make initial consensus from allele frequencies, keep coordinates and masked
        sets: indices and sequence of initial sequence
        '''
        aft = self.get_allele_frequency_trajectories('genomewide')
        # Fill the masked positions with N...
        cons_ind = aft[0].argmax(axis=0)
        cons_ind[aft.mask[0].max(axis=0)] = 5
    
        for af_later in aft[1:]:
            cons_ind_later = af_later.argmax(axis=0)
            cons_ind_later[af_later.mask.max(axis=0)] = 5
            ind_Ns = (cons_ind == 5) & (cons_ind_later != 5)
            cons_ind[ind_Ns] = cons_ind_later[ind_Ns]

        self.initial_indices = cons_ind
        self.initial_sequence = alpha[cons_ind]

    def get_diversity(self, region):
        aft = self.get_allele_frequency_trajectories(region)
        return np.array(map(diversity, aft))

    def get_consensi(self, region):
        aft = self.get_allele_frequency_trajectories(region)
        return [''.join(consensus(x)) for x in aft]

    def get_divergence(self, region):
        aft = self.get_allele_frequency_trajectories(region)
        region_initial_indices = self.initial_indices[self._region_to_indices(region)]
        return np.array([divergence(x,region_initial_indices) for x in aft])

    def map_to_external_reference(self, roi, refname='HXB2', in_patient = True):
        '''
        return a map of positions in the patient to a reference genomewide
        Args:
            roi  --  region of interest given as a string or a tuple (start, end)
            refname --  reference to compare to
            in_patient -- specifies whether the (start, end) refers to reference or patient coordinates
        returns:
            a (len(roi), 3) array with reference coordinates in first column, 
                                        patient coordinates in second 
                                        roi coordinates in third column
        '''
        from .filenames import get_coordinate_map_filename
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

    # TODO: the following is experimental. was meant as a way to easily get an idea what kind of stuff a site is involved in
    def positions_to_features(self):
        '''
        map of positions to features, including the number of proteins, RNA, etc this pos is part of
        '''
        self.pos_to_feature = [{'gene':0, 'RNA':0, 'LTR':0, 'codons':[], 'protein_codon':[]} 
                                for pos in xrange(len(self.reference))]
        for fname, feature in self.annotation.iteritems():
            for ii, pos in enumerate(feature):
                if feature.type=='gene':
                    self.pos_to_feature[pos]['gene']+=1
                    self.pos_to_feature[pos]['codons'].append((fname, ii//3, ii%3))
                elif feature.type=='protein':
                    self.pos_to_feature[pos]['protein_codon'].append((fname, ii//3, ii%3))
                elif 'LTR' in fname:
                    self.pos_to_feature[pos]['LTR']+=1
                elif feature.type=='RNA_structure':
                    self.pos_to_feature[pos]['RNA']+=1

    def get_fragment_depth(self, pad=False, limit_to_dilution = False):
        c = self._annotation_to_fragment_indices('genomewide')
        depth = np.ma.array([s.fragment_depth(c,cov_min = 100, var_min = 0.05, min_points = 10) for s in self.samples])
        if pad:
            for si in xrange(len(self.samples)):
                depth[si][depth.mask[si]] = self.n_templates_dilutions[si]
                depth.mask[si] = False
        if limit_to_dilution:
            for si in xrange(len(self.samples)):
                depth[si] = np.minimum(depth[si], self.n_templates_dilutions[si])
        return depth
