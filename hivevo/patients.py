# vim: fdm=marker
'''
author:     Fabio Zanini/Richard Neher
date:       25/04/2015
content:    Data access module HIV patients.
'''
# Modules
from collections import Counter
from operator import itemgetter
import numpy as np
import pandas as pd
from Bio import SeqIO,Seq
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
    _metadata = ['samples',
                 '_cd4',
                 '_viral_load',
                 '_times',
                 '_n_templates_dilutions',
                 '_n_templates_viral_load',
                 'initial_sequence',
                 'pos_to_feature',
                 'reference',
                 'annotation',
                ]

    def __init__(self, *args, **kwargs):
        '''Initialize a patient with all his samples'''
        include_cell = kwargs.pop('include_cell', False)
        super(Patient, self).__init__(*args, **kwargs)
        self.samples = sorted(load_samples_sequenced(patients=[self.name]),
                              key=itemgetter('days since infection'))
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
        self.annotation = {x.qualifiers['note'][-1]: x for x in self.reference.features}
        self._initial_consensus_noinsertions()

        # TODO: this is experimental
        self.positions_to_features(sources=['annotations', 'shape', 'disorder', 'accessibility', 'structural'])


    @classmethod
    def load(cls, pname):
        from .filenames import get_table_filename
        patients = pd.read_csv(get_table_filename('patients'),
                               sep='\t',
                               index_col=0)
        if pname not in patients.index:
            raise ValueError('Patient '+str(pname)+' not found')
        return cls(patients.loc[pname])


    @property
    def _constructor(self):
        return Patient


    @property
    def _constructor_expanddim(self):
        # TODO: implement pd.DataFrame for several patients at once
        raise NotImplementedError


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
        delta_days = [float(s['days since infection']) for s in self.samples]
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
        fragments = ['F'+str(i) for i in range(1,7)]

        # if requested region is a fragment, return only this fragment
        if anno in fragments:
            coordinates[anno] = (np.arange(coordinates['length']), np.arange(coordinates['length']))

        # loop over fragments and extract the indices of the region on this fragment
        else:
            for frag in fragments:

                # indices of the fragment
                frag_ind = set(self._region_to_indices(frag))

                # intersection of region and fragment positions
                region_indices_on_fragment = sorted(frag_ind.intersection(region_indices))

                # attach indices in region and on fragment
                if len(region_indices_on_fragment):
                    anno_indices_self = np.arange(coordinates['length'])[np.in1d(region_indices, region_indices_on_fragment)]
                    coordinates[frag] = (anno_indices_self,
                                       np.array(region_indices_on_fragment)- int(self.annotation[frag].location.start))
        return coordinates


    def get_coverage_trajectories(self, region, **kwargs):
        '''Get coverage as a function of time'''
        coordinates = self._annotation_to_fragment_indices(region)
        cov = np.ma.array([tmp_sample.get_coverage(coordinates, **kwargs) for tmp_sample in self.samples])
        return cov


    def get_allele_count_trajectories(self, region, safe=False, type='nuc', **kwargs):
        '''Get the allele count trajectories from files

        Args:
          region (str): region to study, a fragment or a genomic feature (e.g. V3)
          type (str): 'nuc' for nucleotides, 'aa' for amino acids
          **kwargs: passed down to the function (VERBOSE, etc.).

        Note: the genomewide counts are currently saved to file.

        Examples:

        1. Nucleotides from a certain region:

        patient.get_allele_count_trajectories('PR', type='nuc')

        2. Amino acids from a certain region:

        patient.get_allele_count_trajectories('PR', type='aa')
        '''
        if type == 'nuc':
            coordinates = self._annotation_to_fragment_indices(region)
        elif type == 'aa':
            coordinates = {region: None}
        else:
            raise ValueError('Data type not understood')

        act = np.ma.array([tmp_sample.get_allele_counts(coordinates, type=type, **kwargs)
                           for tmp_sample in self.samples],
                          hard_mask=True,
                          shrink=False)
        if len(act.mask.shape) < 1:
            act.mask = np.zeros_like(act, dtype=bool)
        return act


    def get_insertion_trajectories(self, region, **kwargs):
        '''Get insertion trajectories

        Returns:
           ict (pd.Series): multiindexed array of insertion counts

        To manipulate the result, pandas methods are recommended, e.g.

        ict_matrix = ict.unstack('DSI').fillna(0).T
        '''
        coordinates = self._annotation_to_fragment_indices(region)
        ict = Counter()
        for tmp_sample in self.samples:
            time = tmp_sample['days since infection']
            ic = tmp_sample.get_insertions(coordinates, **kwargs)
            for (position, insertion), value in ic.items():
                ict[(time, position, insertion)] = value
        ict = pd.Series(ict, name='insertions')
        if len(ict):
            ict.index.names = ['DSI', 'position', 'insertion', "accessibility"]
        return ict


    def get_allele_frequency_trajectories(self, region, safe=False,error_rate = 2e-3, type='nuc',  **kwargs):
        '''Get the allele count trajectories from files

        Args:
          region (str): region to study, a fragment or a genomic feature (e.g. V3)
          type (str): 'nuc' for nucleotides, 'aa' for amino acids
          **kwargs: passed down to the function (VERBOSE, etc.).

        Note: the genomewide counts are currently saved to file.
        Examples:

        1. Nucleotides from a certain region:

        patient.get_allele_frequencies_trajectories('PR', type='nuc')

        2. Amino acids from a certain region:

        patient.get_allele_frequencies_trajectories('PR', type='aa')

        '''
        if type == 'nuc':
            coordinates = self._annotation_to_fragment_indices(region)
        elif type == 'aa':
            coordinates = {region: None}
        else:
            raise ValueError('Data type not understood')

        aft = np.ma.array([tmp_sample.get_allele_frequencies(coordinates, type=type, **kwargs)
                          for tmp_sample in self.samples], hard_mask=True, shrink=False)
        # set very low frequencies to zero, these are likely sequencing errors
        aft[aft<error_rate]=0
        if len(aft.mask.shape)<1:
            aft.mask = np.zeros_like(aft, dtype=bool)
        return aft


    def get_constrained(self, region):
        if region in self.annotation and self.annotation[region].type in ['gene', 'protein']:
            return np.array([self.pos_to_feature[pos]['RNA']>0 \
                        or self.pos_to_feature[pos]['gene']>1
                        for pos in self.annotation[region]])
        else:
            print(region, "is not a valid protein or gene")
            return None


    def get_gaps_by_codon(self, region, pad=0, threshold=0.1):
        if region in self.annotation and self.annotation[region].type in ['gene', 'protein']:
            aft = self.get_allele_frequency_trajectories(region)
            gap_index = list(alpha).index('-')
            gaps = np.zeros(aft.shape[-1],dtype=bool)
            for ci in range(0, aft.shape[-1],3):
                if np.any(aft[:,gap_index,ci:ci+3]>threshold):
                    gaps[max(0,ci-3*pad):ci+3*(1+pad)]=True
            return gaps
        else:
            print(region,"is not a valid protein or gene")
            return None


    def get_syn_mutations(self, region, mask_constrained=True):
        
        if region in self.annotation and self.annotation[region].type in ['gene', 'protein']:
            try:
                aft = self.get_allele_frequency_trajectories(region)
                if len(aft.mask.shape) == 0:
                    aft_valid = np.ones((aft.shape[0], aft.shape[-1]), dtype=bool)
                else:
                    aft_valid = -np.array([af.mask.sum(axis=0) for af in aft], dtype=bool)
                gaps = self.get_gaps_by_codon(region)
                initial_seq = self.get_initial_sequence(region)
                consensi = []
                for af in aft:
                    tmp = consensus(af)
                    tmp[gaps]='N'
                    consensi.append(tmp)

                cons_aa = np.array([np.fromstring(Seq.translate(''.join(cons)),
                                   dtype='|S1') for cons in consensi])
                no_substitution = np.repeat(np.array([len(np.unique(col[ind]))==1
                                for ind, col in zip(aft_valid.T[::3], cons_aa.T)], dtype=bool), 3)

                syn_muts = np.zeros(aft.shape[1:], dtype=bool)
                for pos in range(aft.shape[-1]):
                    ci = pos//3
                    rf = pos%3
                    codon = ''.join(initial_seq[ci*3:(ci+1)*3])
                    for ni,nuc in enumerate(alpha[:4]):
                        mod_codon = codon[:rf] + nuc + codon[rf+1:]
                        try:
                            syn_muts[ni,pos] = (Seq.translate(codon)==Seq.translate(mod_codon))\
                                                *no_substitution[pos]
                        except:
                            syn_muts[ni,pos] = False
                if mask_constrained:
                    syn_muts[:,self.get_constrained(region)] = False
                return syn_muts
            except:
                import pdb; pdb.set_trace()
        else:
            print(region,"is not a valid protein or gene")
            return None


    def _initial_consensus_noinsertions(self, region='genomewide',VERBOSE=0, type='nuc'):
        '''Make initial consensus from allele frequencies, keep coordinates and masked
        sets: indices and sequence of initial sequence
        '''
        aft = self.get_allele_frequency_trajectories(region, type=type)
        tmp_alpha = alpha if type=='nuc' else alphaa
        # Fill the masked positions with N...
        mask_index = len(tmp_alpha)-1
        cons_ind = aft[0].argmax(axis=0)
        cons_ind[aft.mask[0].max(axis=0)] = mask_index

        for af_later in aft[1:]:
            cons_ind_later = af_later.argmax(axis=0)
            cons_ind_later[af_later.mask.max(axis=0)] = mask_index
            ind_Ns = (cons_ind == mask_index) & (cons_ind_later != mask_index)
            cons_ind[ind_Ns] = cons_ind_later[ind_Ns]

        if region == 'genomewide' and type == 'nuc':
            self.initial_indices = cons_ind
            self.initial_sequence = alpha[cons_ind]
        else:
            return cons_ind

    def get_initial_indices(self, region, type='nuc'):
        if type=='nuc':
            if region == 'genomewide':
                return self.initial_indices
            elif region in self.annotation:
                return np.array([self.initial_indices[pos] for pos in self.annotation[region]])
            else:
                print("Not a valid annotation:", region)
                return None
        elif type=='aa':
            return self._initial_consensus_noinsertions(region, type=type)

    def get_initial_sequence(self, region, type='nuc'):
        tmp_ind = self.get_initial_indices(region, type=type)
        tmp_alpha = alpha if type=='nuc' else alphaa
        if tmp_ind is not None:
            return tmp_alpha[tmp_ind]
        else:
            return None


    def get_diversity(self, region):
        aft = self.get_allele_frequency_trajectories(region)
        return np.array(list(map(diversity, aft)))


    def get_consensi(self, region):
        aft = self.get_allele_frequency_trajectories(region)
        return [''.join(consensus(x)) for x in aft]


    def get_divergence(self, region):
        aft = self.get_allele_frequency_trajectories(region)
        region_initial_indices = self.initial_indices[self._region_to_indices(region)]
        return np.array([divergence(x,region_initial_indices) for x in aft])


    def map_to_external_reference(self, roi, refname='HXB2', in_patient=True):
        '''
        Map positions in the patient to a reference genomewide
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
        coo_fn = get_coordinate_map_filename(self.name, refname=refname)
        genomewide_map = np.loadtxt(coo_fn, dtype=int)

        if roi in self.annotation:
            roi_pos = np.array([x for x in self.annotation[roi]], dtype=int)
            ind = np.in1d(genomewide_map[:,1], roi_pos)
            roi_indices = np.in1d(roi_pos, genomewide_map[:,1]).nonzero()[0]
            return np.vstack((genomewide_map[ind].T, [roi_indices])).T

        elif roi == "genomewide":
            return np.vstack((genomewide_map.T, [genomewide_map[:,1]])).T

        else:
            try:
                start, stop = list(map(int, roi))
                start_ind = np.searchsorted(genomewide_map[:,in_patient], start)
                stop_ind = np.searchsorted(genomewide_map[:,in_patient], stop)
                return np.vstack((genomewide_map[start_ind:stop_ind].T,
                                  [genomewide_map[start_ind:stop_ind, in_patient] - start])).T
            except:
                raise ValueError("ROI not understood")


    def map_to_external_reference_aminoacids(self, roi, refname='HXB2', in_patient=True):
        '''Map positions in the patient to a reference, for amino acids in a protein region

        Args:
            roi -- region of interest, a string or a triple (protein, start, stop)
            refname --  reference to compare to
            in_patient -- specifies whether the (start, end) refers to reference or patient coordinates

        returns:
            a (len(roi), 2) array with reference protein coordinates in first column,
                                       patient protein coordinates in second
        '''
        from .filenames import get_coordinate_map_filename, get_custom_reference_filename
        coo_fn = get_coordinate_map_filename(self.name, refname=refname)
        genomewide_map = np.loadtxt(coo_fn, dtype=int)

        if roi in self.annotation:
            region = roi
        else:
            region = roi[0]

        # translate genbank encoded sequence features into a dictionary
        seq = SeqIO.read(get_custom_reference_filename(refname, format='gb'), format='genbank')
        ref_annotation = {x.qualifiers['note'][-1]:x for x in seq.features}

        ref_region_pos = list(ref_annotation[region])
        pat_region_pos = list(self.annotation[region])

        # this is asymmetric because we want the region as defined in the patient
        ind = genomewide_map[np.in1d(genomewide_map[:,1], pat_region_pos)]

        # take protein coordinates
        # NOTE: because the genomewide map may have non-codon gaps, we should check the other two
        # codon positions, but this might give rise to incongruities. It's ok for usual cases
        ind_tmp = []
        for r, p in ind[::3]:
            try:
                r = ref_region_pos.index(r)
                p = pat_region_pos.index(p)
            except ValueError:
                continue
            ind_tmp.append([r, p])
        ind = np.array(ind_tmp, int) // 3

        if roi not in self.annotation:
            start, stop = list(map(int, roi[1:]))
            ind = ind[(ind[:, in_patient] >= start) & (ind[:, in_patient] < stop)]

        return ind


    # TODO: the following is experimental. was meant as a way to easily get an
    # idea what kind of stuff a site is involved in
    def positions_to_features(self, sources=['annotations']):
        '''
        map of positions to features, including the number of proteins, RNA, etc this pos is part of

        Parameters:
           sources: list of sources to take the features from. Allowed sources are
           'annotations', 'shape'.
        '''
        self.pos_to_feature = [{} for nuc in self.reference]

        # Add info from genomic annotations
        if 'annotations' in sources:
            for posfea in self.pos_to_feature:
                posfea['gene'] = 0
                posfea['RNA'] = 0
                posfea['LTR'] = 0
                posfea['codons'] = []
                posfea['protein_codon'] = []

            for fname, feature in self.annotation.items():
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

        # Add info from SHAPE (Siegfried et al. 2014)
        if 'shape' in sources:
            for posfea in self.pos_to_feature:
                posfea['RNA pairing probability'] = None
                posfea['RNA partner'] = None

            from .external import load_pairing_probability_NL43
            pp = load_pairing_probability_NL43()
            m = self.map_to_external_reference('genomewide', 'NL4-3')[:, :2]
            md = dict(m)
            for pos_ref, pos_pat in m:
                if pos_ref in pp.index:
                    #NOTE: more than one pairing is reported
                    tmp = pp.loc[pos_ref]
                    if isinstance(tmp, pd.Series):
                        partner = tmp['partner']
                        prob = tmp['probability']
                    else:
                        tmp = tmp.set_index('partner')['probability']
                        partner = tmp.argmax()
                        prob = tmp.loc[partner]
                    self.pos_to_feature[pos_pat]['RNA pairing probability'] = prob
                    if partner in md:
                        self.pos_to_feature[pos_pat]['RNA partner'] = md[partner]
                    else:
                        self.pos_to_feature[pos_pat]['RNA partner'] = 'missing'

        if 'disorder' in sources:
            from .external import load_disorder_scores_HXB2
            # the following returns a dictionary for each protein that contains a list of positions
            # and values extracted from the tables by Li et al, Positions are codons (hopefully) from hxb2, zero numbering
            dscores = load_disorder_scores_HXB2()
            for prot in dscores:
                # ignore tat and rev since they are split
                if prot in self.annotation and prot not in ['tat', 'rev']:
                    m = self.map_to_external_reference(prot, 'HXB2')
                    # remove trailing stop codon from the protein
                    if prot in ['p6', 'IN', 'gp41', 'vif', 'nef', 'vpu','vpr']:
                        m=m[:-3]
                    hxb2_codons = (m[:,0] - m[0,0])[::3]//3
                    try:
                        # loop over every codon that aligns to HXB2
                        for pos, val in zip(m[0::3,1], dscores[prot]['val'][hxb2_codons]):
                            for ii in range(3): # add the same score to every position in the codon
                                if 'disorder' not in self.pos_to_feature[pos+ii]:
                                    self.pos_to_feature[pos+ii]['disorder']={}
                                self.pos_to_feature[pos+ii]['disorder'][prot] = val
                    except:
                        print("Can't import disorder scores")
                        #import ipdb; ipdb.set_trace()


        if 'accessibility' in sources:
            from .external import load_accessibility
            ascores = load_accessibility()
            for prot in ascores:
                # ignore tat and rev since they are split
                if prot in self.annotation and prot not in ['tat', 'rev']:
                    m = self.map_to_external_reference(prot, 'HXB2')
                    # remove trailing stop codon from the protein
                    if prot in ['p6', 'IN', 'gp41', 'vif', 'nef', 'vpu','vpr']:
                        m=m[:-3]
                    try:
                        # loop over position value pairs
                        for pos, val in ascores[prot]:
                            for ii in range(3): # loop over positions in codon
                                nuc_pos = m[0,0] + pos*3 + ii #nucleotide position in hxb2
                                if nuc_pos in m[:,0]: #if maps to patient
                                    # find index and corresponding position in patient
                                    nuc_ii = np.searchsorted(m[:,0], nuc_pos)
                                    pat_pos = m[nuc_ii,1]
                                    if 'accessibility' not in self.pos_to_feature[pat_pos]:
                                        self.pos_to_feature[pat_pos]['accessibility']={}
                                    self.pos_to_feature[pat_pos]['accessibility'][prot] = val
                    except:
                        import ipdb; ipdb.set_trace()

        if 'structural' in sources:
            from .external import load_structural_effects_NL43
            struct_scores, cons_seqs = load_structural_effects_NL43()
            if struct_scores is not None:
                for prot in struct_scores:
                    # ignore tat and rev since they are split
                    if prot in self.annotation and prot not in ['tat', 'rev']:
                        m = self.map_to_external_reference(prot, 'NL4-3')
                        # remove trailing stop codon from the protein
                        if prot in ['p6', 'IN', 'gp41', 'vif', 'nef', 'vpu','vpr']:
                            m=m[:-3]
                        try:
                            # loop over position value pairs
                            for pi, (pos, val) in enumerate(struct_scores[prot]):
                                for ii in range(3): # loop over positions in codon
                                    nuc_pos = m[0,0] + pos*3 + ii #nucleotide position in hxb2
                                    if prot=='pol': # they start numbering at the start of the PR
                                        nuc_pos+=56*3
                                    if nuc_pos in m[:,0]: #if maps to patient
                                        # find index and corresponding position in patient
                                        nuc_ii = np.searchsorted(m[:,0], nuc_pos)
                                        pat_pos = m[nuc_ii,1]
                                        if 'structural' not in self.pos_to_feature[pat_pos]:
                                            self.pos_to_feature[pat_pos]['structural']={}
                                            self.pos_to_feature[pat_pos]['ref']={}
                                        self.pos_to_feature[pat_pos]['structural'][prot] = val
                                        self.pos_to_feature[pat_pos]['ref'][prot] = cons_seqs[prot][pi]
                        except:
                            import ipdb; ipdb.set_trace()
            else:
                print("couldn load structural scores")


    def get_fragment_depth(self, pad=False, limit_to_dilution = False):
        c = self._annotation_to_fragment_indices('genomewide')
        depth = np.ma.array([s.fragment_depth(c,cov_min=100, var_min=0.05, min_points=10)
                             for s in self.samples])
        if pad:
            for si in range(len(self.samples)):
                depth[si][depth.mask[si]] = self.n_templates_dilutions[si]
                depth.mask[si] = False
        if limit_to_dilution:
            for si in range(len(self.samples)):
                depth[si] = np.minimum(depth[si], self.n_templates_dilutions[si])
        return depth


    def get_hla_type(self, MHC=1):
        '''Get a list with all HLA loci

        Parameters:
           MHC (None/1/2): MHC class I/II only, or all
        '''
        if MHC == 1:
            loci = ('A', 'B', 'C')
        elif MHC == 2:
            loci = ('DRB1', 'DRQ1')
        else:
            loci = ('A', 'B', 'C', 'DRB1', 'DRQ1')

        hla = np.concatenate([[locus+self['HLA-'+locus],
                               locus+self['HLA-'+locus+'-2']]
                              for locus in loci]).tolist()
        return hla


    def get_ctl_epitopes(self,
                         regions=['gag', 'pol',
                                  'gp120', 'gp41',
                                  'vif', 'vpr', 'vpu', 'nef'],
                         kind='mhci=80',
                        ):
        '''Get list of CTL epitopes

        Parameters:
           regions (list): restrict to epitopes within these regions
           kind (str): LANL/epitoolkit/mhci=<n>, where <n> is the cutoff for
           the MHCi predicted list: the first <n> entries are taken.
        '''
        # Get epitope table for patient HLA
        if kind == 'LANL':
            from hivwholeseq.cross_sectional.ctl_epitope_map import (get_ctl_epitope_map,
                                                           get_ctl_epitope_hla)
            ctl_table_main = get_ctl_epitope_map(species='human')
            hla = self.get_hla_type(MHC=1)
            ctl_table_main = get_ctl_epitope_hla(ctl_table_main, hla)
            del ctl_table_main['HXB2 start']
            del ctl_table_main['HXB2 end']
            del ctl_table_main['HXB2 DNA Contig']
            del ctl_table_main['Protein']
            del ctl_table_main['Subprotein']

        elif 'mhci=' in kind:
            n_entries = int(kind[5:])
            from .filenames import get_ctl_epitope_map_filename
            ctl_table_main = pd.read_csv(get_ctl_epitope_map_filename(self.name),
                                         skiprows=3,
                                         sep='\t',
                                         usecols=['peptide'],
                                         # NOTE: top epitopes only, this is a parameter
                                         nrows=n_entries,
                                        )
            ctl_table_main.drop_duplicates(inplace=True)
            ctl_table_main.rename(columns={'peptide': 'Epitope'}, inplace=True)

        else:
            raise ValueError('kind of CTL table not understood')

        data = []
        for region in regions:

            # Restrict epitope table to founder virus sequence
            fea = self.annotation[region]
            regpos = fea.location.nofuzzy_start
            seq = fea.extract(self.reference)
            prot = str(seq.seq.translate())
            ind = [i for i, epi in enumerate(ctl_table_main['Epitope']) if epi in prot]
            ctl_table = ctl_table_main.iloc[ind].copy()

            # Set position in region
            # NOTE: the same epitope could be there twice+ in a protein, so we use
            # regular expressions for that
            import re
            tmp = []
            for epi in ctl_table['Epitope']:
                for match in re.finditer(epi, prot):
                    pos = match.start()
                    tmp.append({'Epitope': epi,
                                'start_region': 3 * pos,
                                'end_region': 3 * (pos + len(epi)),
                               })
            ctl_table = pd.DataFrame(tmp)
            if not len(ctl_table):
                continue

            # Set position genomewide
            ctl_table['start'] = ctl_table['start_region'] + regpos
            ctl_table['end'] = ctl_table['end_region'] + regpos

            # Set start/end positions in HXB2 coordinates
            comap = dict(self.map_to_external_reference(region)[:, ::-2])
            poss = []
            for x in ctl_table['start_region']:
                while True:
                    if x in comap:
                        poss.append(comap[x])
                        break
                    elif x < 0:
                        poss.append(-1)
                        break
                    x -= 1
            ctl_table['start_HXB2'] = np.array(poss, int)
            poss = []
            for x in ctl_table['end_region']:
                while True:
                    if x in comap:
                        poss.append(comap[x])
                        break
                    elif x > 10000:
                        poss.append(-1)
                        break
                    x += 1
            ctl_table['end_HXB2'] = np.array(poss, int)

            # Filter out epitopes for which we cannot find an HXB2 position
            ctl_table = ctl_table.loc[(ctl_table[['start_HXB2', 'end_HXB2']] != -1).all(axis=1)]

            ctl_table['region'] = region

            data.append(ctl_table)

        ctl_table = pd.concat(data).sort('start_HXB2')
        ctl_table.index = list(range(len(ctl_table)))

        return ctl_table


    def get_haplotype_alignment(self, region):
        '''Get alignment of minor haplotypes from whole infection

        Parameters:
           region (str or ROI): if a string, it indicates a genomic region
           (e.g. V3) and the haplotype alignment is read from file. If a
           ROI, i.e. a triple (region, start, end), it extracts the haplotypes
           de novo from the BAM files (TODO).

        TODO: ROI de novo computation not implemented yet.
        '''
        if isinstance(region, str):
            from Bio import AlignIO
            from .filenames import get_haplotype_alignment_filename
            fn = get_haplotype_alignment_filename(self.name, region, 'fasta')
            return AlignIO.read(fn, 'fasta')

        else:
            raise NotImplementedError


    def get_haplotype_alignment_for_insertion(self, position, length=40):
        '''Get alignment of minor haplotypes close to an insertion'''
        def find_closest_alignment_filename():
            import os
            from .filenames import get_haplotype_alignment_filename
            for dist in range(10000):
                for start in [position - dist, position + dist]:
                    region = 'insertion_'+str(start)+'-'+str(start+length)
                    fn = get_haplotype_alignment_filename(self.name, region, 'fasta')
                    if os.path.isfile(fn):
                        return {'start': start,
                                'filename': fn}
            else:
                raise IOError('No alignment found')

        from Bio import AlignIO
        datum = find_closest_alignment_filename()
        datum['ali'] = AlignIO.read(datum['filename'], 'fasta')

        return datum
