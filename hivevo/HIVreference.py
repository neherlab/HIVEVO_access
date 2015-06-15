# vim: fdm=marker
'''
author:     Fabio Zanini/Richard Neher
date:       25/04/2015
content:    Data access module HIV patients.
'''
# Modules
import numpy as np
import pandas as pd
from Bio import SeqIO, AlignIO
from .sequence import alpha
from .filenames import get_custom_reference_filename


class HIVreference(object):
    """docstring for HIVreference"""
    def __init__(self, refname='HXB2', subtype='B', load_alignment=True):
        self.refname = refname
        self.seq = SeqIO.read(get_custom_reference_filename(self.refname, format = 'gb'), format='genbank')
        # translate genbank encoded sequence features into a dictionary
        self.annotation = {x.qualifiers['note'][-1]:x for x in self.seq.features}

        if load_alignment:
            from .filenames import get_subtype_reference_alignment_filename
            self.aln = np.array(AlignIO.read(get_subtype_reference_alignment_filename(subtype=subtype), 'fasta'))
            self.calc_nucleotide_frequencies()

        else:
            from .filenames import get_subtype_reference_allele_frequencies_filename
            self.af = np.load(get_subtype_reference_allele_frequencies_filename(subtype=subtype))

        self.consensus_indices = np.argmax(self.af, axis=0)
        self.consensus = alpha[self.consensus_indices]
        self.calc_entropy()


    def calc_nucleotide_frequencies(self):
        self.af = np.zeros((len(alpha)-1, self.aln.shape[1]), dtype = float)
        for ni, nuc in enumerate(alpha[:-1]):
            self.af[ni,:] = np.sum(self.aln==nuc, axis=0)
        cov = np.sum(self.af, axis=0)
        self.af/=cov


    def calc_entropy(self):
        self.entropy = np.maximum(0,-np.sum(self.af*np.log(1e-10+self.af), axis=0))


    def map_to_sequence_collection():
        pass


    def get_ungapped(self, threshold = 0.05):
        return self.af[-1,:]<0.05


    def get_entropy_quantiles(self, q):
        from scipy.stats import scoreatpercentile
        thresholds = [scoreatpercentile(self.entropy, 100.0*i/q) for i in range(q+1)]
        return {i: {'range':(thresholds[i],thresholds[i+1]), 
                    'ind':np.where((self.entropy>=thresholds[i])*(self.entropy<thresholds[i+1]))[0]}
               for i in range(q)}


    def get_entropy_in_patient_region(self, map_to_ref):
        '''
        returns entropy in a specific regions defined by a set of indices in the reference
        params:
        map_to_ref  --  either a one dimensional vector specifying indices in the reference
                        or a (3, len(region)) array with the reference coordinates in the first column
                        this is the output of Patient.map_to_external_reference
        '''
        if len(map_to_ref.shape)==2:
            return self.entropy[map_to_ref[:,0]]
        elif len(map_to_ref.shape)==1:
            return self.entropy[map_to_ref]


    def get_consensus_in_patient_region(self, map_to_ref):
        '''
        returns consensus in a specific regions defined by a set of indices in the reference
        params:
        map_to_ref  --  either a one dimensional vector specifying indices in the reference
                        or a (3, len(region)) array with the reference coordinates in the first column
                        this is the output of Patient.map_to_external_reference
        '''
        if len(map_to_ref.shape)==2:
            return self.consensus[map_to_ref[:,0]]
        elif len(map_to_ref.shape)==1:
            return self.consensus[map_to_ref]


    def get_consensus_indices_in_patient_region(self, map_to_ref):
        '''
        returns consensus_indices in a specific regions defined by a set of indices in the reference
        params:
        map_to_ref  --  either a one dimensional vector specifying indices in the reference
                        or a (3, len(region)) array with the reference coordinates in the first column
                        this is the output of Patient.map_to_external_reference
        '''
        if len(map_to_ref.shape)==2:
            return self.consensus_indices[map_to_ref[:,0]]
        elif len(map_to_ref.shape)==1:
            return self.consensus_indices[map_to_ref]
