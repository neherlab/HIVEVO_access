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
from .filenames import get_custom_reference_filename


class HIVreference(object):
    """docstring for HIVreference"""
    def __init__(self, ref_name='HXB2'):
        self.ref_name = ref_name
        self.seq = SeqIO.read(get_custom_reference_filename('HXB2', format = 'gb'), format='genbank')
        # translate genbank encoded sequence features into a dictionary
        self.annotation = {x.qualifiers['note'][-1]:x for x in self.seq.features}

    def map_to_sequence_collection():
        pass

    @property   
    def entropy(self):
        pass

    def consensus(self):
        pass
        
    def get_entropy_quantiles(self, q):
        from scipy.stats import score_at_percentile
        thresholds = [score_at_percentile(self.entropy, 100.0*i/q) for i in range(q+1)]
        return {i: {'range':(thresholds[i],thresholds[i+1]), 
                    'ind':np.where((self.entropy>=thresholds[i])*(self.entropy<thresholds[i+1]))[0]}
               for i in range(q)}
        