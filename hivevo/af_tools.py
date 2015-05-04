# vim: fdm=marker
'''
author:     Fabio Zanini/Richard Neher
date:       25/04/2015
content:    Data access module HIV patients.
'''
# Modules
import numpy as np
from hivwholeseq.utils.sequence import alpha, alphaa

def diversity(af):
    return 1 - np.mean(np.sum(af**2, axis=0))

def divergence(af, initial):
    return 1 - np.mean(af[initial,np.arange(len(initial))])

def majority_frequency(af):
    return np.max(af, axis=0)

def majority_indices(af):
    return np.argmax(af, axis=0)

def consensus(af):
    return alpha[majority_indices(af)]
