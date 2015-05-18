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
    return np.mean(af.sum(axis=0)) - np.mean(np.sum(af**2, axis=0))

def divergence(af, initial):
    return np.mean(af.sum(axis=0)) - np.mean(af[initial,np.arange(len(initial))])

def majority_frequency(af):
    return np.max(af, axis=0)

def majority_indices(af):
    return np.argmax(af, axis=0)

def consensus(af):
    return alpha[majority_indices(af)]

def LD(af2p, af1p, cov, cov_min = 100):
    p = af1p.max(axis=0)
    pi = af1p.argmax(axis=0)
    q=1-p

    p12 = np.ma.zeros(cov.shape, dtype = float)
    ind1 = np.arange(pi.shape[0])
    ind2 = np.ones(pi.shape[0], dtype=int)   
    #import pdb; pdb.set_trace()
    for ii, nuci in enumerate(pi):
        p12[ii,:] = af2p[nuci][(pi,ii*ind2,ind1)]
    p12.mask = cov<cov_min
    np.fill_diagonal(p12.mask, True)
    p1p2 = np.outer(p,p)
    p1q1p2q2 = np.outer(p*q,p*q)
    p1q2 = np.outer(p,q)

    D = p12 - p1p2
    LD = np.sqrt(D**2/(1e-10+p1q1p2q2))

    Dp = D
    Dp[Dp>0] /= np.minimum(p1q2, p1q2.T)[Dp>0]
    Dp[Dp<0] /= np.minimum(p1p2, p1p2.T)[Dp<0]
    Dp = np.abs(Dp)
    np.fill_diagonal(LD,0)
    np.fill_diagonal(Dp,0)

    return LD, Dp, p12

