# vim: fdm=indent
'''
author:     Fabio Zanini
date:       20/05/15
content:    Support module for nucleid and amino acid sequences.
'''
# Modules
from numpy import array



# Globals
# Alphabet of nucleotides
alphas = 'ACGT-N'
alphal = list(alphas)
alpha = array(alphal, 'S1')

# Alphabet of amino acids
alphaas = 'ACDEFGHIKLMNPQRSTVWY*-X'
alphaal = list(alphaas)
alphaa = array(alphaal, 'S1')
