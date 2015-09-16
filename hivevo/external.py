# vim: fdm=indent
'''
author:     Fabio Zanini
date:       16/09/15
content:    Module to interface with external data.
'''
# Modules
from pandas import read_csv

from .filenames import local_data_folder


def load_pairing_probability_NL43():
    data_filename = local_data_folder+'external/nmeth.3029-S4.txt'
    data = read_csv(data_filename, header=1, sep='\t', index_col=0)
    data.rename(columns={'j': 'partner'}, inplace=True)
    data.index.name = 'position'

    # Use positive pairing probabilities
    data.rename(columns={'-log10(Probability)': 'probability'}, inplace=True)
    data['probability'] = 10**(-data['probability'])

    # Indices start from 0 in Python
    data.index -= 1
    data['partner'] -= 1

    # The sequence of Siegfried et al. 2014 starts later than standard NL4-3
    start = 454
    data.index += start
    data['partner'] += start

    return data

