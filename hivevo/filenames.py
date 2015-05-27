# vim: fdm=indent
'''
author:     Fabio Zanini
date:       20/05/15
content:    Support module for file names and paths.
'''
# Modules
import os

# FIXME: using existing factory, copy/paste this (but it reveals PATIENT DATA!)
from hivwholeseq.patients.filenames import *



# Globals
if os.path.isdir('/media/FZ_MPI/HIV_Sweden/'):
    root_data_folder = '/media/FZ_MPI/HIV_Sweden/'
elif os.path.isdir('/var/www/hivwholeweb/'):
    root_data_folder = '/var/www/hivwholeweb/app/hiv/static/data/'
elif os.path.isdir('/home/fabio/') and (not os.path.isdir('/ebio/ag-neher/share/data')):
    root_data_folder = '/home/fabio/university/phd/sequencing/data/'
else:
    root_data_folder = '/ebio/ag-neher/share/data/MiSeq_HIV_Karolinska/'

reference_folder = root_data_folder+'reference/'

self = os.path.dirname(os.path.abspath(__name__.replace('.', '/')))+'/'
table_folder = self + '../data/tables/'



# Functions
def get_table_filename(kind, format='xlsx'):
    return table_folder+kind+'.'+format


def get_custom_reference_filename(reference, format='fasta'):
    '''Get the filename of a custom reference sequence'''
    filename = reference
    filename = filename+'.'+format
    return reference_folder+filename


def get_subtype_reference_alignment_filename(subtype='B', format='fasta', refname='HXB2'):
    '''Get the filename of a reference alignment'''
    filename = 'genomewide.'+subtype+'.nuc.aligned.'+format
    return reference_folder+'reference/alignments/pairwise_to_'+refname+'/'+filename
