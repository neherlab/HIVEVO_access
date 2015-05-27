# vim: fdm=indent
'''
author:     Fabio Zanini
date:       20/05/15
content:    Support module for file names and paths.
'''
# FIXME: using existing factory, copy/paste this (but it reveals PATIENT DATA!)
from hivwholeseq.patients.filenames import *
from hivwholeseq.filenames import table_filename
from hivwholeseq.filenames import get_custom_reference_filename

def get_subtype_alignment_filename(subtype='B'):
	if subtype=='B':
		return '/ebio/ag-neher/share/data/MiSeq_HIV_Karolinska/reference/alignments/pairwise_to_HXB2/genomewide.B.nuc.aligned.fasta'
	else:
		print('subtype not availabel')
		return None
