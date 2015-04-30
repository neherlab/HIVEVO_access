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


class HIVreference(object):
	"""docstring for HIVreference"""
	def __init__(self, ref_name='HXB2'):
		super(HIVreference, self).__init__()
		self.ref_name = ref_name
		self.seq = SeqIO.read()

	def map_to_sequence_collection():
		pass

	def entropy(self):
		pass

	def consensus(self):
		pass
		

		