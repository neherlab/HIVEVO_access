# vim: fdm=marker
'''
author:     Fabio Zanini
date:       07/10/14
content:    Description module for HIV patient samples.
'''
# Modules
import numpy as np
import pandas as pd
import sys,os
sys.path.append('/ebio/ag-neher/share/users/rneher/mapping')
from hivwholeseq.sequencing.filenames import table_filename



# Classes
class SamplePat(pd.Series):
    '''Patient sample'''

    def __init__(self, *args, **kwargs):
        '''Initialize a patient sample'''
        super(SamplePat, self).__init__(*args, **kwargs)
        self._sequenced_samples = None

    @property
    def _constructor(self):
        return SamplePat


    def get_n_templates_dilutions(self):
        '''Get the time course of the number of templates to PCR, limiting depth'''
        from hivwholeseq.patients.get_template_number import get_template_number
        return get_template_number(self.dilutions)

    def _get_allele_counts(self, fragment, use_PCR1, VERBOSE):
        '''Get allele counts for a single patient sample'''
        if VERBOSE >= 1:
            print 'Getting allele counts:', self.patient, self.name, fragment
    
        # FIXME: add depth min
        from hivwholeseq.patients.filenames import get_allele_counts_filename
        # PCR1 filter here
        fn1 = get_allele_counts_filename(self.patient, self.name, fragment, PCR=1)
        fn2 = get_allele_counts_filename(self.patient, self.name, fragment, PCR=2)
        if os.path.isfile(fn1):
            fname = fn1
            if VERBOSE >= 3:
                print self.name, 1
        elif os.path.isfile(fn2) and use_PCR1<2:
            fname = fn2
            if VERBOSE >= 3:
                print self.name, 2
        else:
            return None, None
        ac = np.load(fname).sum(axis=0)
        return ac

    def get_allele_counts(self, coordinates, add = True, cov_min = 100, use_PCR1=1, VERBOSE=0, **kwargs):
        ac = np.ma.zeros((6, coordinates['length']), dtype = int)
        for fragment in ['F'+str(i) for i in range(1,7)]:
            if fragment in coordinates:
                tmp_ac = self._get_allele_counts(fragment, use_PCR1, VERBOSE)
                if tmp_ac is not None:
                    if add:
                        ac[:,coordinates[fragment][0]]+=tmp_ac[:,coordinates[fragment][1]]
                    else:
                        ac[:,coordinates[fragment][0]]=np.maximum(ac[:,coordinates[fragment][0]], 
                                                                  tmp_ac[:,coordinates[fragment][1]])
        if cov_min is not None:
            cov = ac.sum(axis=0)
            ac.mask = np.repeat([cov<cov_min], ac.shape[0], axis=0)
        return ac

    def get_allele_frequencies(self, coordinates, add=True, cov_min = 100, use_PCR1=1, VERBOSE=0, **kwargs):
        ac = self.get_allele_counts(coordinates, add=add, cov_min=cov_min, use_PCR1=use_PCR1, VERBOSE=VERBOSE, **kwargs)
        cov = ac.sum(axis=0)
        af = np.ma.masked_array(ac, dtype=float)
        af/=cov+1e-10
        af.mask += np.repeat([cov_min<0], aft.shape[0], axis=0)
        return af
        
    def get_coverage(self, coordinates, add=True,use_PCR1=1, VERBOSE=0, **kwargs):
        ac = self.get_allele_counts(coordinates, add=add, cov_min=None, use_PCR1=use_PCR1, VERBOSE=VERBOSE, **kwargs)
        return ac.sum(axis=0)

    def get_local_haplotypes(self, fragment, start, stop):
                             VERBOSE=0, maxreads=-1, filters=None, PCR=1):

        
# Functions
def load_samples_sequenced(patients=None, include_wrong=False, include_cell=False):
    '''Load patient samples sequenced from general table'''
    sample_table = pd.read_excel(table_filename, 'Samples timeline sequenced',
                                 index_col=0)

    sample_table.index = pd.Index(map(str, sample_table.index))
    sample_table.loc[:, 'patient'] = map(str, sample_table.loc[:, 'patient'])
    # Note: this refers to the TOTAL # of templates, i.e. the factor 2x for
    # the two parallel RT-PCR reactions
    sample_table['n templates'] = sample_table['viral load'] * 0.4 / 12 * 2

    if not include_wrong:
        sample_table = sample_table.loc[sample_table.loc[:, 'wrong'] != 'x']
        del sample_table['wrong']

    if not include_cell:
        sample_table = sample_table.loc[sample_table.loc[:, 'sample type'] == 'RNA']

    if patients is not None:
        sample_table = sample_table.loc[sample_table.loc[:, 'patient'].isin(patients)]

    return [SamplePat(val) for x, val in sample_table.iterrows()]

