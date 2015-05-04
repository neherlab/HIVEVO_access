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
from collections import defaultdict


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
        dilutions = np.array([1, 10,100,1000, 10000, 100000], dtype = int)
        successful = []
        dilstr = self.dilutions
        # parse the dilutions string expected to be of type 1:100 (1/2), where (1/2) mean one successful
        # amplification out of 2 trials
        if not isinstance(dilstr, basestring):
            print "expecting a string"
        else:
            dil_factor = float(dilstr.split()[0].split(':')[1])
            n_positive = map(int, dilstr.split()[1][1:-1].split('/'))
        for dil in dilutions:
            if dil<dil_factor:
                successful.append([2,2])
            elif dil==dil_factor:
                successful.append(n_positive)
            else:
                successful.append([0,2])
        successful = np.array(successful) # -> successful now contains [0] the number of successful dilutations in [1] trials

        def prob(logcn, dil, suc):
            # no amplifiable molecule: p=exp(-cn/dil*0.5) -> 2/2 (1-p)^2; 1/2: 2p(1-p); 0/2: p^2
            # the extra 0.5 in the exponent accounts for the fact that we do to series starting with
            # a total of 1/10th of the material for fragment 4
            p = np.exp(-np.exp(logcn)/dil*0.5)
            return -np.sum( np.log((suc[:,0]==2)*(1-p)**2 + (suc[:,0]==1)*2*p*(1-p) + (suc[:,0]==0)*p**2 ))

        from scipy.optimize import minimize
        x_opt = minimize(prob, x0=2.0, args = (dilutions, successful), method='Powell')
        val = np.exp(x_opt.x)
        return val

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
                if tmp_ac[0] is not None:
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
        af.mask += np.repeat([cov<cov_min], af.shape[0], axis=0)
        return af
        
    def get_coverage(self, coordinates, add=True,use_PCR1=1, VERBOSE=0, **kwargs):
        ac = self.get_allele_counts(coordinates, add=add, cov_min=None, use_PCR1=use_PCR1, VERBOSE=VERBOSE, **kwargs)
        return ac.sum(axis=0)

    def haplotypes(self, fragment, start, stop, VERBOSE=0, maxreads=-1, filters=None, PCR=1):
        from hivwholeseq.patients.get_local_haplotypes import get_local_haplotypes
        from hivwholeseq.patients.filenames import get_mapped_filtered_filename
        bam_fname = get_mapped_filtered_filename(self.patient, self.name, fragment, type='bam', PCR=PCR, decontaminated=True)
        #try:
        return get_local_haplotypes(bam_fname, start, stop)
    #    except:
     #       raise ValueError("can't read haplotypes")

    def tile_region(self, coordinates, length, padding):
        haps = defaultdict(list)
        for frag in ['F'+str(i) for i in range(1,7)]:
            if frag in coordinates:
                region_indices = coordinates[frag][0]
                fragment_indices = coordinates[frag][1]
                global_indices = coordinates[frag][0] + coordinates['start']
                start = fragment_indices[0]
                stop = start + length + padding
                while stop < fragment_indices[-1]:
                    if True: #try:
                        res = (start, stop, self.haplotypes(frag, start, stop))
                    else: #except:
                        print "can't read haplotypes"
                        res = (start, stop, None)
                    haps[frag].append(res)
                    start += length
                    stop = start + length + padding
        return haps

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

