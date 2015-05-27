# vim: fdm=marker
'''
author:     Fabio Zanini & Richard Neher
date:       01/05/15
content:    Description module for HIV patient samples.
'''
# Modules
import numpy as np
import pandas as pd
import os
import sys
from collections import defaultdict
from itertools import izip

all_fragments = ['F'+str(i) for i in range(1,7)]


# Classes
class Sample(pd.Series):
    '''
    Class providing access to sample specific quantities such as as depth, counts etc
    '''

    def __init__(self, *args, **kwargs):
        '''Initialize a patient sample'''
        super(Sample, self).__init__(*args, **kwargs)
        self._sequenced_samples = None

    @property
    def _constructor(self):
        return Sample

    def get_n_templates_dilutions(self):
        '''Get the number of templates to PCR calculated by interpolating the limiting dilution'''

        dilutions = np.array([1, 10,100,1000, 10000, 100000], dtype = int)
        successful = []
        dilstr = self.dilutions
        # parse the dilutions string expected to be of type 1:100 (1/2), where (1/2) mean one successful
        # amplification out of 2 trials
        if not isinstance(dilstr, basestring):
            print "expecting a string, got", dilstr, "falling back on viral load /60:", self['viral load']/60.0
            val = self['viral load']/60.0 
        else:
            try:
                # extract the last positive dilution (dil_factor) and the number of successful dilutions 
                dil_factor = float(dilstr.split()[0].split(':')[1])
                n_positive = map(int, dilstr.split()[1][1:-1].split('/'))
                for dil in dilutions:
                    if dil<dil_factor: #assume all previous dilutions were successful
                        successful.append([2,2])
                    elif dil==dil_factor:
                        successful.append(n_positive)
                    else:
                        successful.append([0,2])  # assume all subsequent dilutions failed
                successful = np.array(successful) # -> successful now contains [:,0] the number of successful dilutions in [:, 1] trials
        
                def prob(logcn, dil, suc):
                    # no amplifiable molecule: p=exp(-cn/dil*0.5) -> 2/2 (1-p)^2; 1/2: 2p(1-p); 0/2: p^2
                    # the extra 0.5 in the exponent accounts for the fact that we do to series starting with
                    # a total of 1/10th of the material for fragment 4
                    p = np.exp(-np.exp(logcn)/dil*0.5)
                    return -np.sum( np.log((suc[:,0]==2)*(1-p)**2 + (suc[:,0]==1)*2*p*(1-p) + (suc[:,0]==0)*p**2 ))
        
                # numerically optimize the log probability of observing the pattern of dilutions
                from scipy.optimize import minimize
                x_opt = minimize(prob, x0=2.0, args = (dilutions, successful), method='Powell')
                val = np.exp(x_opt.x)
            except:
                print "Parsing", dilstr, "didn't work. falling back on viral load /60:", self['viral load']/60.0
                val = self['viral load']/60.0 
        return val

    def _get_fragment_counts(self, fragment, use_PCR1, VERBOSE):
        '''
        Get allele counts for a specific fragment of a specified PCR
        this method is not meant for general use
        '''
        if VERBOSE >= 1:
            print 'Getting allele counts:', self.patient, self.name, fragment

        from .filenames import get_allele_counts_filename
    
        # FIXME: add depth min
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
            return None
        ac = np.load(fname).sum(axis=0)
        return ac

    def get_allele_counts(self, coordinates, add = True, cov_min = 100, use_PCR1=1, VERBOSE=0, **kwargs):
        '''
        get counts at positions specified by coordinates
        parameters:
        coordinates  -- a dictionary with  fields 'F1': (array(pos in region of interest), array(pos on fragment))
                        in addition, a field 'length': total length of the region of interest is required
        add          -- if True, add counts when fragments overlap (default), other wise take max
        use_PCR1     -- 1 == use PCR1 is available (default), 2 == force PCR1 or fail
        cov_min      -- mask values where the final coverage is below that threshold
        '''
        ac = np.ma.zeros((6, coordinates['length']), dtype = int)
        for fragment in all_fragments:
            if fragment in coordinates:
                tmp_ac = self._get_fragment_counts(fragment, use_PCR1, VERBOSE)
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
        '''
        get counts at positions specified by coordinates -- all options as get_allele_counts
        '''
        ac = self.get_allele_counts(coordinates, add=add, cov_min=cov_min, use_PCR1=use_PCR1, VERBOSE=VERBOSE, **kwargs)
        cov = ac.sum(axis=0)
        af = np.ma.masked_array(ac, dtype=float)
        af/=cov+1e-10
        af.mask += np.repeat([cov<cov_min], af.shape[0], axis=0)
        return af
        
    def get_coverage(self, coordinates, add=True,use_PCR1=1, VERBOSE=0, **kwargs):
        ac = self.get_allele_counts(coordinates, add=add, cov_min=None, use_PCR1=use_PCR1, VERBOSE=VERBOSE, **kwargs)
        return ac.sum(axis=0)

    def get_cocounts(self, fragment, use_PCR1=True, compressed=True):
        '''
        get joint counts of nuc1 at pos1 and nuc2 at pos2 for a fragment
        returns:
        ac     --  array of dimension (6,6,L,L) where L is the length of the fragment
        '''
        from hivwholeseq.patients.filenames import get_allele_cocounts_filename
        try: # TODO: compressed or not needs cleaning up
            fname = get_allele_cocounts_filename(self.patient, self.name, fragment, PCR=1 if use_PCR1 else 2, compressed=compressed)
            ac = np.load(fname)['cocounts']
        except:
            fname = get_allele_cocounts_filename(self.patient, self.name, fragment, PCR=1 if use_PCR1 else 2, compressed=False)
            ac = np.load(fname)

        return ac

    def get_pair_frequencies(self, fragment, var_min=0, use_PCR1=True, compressed=True):
        '''
        return the fraction of observations of nuc1, nuc2 and pos1, pos2 reduced to variable positions
        parameters:
        fragment -- the fragments to obtain the frequencies of
        var_min  -- the minimal on site variability to be included in the matrix. 1-\sum_i x_i**2 >var_min
        returns:
        pos      -- variable positions
        af2p     -- pair frequencies at variable positions
        cov      -- coverage at the pairs
        af1p     -- single nucleotide frequencies
        '''
        import gc
        try:
            acc = self.get_cocounts(fragment, use_PCR1=use_PCR1, compressed=compressed)
        except:
            return_args = (None, None, None, None)
        else:
            # calculate single nucleotide frequencies and determine variable sites
            af1p = np.array(acc[:,:,np.arange(acc.shape[2]), np.arange(acc.shape[3])].sum(axis=1),dtype=float)
            af1p = af1p/(1e-10+af1p.sum(axis=0))
            variable_sites = np.sum(af1p, axis=0)-np.sum(af1p**2, axis=0)>var_min
            if variable_sites.sum()>1: # check whether there are at least 2 variable sites
                reduced_af1p = af1p[:,variable_sites]
                positions = np.where(variable_sites)[0]
                n_variable_sites = reduced_af1p.shape[-1]
                # make a reduced matric
                reduced_af2p = np.zeros((acc.shape[0], acc.shape[1], n_variable_sites, n_variable_sites), dtype = float)
                for di,dsite in enumerate(positions):
                    reduced_af2p[:,:,di,:] = acc[:,:,dsite,variable_sites]

                reduced_acc_cov = np.array(reduced_af2p.sum(axis=1).sum(axis=0), dtype=int)
                reduced_af2p /= (1e-10+reduced_acc_cov) 
                return_args = (positions, reduced_af2p, reduced_acc_cov, af1p[:,variable_sites])
            else:
                print "no variable sites"
                return_args = (None, None, None, None)
            del acc, af1p
            gc.collect()
        return return_args

    def fragment_depth(self, coordinates, var_min = 0.03, cov_min = 100, min_points = 5, pseudo_counts = 3):
        '''
        returns estimates of the fragment specific depth estimated from overlap frequencies
        parameters:
        coordinates   --  coordinates of all fragments as obtained for allele counts of the entire genome
        var_min       --  minimal variance to include into the estimation
        cov_min       --  minimal coverage to include into the estimation
        min_points    --  the minimal number of valid frequencies required to attempt estimation
        pseudo_counts --  number of additional pseudo_counts to include, these are equal to the depth estimate from dilutions
        '''
        overlaps = zip(all_fragments[:-1], all_fragments[1:])
        # extract variant frequencies in overlaps
        overlap_frequencies = []
        for f1,f2 in overlaps:
            overlap_frequencies.append([])
            for f,fother in [(f1,f2), (f2,f1)]:
                tmp_ind = np.in1d(coordinates[f][0], coordinates[fother][0])
                frag_ind = coordinates[f][1][tmp_ind]
                af = self.get_allele_frequencies({f:(coordinates[f][1],coordinates[f][1]), 
                                                 'length':len(coordinates[f][0])}, cov_min=100)[:,frag_ind]
                overlap_frequencies[-1].append(af)

        # subset frequencies in overlaps to variable positions
        for junction in overlap_frequencies:
            mean_freq, delta_freq = 0.5*(junction[0]+junction[1]), junction[0]-junction[1]
            tmp_variable = (mean_freq>var_min)*(mean_freq<1-var_min)
            good_positions = (tmp_variable * (mean_freq.mask==False)).nonzero()
            if len(good_positions[0])>min_points:
                junction[0] = np.array(mean_freq[good_positions])
                junction[1] = np.array(delta_freq[good_positions])
            else:
                junction[0], junction[1] = None, None

        # calculate overlap variances
        neff_overlaps = [1.0/np.mean(np.concatenate([0.5*y**2/(x*(1-x)), 
                                     np.ones(pseudo_counts, dtype=float)/self.get_n_templates_dilutions()]))
                         if x is not None else None for x,y in overlap_frequencies]

        # split into chains
        chains = [[]]
        for (f1,f2), neff in izip(overlaps, neff_overlaps):
            if neff is None:
                chains.append([])
            else:
                chains[-1].append((f1,f2,neff))
        chains = filter(lambda x:len(x), chains) # clean empty chains out

        # distribute variation to fragments
        neff_fragments = np.ma.masked_all(len(all_fragments))
        for chain in chains:
            for f1,f2,neff in chain:
                neff_fragments[all_fragments.index(f1)] = max(neff,neff_fragments[all_fragments.index(f1)])
                neff_fragments[all_fragments.index(f2)] = max(neff,neff_fragments[all_fragments.index(f2)])

        return neff_fragments

    # TODO: the following doesn't work 
    def haplotypes(self, fragment, start, stop, VERBOSE=0, maxreads=-1, filters=None, PCR=1):
        from hivwholeseq.patients.get_local_haplotypes import get_local_haplotypes
        from .filenames import get_mapped_filtered_filename

        bam_fname = get_mapped_filtered_filename(self.patient, self.name, fragment, type='bam', PCR=PCR, decontaminated=True)
        try:
            return get_local_haplotypes(bam_fname, start, stop)
        except:
            raise ValueError("can't read haplotypes")

    def tile_region(self, coordinates, length, padding):
        haps = defaultdict(list)
        for frag in all_fragments:
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
def load_samples_sequenced(patients=None, include_empty=False):
    '''Load patient samples sequenced from general table'''
    from .filenames import get_table_filename
    sample_table = pd.read_excel(get_table_filename('samples'),
                                 'Samples',
                                 index_col=0)

    # Reindex DataFrame
    sample_table.index = pd.Index(map(str, sample_table.index))
    sample_table.loc[:, 'patient'] = map(str, sample_table.loc[:, 'patient'])

    # Note: this refers to the TOTAL # of templates, i.e. the factor 2x for
    # the two parallel RT-PCR reactions
    sample_table['n templates'] = sample_table['viral load'] * 0.4 / 12 * 2

    if not include_empty:
        ind = [i for i, sample in sample_table.iterrows()
               if (sample[['F1', 'F2', 'F3', 'F4', 'F5', 'F6']] != 'miss').any()]
        sample_table = sample_table.loc[ind]

    if patients is not None:
        sample_table = sample_table.loc[sample_table.loc[:, 'patient'].isin(patients)]

    return [Sample(val) for x, val in sample_table.iterrows()]
