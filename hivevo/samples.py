# vim: fdm=marker
'''
author:     Fabio Zanini & Richard Neher
date:       01/05/15
content:    Description module for HIV patient samples.
'''
# Modules
from collections import Counter
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
            print (self.name+": Templates from dilutions expecting a string, got", dilstr,
                   "falling back on viral load /60:", self['viral load']/60.0)
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


    def get_allele_counts(self, coordinates, add=True, cov_min=100, VERBOSE=0,
                          type='nuc',
                          **kwargs):
        '''
        get counts at positions specified by coordinates
        parameters:
        coordinates  -- a dictionary with different contents for nucleotides and amino acids.
          - nucleotides: fields 'F1': (array(pos in region of interest), array(pos on fragment))
                         a field 'length': total length of the region of interest is required
          - amino acids: fields 'PR': array(pos in region of interest)
                         if the array is None, all positions are taken
        add          -- if True, add counts when fragments overlap (default), other wise take max
                        this parameter is only used for nucleotides
        cov_min      -- mask values where the final coverage is below that threshold
        type         -- 'nuc' for nucleotides, 'aa' for amino acids, 'cod' for codons

        NOTE: the difference in the coordinate data structure for nucleotides and
        amino acids stems from the different storage strategies. Nucleotides are
        saved fragment by fragment to allow a separation of the PCR of origin.
        Amino acid may merge fragments because they need a frame of reference
        and it is confusing to separate by fragment AND protein at the same time.
        If you do need such surgical tayloring, you can recompute the amino acid
        counts from the reads (which are split by fragment).

        Examples:

        1. Nucleotides from a certain fragment

        sample.get_allele_counts({'F3': [[0, 1, 2], [56, 57, 58]], 'length': 100},
                                 type='nuc')


        2. Amino acids from a certain protein

        sample.get_allele_counts({'PR': [0, 1, 2]}, type='aa')
        '''
        from .filenames import get_allele_counts_filename

        if type == 'nuc':
            from .sequence import alpha
            ac = np.ma.zeros((len(alpha), coordinates['length']), dtype=int)
            for fragment, coord in coordinates.iteritems():
                fname = get_allele_counts_filename(self.name, fragment, type=type)
                if not os.path.isfile(fname):
                    continue
                tmp_ac = np.load(fname).sum(axis=0)
                try:
                    if add:
                        ac[:,coord[0]] += tmp_ac[:, coordinates[fragment][1]]
                    else:
                        ac[:,coord[0]] = np.maximum(ac[:, coord[0]], tmp_ac[:, coord[1]])
                except:
                    print("can't load allele counts:", fname)

        elif type == 'aa':
            # The dict should have a single protein not the best data structure
            for protein, coord in coordinates.iteritems():
                fname = get_allele_counts_filename(self.name, protein, type=type)
                ac = np.load(fname)
                if coord is not None:
                    ac = ac[:, coord]
                break

        else:
            raise ValueError('Data type not understood')

        ac = np.ma.asarray(ac)

        if cov_min is not None:
            cov = ac.sum(axis=0)
            ac.mask = np.repeat([cov<cov_min], ac.shape[0], axis=0)
        return ac


    def get_insertions(self, coordinates, add=True, VERBOSE=0, **kwargs):
        '''Get insertions from this sample'''
        from .filenames import get_insertions_filename
        ic = Counter()
        for fragment, coord in coordinates.iteritems():
            fname = get_insertions_filename(self.name, fragment)
            if not os.path.isfile(fname):
                continue
            coordd = dict(np.array(coord).T[:, ::-1])
            try:
                tmp_ic = pd.read_pickle(fname)
            except IOError:
                continue
            for (position, insertion), value in tmp_ic.iteritems():
                if position not in coordd:
                    continue
                key = (coordd[position], insertion)
                if add:
                    ic[key] += value
                else:
                    ic[key] = max(ic[key], value)
        ic = pd.Series(ic, name='insertions')
        if len(ic):
            ic.index.names = ['position', 'insertion']
        return ic


    def get_allele_frequencies(self, coordinates, add=True, cov_min = 100, VERBOSE=0, type='nuc', **kwargs):
        '''
        get counts at positions specified by coordinates -- all options as get_allele_counts
        '''
        ac = self.get_allele_counts(coordinates, add=add, cov_min=cov_min, VERBOSE=VERBOSE, type=type, **kwargs)
        cov = ac.sum(axis=0)
        af = np.ma.masked_array(ac, dtype=float)
        af/=cov+1e-10
        af.mask += np.repeat([cov<cov_min], af.shape[0], axis=0)
        return af


    def get_coverage(self, coordinates, add=True, VERBOSE=0, **kwargs):
        ac = self.get_allele_counts(coordinates, add=add, cov_min=None, VERBOSE=VERBOSE, **kwargs)
        return ac.sum(axis=0)


    def get_cocounts(self, fragment, compressed=True, type='nuc'):
        '''
        get joint counts of nuc1 at pos1 and nuc2 at pos2 for a fragment
        returns:
        ac     --  array of dimension (6,6,L,L) where L is the length of the fragment
        '''
        # FIXME: get_allele_counts takes a "coordinates" object, whereas cocounts
        # takes a string... this should be homogenized
        from .filenames import get_allele_cocounts_filename
        if compressed:
            fname = get_allele_cocounts_filename(self.name, fragment, type=type, format='npz')
            ac = np.load(fname)['cocounts']
        else:
            fname = get_allele_cocounts_filename(self.name, fragment, type=type, format='npy')
            ac = np.load(fname)

        return ac


    def get_pair_frequencies(self, fragment, var_min=0, compressed=True):
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
        # NOTE: because the files are very big, we manually control garbage collection
        import gc

        try:
            acc = self.get_cocounts(fragment, compressed=compressed)

        # FIXME: all-catching is a bad practice that creates bugs, need to change
        # the following line
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

        finally:
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
    def haplotypes(self, fragment, start, stop, VERBOSE=0, maxreads=-1, filters=None):
        from hivwholeseq.patients.get_local_haplotypes import get_local_haplotypes
        from .filenames import get_mapped_filtered_filename

        bam_fname = get_mapped_filtered_filename(self.patient, self.name, fragment, type='bam', decontaminated=True)
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
def load_samples_sequenced(patients=None):
    '''Load patient samples sequenced from general table'''
    from .filenames import get_sample_table_filenames

    sample_table = []
    for fn in get_sample_table_filenames(pnames=patients):
        sample_table.append(pd.read_csv(fn, sep='\t', index_col=0))
    sample_table = pd.concat(sample_table)

    # Note: this refers to the TOTAL # of templates, i.e. the factor 2x for
    # the two parallel RT-PCR reactions
    sample_table['n templates'] = sample_table['viral load'] * 0.4 / 12 * 2

    return [Sample(val) for x, val in sample_table.iterrows()]
