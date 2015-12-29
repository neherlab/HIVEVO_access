# vim: fdm=indent
'''
author:     Fabio Zanini
date:       16/09/15
content:    Module to interface with external data.
'''
# Modules
from pandas import read_csv
import numpy as np
from .HIVreference import HIVreference
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


Li_et_al_protein_translation = {"Matrix":'p17', "Capsid":"p24", "Nucleocapsid":"p7", "p6":"p6",
                                "Protease":"PR","RT":["RT", 'p15'],"Integrase":"IN", "Vif":"vif",
                                "Vpr":"vpr", "Tat":"tat", "Rev":"rev", "Vpu":"vpu",
                                "GP120":"gp120", "GP41":"gp41", "Nef":"nef"}

def load_disorder_scores_HXB2():
    disorder_filename = local_data_folder+'external/Li_Retrovirology_2015/HIVGenome_DisorderScore.csv'
    hxb2 = HIVreference("HXB2")
    dscores = {}

    with open(disorder_filename) as dfile:
        for line in dfile:
            if line.split()[0] in Li_et_al_protein_translation:
                pname = line.split()[0]
                our_pname = Li_et_al_protein_translation[pname]
                npos = int(line.split(',')[0].split()[1][1:])
                print pname, our_pname
                if type(our_pname)==str:
                    print len(hxb2.annotation[our_pname]), len(hxb2.annotation[our_pname])/3, npos
                else:
                    for n in our_pname:
                        print len(hxb2.annotation[n]), len(hxb2.annotation[n])/3, npos
                    our_pname = "_".join(our_pname)
                dscores[our_pname]={}

                line = dfile.next()
                pos = np.array(map(lambda x:int(x[3:]), filter(lambda x:len(x)>3, line.strip().split(','))))-1
                dscores[our_pname]['pos'] = pos

                line = dfile.next()
                val = np.array(map(float, filter(lambda x:len(x)>2, line.strip().split(','))))
                dscores[our_pname]['val'] = val
    # fix RT_p15
    our_pname = "RT_p15"
    if our_pname in dscores:
        dscores['RT'] =  {'pos':dscores[our_pname]['pos'][:len(hxb2.annotation['RT'])/3],
                          'val':dscores[our_pname]['val'][:len(hxb2.annotation['RT'])/3] }
        dscores['p15'] =  {'pos':dscores[our_pname]['pos'][len(hxb2.annotation['RT'])/3:]-len(hxb2.annotation['RT'])/3,
                          'val':dscores[our_pname]['val'][len(hxb2.annotation['RT'])/3:] }

    return dscores


def protein_areaSAS(fnames):
    from collections import defaultdict
    vals = defaultdict(list)
    for fname in fnames:
        with open(fname) as infile:
            for line in infile:
                if line[0]=='\t':
                    pos, val = line.strip().split('\t')
                    pos = int(pos[1:].split('.')[0])
                    val = float(val)
                    vals[pos].append(val)
    flat_vals = sorted([(p-1, np.mean(val)) for p, val in vals.iteritems()])

    return np.array(flat_vals)

def load_accessibility():
    import glob
    acc = {}
    for prot, our_pname in Li_et_al_protein_translation.iteritems():
        fnames = glob.glob('data/external/Li_Retrovirology_2015/SurfaceAreaData/'+prot+'*.txt')
        if fnames:
            if type(our_pname)==list:
                our_pname='_'.join(our_pname)
            acc[our_pname] = protein_areaSAS(fnames)
    return acc






