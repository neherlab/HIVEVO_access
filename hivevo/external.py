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


def load_structural_effects_NL43():
    from collections import defaultdict
    import csv
    from .sequence import alphaa
    data_filename = local_data_folder+'external/Carlson_allWithP17Mono.energyPlusPlus.mutationCountCorrection.N.TXT'
    with open(data_filename) as dfile:
        data = csv.DictReader(dfile, delimiter='\t')
        mutations = {}
        for mut in data:
            prot = mut["Protein"].lower()
            if prot not in mutations: mutations[prot] = defaultdict(list)
            mutations[prot][int(mut['Position'])-1].append(map(float, [mut[x] for x in
                            ['RawAbs', 'Abs', 'ExpAbsDDE',  'ExpFreq', 'IhacFreq']]))

    constraint = {}
    cons_seqs = {}
    for prot in mutations:
        cons_seq = []
        tmp = []
        for pos in sorted(mutations[prot].keys()):
            M = np.array(mutations[prot][pos])
            pvec = M[:,-2]
            cons_seq.append(alphaa[M[:,-1].argmax()])
            tmp.append( (pos, -np.sum(pvec*np.log(pvec+1e-20))))
        constraint[prot] = np.array(tmp)
        cons_seqs[prot]=cons_seq
    return constraint, cons_seqs

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


# they use different names for proteins, this maps their names to our names
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
            # following the name, their file contains one line with the positions and one line with the values
            if line.split()[0] in Li_et_al_protein_translation:
                pname = line.split()[0]
                our_pname = Li_et_al_protein_translation[pname]
                npos = int(line.split(',')[0].split()[1][1:])
                if type(our_pname)==list: # deal with multipart maps
                    our_pname = "_".join(our_pname)
                dscores[our_pname]={}

                line = dfile.next()
                # make a ist of the positions, subtract one in the end since they start numbering at 1
                pos = np.array(map(lambda x:int(x[3:]), filter(lambda x:len(x)>3, line.strip().split(','))))-1
                dscores[our_pname]['pos'] = pos

                line = dfile.next()
                # parse the values. the lines are longer than the number of position, hence filter for non-empty
                val = np.array(map(float, filter(lambda x:len(x)>2, line.strip().split(','))))
                dscores[our_pname]['val'] = val
    # RT needs special attention, since their RT corresponds to our RT and p15.
    our_pname = "RT_p15"
    if our_pname in dscores:
        dscores['RT'] =  {'pos':dscores[our_pname]['pos'][:len(hxb2.annotation['RT'])/3],
                          'val':dscores[our_pname]['val'][:len(hxb2.annotation['RT'])/3] }
                          # translate into p15 coordinates by subtracting the length of the RT
        dscores['p15'] =  {'pos':dscores[our_pname]['pos'][len(hxb2.annotation['RT'])/3:]-len(hxb2.annotation['RT'])/3,
                          'val':dscores[our_pname]['val'][len(hxb2.annotation['RT'])/3:] }

    return dscores


def protein_areaSAS(fnames):
    from collections import defaultdict
    vals = defaultdict(list)
    for fname in fnames: #loop over filenames
        with open(fname) as infile:
            for line in infile:
                if line[:2]=='\t:': #data lines start with tab, followed by a colon
                    pos, val = line.strip().split('\t')
                    pos = int(pos[1:].split('.')[0]) # positions contain the chain as :17.A
                    val = float(val)
                    vals[pos].append(val)
    # subtract one from position to obtain zero numbering and take the mean of values from different chains
    flat_vals = sorted([(p-1, np.mean(val)) for p, val in vals.iteritems()])

    return np.array(flat_vals)

def load_accessibility():
    import glob
    acc = {}
    # the accessibility data is provided as many different files, at times several for each
    # protein when structures have been determined piecewise.
    for prot, our_pname in Li_et_al_protein_translation.iteritems():
        fnames = glob.glob(local_data_folder+'/external/Li_Retrovirology_2015/SurfaceAreaData/'+prot+'*.txt')
        if fnames:
            if type(our_pname)==list: #fix multipartite names
                our_pname='_'.join(our_pname)
            #this returns position value pairs which are added to the dictionary
            acc[our_pname] = protein_areaSAS(fnames)

    # fix RT_p15
    our_pname = "RT_p15"
    if our_pname in acc:
        RT = []
        p15 = []
        for pos, val in acc["RT_p15"]:
            if pos<440: #440 is the start of p15 when counting from the beginning of RT
                RT.append((pos,val))
            else:
                p15.append((pos-440,val))

        acc['RT']=np.array(RT)
        acc['p15']=np.array(p15)
    return acc






