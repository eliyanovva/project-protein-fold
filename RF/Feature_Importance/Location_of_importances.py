import sys
sys.path.append('../../project-protein-fold/RF/')
import Globals
from collections import Counter

AA_seqs = Globals.initialize_AA_dict(Globals.initialize_protein_list())
Di_seqs = Globals.initialize_3Di_dict(Globals.initialize_protein_list())
tmdict = {'3':0, '5':1, '6':2, '7':3}

def find_feature(tm, seq, dictionary):
    ret = []
    tmind = tmdict[tm]
    for AA in list(dictionary.keys()):
        index = 0
        while index >= 0:
            index = dictionary[AA][tmind].find(seq, index)
            if index >= 0:
                ret.append(index)
                index += 1
    residues = []
    retcount = Counter(ret)
    temp = retcount.most_common(1)[0][1]
    for num in ret:
        if ret.count(num) == temp:
            residues.append(num)
    return set(residues)

with open('Feature_Importance/sulfur_importance.txt') as f:
    lines = f.readlines()
    i = 0
    ret = {}
    for line in lines:
        if i == 50:
            break
        i+=1
        line = line.replace('\n', "")
        if 'TM' in line:
            if line[-1] not in ret:
                ret[line[-1]] = []
            if line[0].isupper():
                ret[line[-1]].append(find_feature(line[-1], line[0:5], Di_seqs))
            else:
                ret[line[-1]].append(find_feature(line[-1], line[0:5], AA_seqs))
    print(ret)