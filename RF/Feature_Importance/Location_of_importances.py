import sys
sys.path.append('../../project-protein-fold/RF/')
import Globals
from statistics import mode

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
    return mode(ret)

print(find_feature('7', 'aacab', AA_seqs))

print(find_feature('6', 'VVLVV', Di_seqs))


with open('Feature_Importance/sulfur_importance.txt') as f:
    lines = f.readlines()
    for line in lines:
        line = line.replace('\n', "")
        if 'TM' in line:
            if line[0].isupper():
                find_feature(line[-1], line[0:5], Di_seqs)
            else:
                find_feature(line[-1], line[0:5], AA_seqs)