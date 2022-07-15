import sys
sys.path.append('../../project-protein-fold/RF/')
import Globals
from statistics import mode

AA_seqs = Globals.initialize_AA_dict(Globals.initialize_protein_list())
Di_seqs = Globals.initialize_3Di_dict(Globals.initialize_protein_list())

def findAA(tm, seq):
    ret = []
    for AA in list(AA_seqs.keys()):
        index = 0
        while index >= 0:
            index = AA_seqs[AA][tm].find(seq, index)
            if index >= 0:
                ret.append(index)
                index += 1
    return mode(ret)

print(findAA(3, 'aacab'))

def find3Di(tm, seq):
    ret = []
    for Di in list(Di_seqs.keys()):
        index = 0
        while index >= 0:
            index = Di_seqs[Di][tm].find(seq, index)
            if index >= 0:
                ret.append(index)
                index += 1
    return mode(ret)

print(find3Di(2, 'VVLVV'))