import pickle
import os
from os.path import exists

symb_to_uniprot = pickle.load(open("symb_to_uniprot.p", "rb"))
print('OR10A2' in symb_to_uniprot.keys())
print(symb_to_uniprot['OR10A2'])
uniprot_to_seq = {}
with open('sequences.fasta', 'r') as file:
    sequence = ''
    uniprot = ''
    for line in file:
        if line[0] == '>': #new prot
            if (uniprot != ''):
                uniprot_to_seq.update({uniprot: sequence})
            sequence = ''
            uniprot = line.split('|')[1]
        else:
            sequence += line
            if sequence[-1] == '\n':
                sequence = sequence[0:-1]

print(uniprot_to_seq['Q9H208'])


with open('./no_pdb.txt') as file:
    for line in file:
        print(line.rstrip())


def mutate(seq, mut):
    ret = []
    for c in seq:
        ret.append(c)
    muts = mut.split('_')
    for m in muts:
        index = int(m[1:-1]) - 1
        if len(seq) <= index:
            return 'invalid'
        ret[index] = m[-1]
    return ''.join(ret)
