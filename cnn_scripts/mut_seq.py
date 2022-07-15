import csv
import os
from os.path import exists
import pickle

symb_to_uniprot = {}
human_path1 = '/home/users/bmp40/share/deorphanOR_human/'
human_path2 = '/home/users/bmp40/share/orphanOR_human/'

with open('results.txt', 'r') as file:
    for line in file:
        data = line.split()
        if (data[0] == 'Status'):
            continue
        symb_to_uniprot.update({data[1]: data[10]})

symb_to_uniprot.update({'OR13C7': 'P0DN81'})
symb_to_uniprot.update({'OR52A4P': 'A6NMU1'})
symb_to_uniprot.update({'OR3A4P': 'P47883'})
symb_to_uniprot.update({'OR5AL1': 'P0C617'})
symb_to_uniprot.update({'OR11H7': 'Q8NGC8'})
symb_to_uniprot.update({'OR4F14': 'A0A126GV89'})
symb_to_uniprot.update({'OR8U3': 'Q8NH85'})
symb_to_uniprot.update({'OR1F12P': 'Q8NHA8'})
symb_to_uniprot.update({'OR8K3': 'Q8NH51'})
symb_to_uniprot.update({'OR14L1': 'A0A126GV55'})
symb_to_uniprot.update({'OR5AK3P': 'Q8NH89'})
symb_to_uniprot.update({'OR4A4': 'A0A126GWA7'})
symb_to_uniprot.update({'OR4K3': 'Q96R72'})
symb_to_uniprot.update({'OR52P1': 'Q8NH57'})
symb_to_uniprot.update({'OR52L2P': 'Q8NGH6'})
symb_to_uniprot.update({'OR5S1': 'A0A126GVL2'})
symb_to_uniprot.update({'OR2B8P': 'P59922'})
symb_to_uniprot.update({'OR11H13': 'A0A126GW28'})
symb_to_uniprot.update({'OR52E1': 'Q8NGJ3'})
symb_to_uniprot.update({'OR2T32': 'A0A126GWP8'})
symb_to_uniprot.update({'OR1E3': 'Q8WZA6'})
symb_to_uniprot.update({'OR4C45': 'A6NMZ5'})
symb_to_uniprot.update({'OR10A6': 'Q8NH74'})
symb_to_uniprot.update({'OR2W5P': 'A6NFC9'})

set = set()
uniprot_to_seq = {}
uniprot_to_seq_new = {}

pickle.dump(symb_to_uniprot, open( "symb_to_uniprot.p", "wb" ))


with open('sequences.fasta', 'r') as file:
    sequence = ''
    uniprot = ''
    for line in file:
        if line[0] == '>': #new prot
            if (uniprot != ''):
                uniprot_to_seq.update({uniprot: sequence})
                if not exists(os.path.join(human_path1,'AF-' + uniprot + '-F1-model_v2.pdb')) and not exists(os.path.join(human_path2,'AF-' + uniprot + '-F1-model_v2.pdb')):
                    uniprot_to_seq_new.update({uniprot: sequence})

            sequence = ''
            uniprot = line.split('|')[1]
        else:
            sequence += line
            if sequence[-1] == '\n':
                sequence = sequence[0:-1]

print(len(uniprot_to_seq))
print(len(uniprot_to_seq_new))

def dict_to_fasta(d, name):
    """
    Writes a fasta file from a dictionary where keys are sequence names and
    values are sequences
    :param d: dictionary
    :param dir: directory for fasta file to be stored
    :param name: filename to write
    :return:
    """
    with open(f'{name}.fa', 'w') as f:
        for k, v in d.items():
            f.write(f'>{k}\n')
            f.write(f'{v}\n')

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

with open('odorants_ORs_paper.csv') as csv_file:
    for row in csv.reader(csv_file, delimiter=','):
        mutation = ''
        if '_' in row[0]:
           symbol = row[0][1:row[0].index('_')]
           mutation = row[0][row[0].index('_')+1:]
        else:
            symbol = row[0][1:]
        #print(symbol)
        if (symbol == 'ene'):
            continue
        
        #print(mutation)

        if symbol in symb_to_uniprot:
            unip = symb_to_uniprot[symbol]
            if unip in uniprot_to_seq:
                if mutation != '' and mutation != 'Indel' and (unip + '_' + mutation) not in uniprot_to_seq and mutate(uniprot_to_seq[unip], mutation) != 'invalid':
                    uniprot_to_seq.update({unip + '_' + mutation: mutate(uniprot_to_seq[unip], mutation)})
                    if unip in uniprot_to_seq_new:
                        uniprot_to_seq_new.update({unip + '_' + mutation: mutate(uniprot_to_seq[unip], mutation)})



            #print(exists(path1) or exists(path2))
print(len(uniprot_to_seq))
print(len(uniprot_to_seq_new))
dict_to_fasta(uniprot_to_seq, 'seq_dict')
dict_to_fasta(uniprot_to_seq_new, 'new_seq_dict')

#print(symb_to_uniprot)


