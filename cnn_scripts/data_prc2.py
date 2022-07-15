import numpy as np
import pandas as pd
import pickle
from moleculekit.molecule import Molecule
from moleculekit.tools.voxeldescriptors import getVoxelDescriptors, viewVoxelFeatures
from moleculekit.tools.atomtyper import prepareProteinForAtomtyping
from moleculekit.smallmol.smallmol import SmallMol
from moleculekit.home import home
import os
from os.path import exists
import torch
import gzip
import statistics

file = pd.read_csv('/home/users/bmp40/project-protein-fold/cnn_scripts/odorants_ORs_paper.csv')
cid_file = pd.read_csv('/home/users/bmp40/project-protein-fold/cnn_scripts/odorants_paper.csv', encoding='latin-1')
all_y_data = []
all_x_data = []
vox_dict = dict()
symb_to_uniprot = pickle.load(open( "symb_to_uniprot.p", "rb" ))
cid_dict = dict()

def process_gene(gene_id, y): #appends voxelization to x and bind to y

    mut = ''

    if gene_id.find('_') != -1:
        mut = gene_id[gene_id.find('_'):]
        gene_id = gene_id[:gene_id.find('_')]

    human_path_1 = '/home/users/bmp40/share/deorphanOR_human/'
    human_path_2 = '/home/users/bmp40/share/orphanOR_human/'
    if gene_id[1:] not in symb_to_uniprot.keys():
        return False
    uniprot = symb_to_uniprot[gene_id[1:]]

    #print(uniprot)

    if y == 'no':
        bind = 0
    else:
        bind = 1
    if gene_id in vox_dict.keys():
        all_x_data.append(dict)
        all_y_data.append(bind)
    #elif exists(os.path.join(human_path_1,'AF-' + gene_id + '-F1-model_v2.pdb')):
        #print('exists')
    #else:
        #print('not exists')

print(cid_file.columns)
cid_file = cid_file.rename(columns={'Number ': 'Number'})

for i in range(len(cid_file.Odorant)):
    cid_dict.update({file.Odorant[i]: file.PubChem_CID[i]})

for i in range(len(file.Gene)):
    process_gene(file.Gene[i], file.Responsive[i])
