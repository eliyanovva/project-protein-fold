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
all_y_data = []
all_x_data = []
vox_dict = dict()
symb_to_uniprot = pickle.load(open( "symb_to_uniprot.p", "rb" ))

def get_lig(cid):
    print('a')

def process_gene(gene_id, cid, y): #appends voxelization to x and bind to y
    mut = ''

    if gene_id.find('_') != -1:
        mut = gene_id[gene_id.find('_'):]
        gene_id = gene_id[:gene_id.find('_')]

    human_path_1 = '/home/users/bmp40/share/deorphanOR_human/'
    human_path_2 = '/home/users/bmp40/share/orphanOR_human/'
    alphafold_path = '/usr/xtmp/derf/alphafold/output/'
    
    if gene_id[1:] not in symb_to_uniprot.keys():
        return False

    uniprot = symb_to_uniprot[gene_id[1:]]

    #print(os.path.join(alphafold_path, uniprot + '/ranked_0.pdb'))

    #print(uniprot)

    if y == 'no':
        bind = 0
    else:
        bind = 1
    if gene_id in vox_dict.keys():
        all_x_data.append(dict)
        all_y_data.append(bind)
    elif exists(os.path.join(human_path_1,'AF-' + uniprot + '-F1-model_v2.pdb')):
        #print('path 1 exists')
        pass
    elif exists(os.path.join(human_path_2,'AF-' + uniprot + '-F1-model_v2.pdb')):
        #print('path 2 exists')
        pass
    elif exists(os.path.join(alphafold_path + uniprot)):
        print('path 3 exists')
        
    else:
        print(uniprot + mut + ' not exists')

for i in range(len(file.Gene)):
    process_gene(file.Gene[i], file.CID[i], file.Responsive[i])
