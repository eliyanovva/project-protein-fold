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
import shutil

file = pd.read_csv('/home/users/bmp40/project-protein-fold/cnn_scripts/odorants_ORs_paper.csv')
all_y_data = []
all_x_data = []
vox_dict = dict()
symb_to_uniprot = pickle.load(open( "symb_to_uniprot.p", "rb" ))

print(symb_to_uniprot['OR10A2'])

symb_to_filepath = dict()

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
        symb_to_filepath.update({gene_id : human_path_1 +'AF-' + uniprot + '-F1-model_v2.pdb'})
        shutil.copyfile(os.path.join(human_path_1,'AF-' + uniprot + '-F1-model_v2.pdb'), './all_pdb/AF-' + uniprot + '-F1-model_v2.pdb')
        pass
    elif exists(os.path.join(human_path_2,'AF-' + uniprot + '-F1-model_v2.pdb')):
        shutil.copyfile(os.path.join(human_path_2,'AF-' + uniprot + '-F1-model_v2.pdb'), './all_pdb/AF-' + uniprot + '-F1-model_v2.pdb')
        #print('path 2 exists')
        symb_to_filepath.update({gene_id : human_path_2 +'AF-' + uniprot + '-F1-model_v2.pdb'})
        pass
    if exists(os.path.join(alphafold_path + uniprot.lower() + mut.lower() + '/ranked_0.pdb')):
        #print('path 3 exists')
        symb_to_filepath.update({gene_id + mut : alphafold_path + uniprot.lower() + mut.lower() + '/ranked_0.pdb'})
        if (mut == '_H43R_H207R_K258T'):
            print('exists')
        shutil.copyfile(os.path.join(alphafold_path + uniprot.lower() + mut.lower()) + '/ranked_0.pdb', './all_pdb/AF-' + uniprot + mut + '-F1-model_v2.pdb')
        pass
    else:
        #print(uniprot + mut + ' not exists')
        pass

for i in range(len(file.Gene)):
    process_gene(file.Gene[i], file.CID[i], file.Responsive[i])

#print(symb_to_filepath)
print (symb_to_filepath.keys())
pickle.dump(symb_to_filepath, open( "symb_to_filepath.p", "wb" ) )