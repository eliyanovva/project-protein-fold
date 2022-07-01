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

file = pd.read_csv('/home/users/bmp40/project-protein-fold/cnn_scripts/olfr_de/pS6_DE_1p_heptanal.csv')
mapping = pickle.load(open('ens_to_uniprot.txt', 'rb'))

x_data = []
y_data = []
lig_data = []

lig = SmallMol('/home/users/bmp40/project-protein-fold/mol_data_files/heptanal.mol', force_reading=True)
lig_vox, lig_centers, lig_N = getVoxelDescriptors(lig)

datapath = home(dataDir='/home/users/bmp40/mouse')
box = [20, 20, 20]
cent = [10, 10, 10]
#print(file.logFC)
for i in range(len(file.ensembl_gene_id)):
#for i in range(100):
    if file.ensembl_gene_id[i] in mapping:
        path = os.path.join(datapath,
         'AF-' + mapping[file.ensembl_gene_id[i]] + '-F1-model_v2.pdb')
        if (not exists(path)):
            os.system('gunzip --keep ' + path + '.gz')
            if (exists(path)):
                lig_data.append(lig_vox.transpose().reshape([lig_vox.shape[1], lig_N[0], lig_N[1], lig_N[2]]))
                x_data.append(prepareProteinForAtomtyping(Molecule(os.path.join(datapath, 'AF-' + mapping[file.ensembl_gene_id[i]] + '-F1-model_v2.pdb'))))
                y_data.append(file.logFC[i])
        else:
            lig_data.append(lig_vox.transpose().reshape([lig_vox.shape[1], lig_N[0], lig_N[1], lig_N[2]]))
            x_data.append(prepareProteinForAtomtyping(Molecule(os.path.join(datapath, 'AF-' + mapping[file.ensembl_gene_id[i]] + '-F1-model_v2.pdb'))))
            y_data.append(file.logFC[i])
#print(file.ensembl_gene_id)

#classify into 0 or 1 for now


"""
y_min = min(y_data)
y_max = max(y_data)
for i in range(len(y_data)):
    y_data[i] = round((y_data[i] - y_min)/(y_max - y_min))
    print(y_data[i])
"""


y_med = statistics.median(y_data)
for i in range(len(y_data)):
    if (y_data[i] >= y_med):
         y_data[i] = 1
    else:
        y_data[i] = 0

prot_vox, prot_centers, prot_N = getVoxelDescriptors(x_data[0], boxsize=box, center=cent)
print(prot_vox.shape[1])
dim1 = prot_N[0]
dim2 = prot_N[1]
dim3 = prot_N[2]

for i in range(len(x_data)):
    prot_vox, prot_centers, prot_N = getVoxelDescriptors(x_data[i], boxsize=box, center=cent)
    print(prot_vox)
    nchannels = prot_vox.shape[1]
    x_data[i] = prot_vox.transpose().reshape([nchannels, prot_N[0], prot_N[1], prot_N[2]])
    print(prot_N[0], prot_N[1], prot_N[2])

#x_data = np.asarray(x_data)
print(len(x_data))
y_np = np.asarray(y_data)
x_np = np.asarray(x_data)
lig_np = np.asarray(lig_data)
np.save('x_data', x_np)
np.save('y_data', y_np)
np.save('lig_data', lig_np)