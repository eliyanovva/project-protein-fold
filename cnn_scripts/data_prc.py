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

file = pd.read_csv('./pS6_DE_1p_2e3mp.csv')
mapping = pickle.load(open('ens_to_uniprot.txt', 'rb'))

x_data = []
y_data = []

datapath = home(dataDir='/usr/project/csplus2/dietrich/datafiles/mouse')

#print(file.logFC)
#for i in range(len(file.ensembl_gene_id)):
for i in range(20):
    if file.ensembl_gene_id[i] in mapping:
        if (exists(os.path.join(datapath,
         'AF-' + mapping[file.ensembl_gene_id[i]] + '-F1-model_v2.pdb.gz'))):
            x_data.append(prepareProteinForAtomtyping(Molecule(os.path.join(datapath,
            'AF-' + mapping[file.ensembl_gene_id[i]] + '-F1-model_v2.pdb.gz'))))
            y_data.append(file.logFC[i])
#print(file.ensembl_gene_id)

#classify into 0 or 1 for now
y_min = min(y_data)
y_max = max(y_data)
for i in range(len(y_data)):
    y_data[i] = round((y_data[i] - y_min)/(y_max - y_min))

x_np = np.array([], np.float64)

for i in range(len(x_data)):
    prot_vox, prot_centers, prot_N = getVoxelDescriptors(x_data[i], buffer=1)
    nchannels = prot_vox.shape[1]
    np.append(x_np, np.asarray(prot_vox.transpose().reshape([nchannels, prot_N[0], prot_N[1], prot_N[2]])))

#x_data = np.asarray(x_data)
y_data = np.asarray(y_data)
pickle.dump(x_np, open('x_data.txt', 'wb'))
pickle.dump(y_data, open('y_data.txt', 'wb'))

