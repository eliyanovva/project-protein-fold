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
import numpy as np

datapath = home(dataDir='/home/users/bmp40/project-protein-fold/pdb_data_files')

path = os.path.join(datapath, 'AF-Q8VGS3-F1-model_v2.pdb')

protein = Molecule(path)

protein = prepareProteinForAtomtyping(protein)

prot_vox, prot_centers, prot_N = getVoxelDescriptors(protein, buffer=1)

np_arr = np.array(prot_vox)
np.set_printoptions(threshold=np.inf)

print(np_arr)
