from moleculekit.molecule import Molecule
from moleculekit.tools.voxeldescriptors import getVoxelDescriptors, viewVoxelFeatures
from moleculekit.tools.atomtyper import prepareProteinForAtomtyping
from moleculekit.smallmol.smallmol import SmallMol
from moleculekit.home import home
import torch
import torch.nn as nn
import torch.nn.functional as F
import numpy as np
import os

tut_data = home(dataDir='/usr/project/csplus2/dietrich/datafiles/mouse')

prot = Molecule(os.path.join(tut_data, 'AF-Q8VGS7-F1-model_v2.pdb.gz'))
prot = prepareProteinForAtomtyping(prot)

prot_vox, prot_centers, prot_N = getVoxelDescriptors(prot, buffer=1)

nchannels = prot_vox.shape[1]

prot_vox_t = prot_vox.transpose().reshape([1, nchannels, prot_N[0], prot_N[1], prot_N[2]])
prot_vox_t = torch.tensor(prot_vox_t.astype(np.float32))

def transform(prot_vox): #reshapes to (nchannels, d, h, w)
    return prot_vox.transpose().reshape([nchannels, prot_N[0], prot_N[1], prot_N[2]])

class Model(nn.Module):
    def __init__(self):
        super(Model, self).__init__()
        self.conv1 = nn.Sequential(
            nn.Conv3d(nchannels, 32, 5),  # (in_channels, out_channels, kernel_size)
            nn.ReLU(inplace=True),
            nn.MaxPool3d(3),
            nn.ReLU(inplace=True),
            nn.Conv3d(32, 64, 3),
            nn.ReLU(inplace=True),
            nn.MaxPool3d(3),
            nn.ReLU(inplace=True),
            nn.Conv3d(64, 128, 3),
            nn.ReLU(inplace=True),
            nn.Flatten(),
            nn.Softmax(1),
        )
        self.fc = nn.Linear(8192, 2)


    def forward(self, x):
       x = self.conv1(x)
       x = self.fc(x)
       return x

transform(prot_vox)

model = Model()
results = model.forward(prot_vox_t)
print(results.shape)
print(nchannels)
