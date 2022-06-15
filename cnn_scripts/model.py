from moleculekit.molecule import Molecule
from moleculekit.tools.voxeldescriptors import getVoxelDescriptors, viewVoxelFeatures
from moleculekit.tools.atomtyper import prepareProteinForAtomtyping
from moleculekit.smallmol.smallmol import SmallMol
from moleculekit.home import home
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import Dataset, DataLoader
import numpy as np
import os
import pickle

x_data = pickle.load(open('x_data.txt', 'rb'))

y_data = pickle.load(open('y_data.txt', 'rb'))

#def transform(prot_vox): #reshapes to (nchannels, d, h, w)
    #return prot_vox.transpose().reshape([nchannels, prot_N[0], prot_N[1], prot_N[2]])

class ProtDataset(Dataset):
    def __init__(self):
        self.x = torch.from_numpy(x_data)
        self.y = torch.from_numpy(y_data)
        self.n_samples = len(x_data)
    
    def __getitem__(self, index):
        return self.x[index], self.y[index]
    
    def __len__(self):
        return self.n_samples

class Model(nn.Module):
    def __init__(self):
        super(Model, self).__init__()
        self.conv1 = nn.Sequential(
            nn.Conv3d(8, 32, 5),  # (in_channels, out_channels, kernel_size)
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


model = Model()
dataset = ProtDataset()
print(y_data)