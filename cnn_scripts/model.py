from moleculekit.molecule import Molecule
from moleculekit.tools.voxeldescriptors import getVoxelDescriptors, viewVoxelFeatures
from moleculekit.tools.atomtyper import prepareProteinForAtomtyping
from moleculekit.smallmol.smallmol import SmallMol
from moleculekit.home import home
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import Dataset, DataLoader
import torch.optim as optim
import numpy as np
import os
import pickle
import sklearn

x_data = np.load('x_data.npy')
y_data = np.load('y_data.npy')
x_data, y_data = sklearn.utils.shuffle(x_data, y_data)

#def transform(prot_vox): #reshapes to (nchannels, d, h, w)
    #return prot_vox.transpose().reshape([nchannels, prot_N[0], prot_N[1], prot_N[2]])

test_train_split = 0.7
split_point = int(test_train_split * len(x_data))

class ProtDataset_Train(Dataset):
    def __init__(self):
        self.x = torch.from_numpy(x_data[0 : split_point])
        self.y = torch.from_numpy(y_data[0 : split_point])
        self.n_samples = len(self.x)
    
    def __getitem__(self, index):
        return self.x[index], self.y[index]
    
    def __len__(self):
        return self.n_samples

class ProtDataset_Test(Dataset):
    def __init__(self):
        self.x = torch.from_numpy(x_data[split_point : len(x_data)])
        self.y = torch.from_numpy(y_data[split_point : len(y_data)])
        self.n_samples = len(self.x)
    
    def __getitem__(self, index):
        return self.x[index], self.y[index]
    
    def __len__(self):
        return self.n_samples

class Model(nn.Module):
    def __init__(self):
        super(Model, self).__init__()
        self.conv1 = nn.Sequential(
            nn.Conv3d(8, 32, 5, stride=1, padding=1),  # (in_channels, out_channels, kernel_size)
            nn.ReLU(inplace=True),
            nn.MaxPool3d(3),
            nn.ReLU(inplace=True),
            nn.Conv3d(32, 64, 3, stride=1, padding=1),
            nn.ReLU(inplace=True),
            nn.MaxPool3d(3),
            nn.ReLU(inplace=True),
            nn.Conv3d(64, 128, 3, stride=1, padding=1),
            nn.ReLU(inplace=True),
            nn.Flatten(),
            nn.Softmax(1),
        )
        self.fc = nn.Linear(1024, 1)


    def forward(self, x):
       x = self.conv1(x)
       x = self.fc(x)
       return torch.sigmoid(x)


model = Model()
dataset_train = ProtDataset_Train()
dataset_test = ProtDataset_Test()
dataloader_train = DataLoader(dataset = dataset_train, batch_size = 4, shuffle = True, num_workers = 2)
dataloader_test = DataLoader(dataset = dataset_test, batch_size = 4, shuffle = False, num_workers = 2)


criterion = nn.BCELoss()
optimizer = optim.SGD(model.parameters(), lr=0.001, momentum=0.9)
model = model.float()

dataiter = iter(dataloader_train)
images, labels = dataiter.next()
#print(labels)

for epoch in range(11):  # loop over the dataset multiple times

    running_loss = 0.0
    for i, data in enumerate(dataloader_train, 0):
        # get the inputs; data is a list of [inputs, labels]
        inputs, labels = data
        # zero the parameter gradients
        optimizer.zero_grad()

        # forward + backward + optimize
        outputs = model(inputs.float())
        labels = labels.unsqueeze(1).float()
        loss = criterion(outputs, labels)
        loss.backward()
        optimizer.step()

        # print statistics
        running_loss += loss.item()
        if i % 10 == 0:    # print every 2000 mini-batches
            print(f'[{epoch + 1}, {i + 1:5d}] loss: {running_loss / 10:.3f}')
            running_loss = 0.0

print('finished training')
# save model
PATH = './cnn.pth'
torch.save(model.state_dict(), PATH)

model.eval()

correct = 0
total = 0

# evaluate model
with torch.no_grad():
    for data in dataloader_test:
        images, labels = data
        # calculate outputs by running images through the network
        outputs = model(images.float())
        # the class with the highest energy is what we choose as prediction
        _, predicted = torch.max(outputs.data, 1)
        total += labels.size(0)
        correct += (predicted == labels).sum().item()

print(f'accuracy of the network: {100 * correct // total} %')

