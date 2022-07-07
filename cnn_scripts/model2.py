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
from sklearn.model_selection import train_test_split
import random
from scipy.ndimage.interpolation import rotate
from tqdm import tqdm

x_data = np.load('x_data.npy')
y_data = np.load('y_data.npy')
lig_data = np.load('lig_data.npy')

#def transform(prot_vox): #reshapes to (nchannels, d, h, w)
    #return prot_vox.transpose().reshape([nchannels, prot_N[0], prot_N[1], prot_N[2]])

train_x, test_x, train_lig, test_lig, train_y, test_y = train_test_split(x_data, lig_data, y_data,  test_size = 0.2, stratify=y_data)

final_train_x = []
final_train_lig = []
final_train_y = []

for i in tqdm(range(train_x.shape[0])):
    final_train_x.append(train_x[i])
    final_train_lig.append(train_lig[i])
    final_train_y.append(train_y[i])
    print('reached')
    for j in range(1):
        final_train_x.append(rotate(train_x[i], angle=random.randrange(-15, 15), reshape=False))
        final_train_lig.append(rotate(train_lig[i], angle=random.randrange(-15, 15), reshape=False))
        final_train_y.append(train_y[i])


class ProtDataset(Dataset):
    def __init__(self):
        self.x = torch.from_numpy(np.load('x_data.npy'))
        #self.lig = torch.from_numpy(np.load('lig_data.npy'))
        self.y = torch.from_numpy(np.load('y_data.npy'))
        self.n_samples = len(self.x)
    
    def __getitem__(self, index):
        #return self.x[index], self.lig[index], self.y[index]
        return self.x[index], self.y[index]
    
    def __len__(self):
        return self.n_samples

class TrainDataset(Dataset):
    def __init__(self):
        self.x = torch.from_numpy(np.asarray(final_train_x, dtype=np.float32))
        self.lig = torch.from_numpy(np.asarray(final_train_lig, dtype=np.float32))
        self.y = torch.from_numpy(np.asarray(final_train_y, dtype=np.float32))
        self.n_samples = len(self.x)
    
    def __getitem__(self, index):
        #return self.x[index], self.lig[index], self.y[index]
        return self.x[index], self.lig[index], self.y[index]
    
    def __len__(self):
        return self.n_samples
    
class TestDataset(Dataset):
    def __init__(self):
        self.x = torch.from_numpy(np.asarray(test_x, dtype=np.float32))
        self.lig = torch.from_numpy(np.asarray(test_lig))
        self.y = torch.from_numpy(np.asarray(test_y, dtype=np.float32))
        self.n_samples = len(self.x)
    
    def __getitem__(self, index):
        #return self.x[index], self.lig[index], self.y[index]
        return self.x[index], self.lig[index], self.y[index]
    
    def __len__(self):
        return self.n_samples

class Model(nn.Module):
    def __init__(self):
        super(Model, self).__init__()
        self.conv1 = nn.Sequential(
            #convolution 1
            nn.Conv3d(16, 32, 5, stride=1, padding=0),  # (in_channels, out_channels, kernel_size)
            nn.MaxPool3d(2),
            nn.ReLU(),
            #convolution 2
            #nn.Conv3d(32, 64, 3, stride=3, padding=3),
            #nn.MaxPool3d(2),
            #nn.ReLU(),
            #convolution 3
            #nn.Conv3d(64, 128, 3, stride=3, padding=3),
            #nn.MaxPool3d(2),
            #nn.ReLU(),
            nn.Flatten()
            #nn.Softmax(1)
        )
        self.fc = nn.Sequential(
            #nn.Linear(576, 1000),
            #nn.Dropout(),
            #nn.ReLU(),
            nn.Linear(186624, 1)
        )


    def forward(self, x, lig):
       x = torch.cat((x, lig), dim=1)
       x = self.conv1(x)
       x = self.fc(x)
       return torch.sigmoid(x)


model = Model()

dataset_full = ProtDataset()
train_test_split = 0.7
train_size = int(train_test_split * len(dataset_full))
test_size = len(dataset_full) - train_size
#dataset_train, dataset_test = torch.utils.data.random_split(dataset_full, [train_size, test_size])
dataset_train = TrainDataset()
dataset_test = TestDataset()

dataloader_train = DataLoader(dataset = dataset_train, batch_size = 4, shuffle = True, num_workers = 2)
dataloader_test = DataLoader(dataset = dataset_test, batch_size = 4, shuffle = False, num_workers = 2)

print('training length: ', len(dataset_train))

criterion = nn.BCELoss()
optimizer = optim.SGD(model.parameters(), lr=0.001, momentum=0.9)
model = model.float()

#print(labels)

total_epochs = 20

for epoch in range(total_epochs):  # loop over the dataset multiple times

    running_loss = 0.0
    for i, data in enumerate(dataloader_train, 0):
        model.train()
        # get the inputs; data is a list of [prot, lig, labels]
        prot, lig, labels = data
        # zero the parameter gradients
        optimizer.zero_grad()

        # forward + backward + optimize
        outputs = model(prot.float(), lig.float())
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
        prot, lig, labels = data
        # calculate outputs by running images through the network
        outputs = model(prot.float(), lig.float())
        #print('outputs:', outputs)
        print('labels:', labels)
        # the class with the highest energy is what we choose as prediction
        outputs = torch.round(outputs)
        print('rounded output:', outputs)
        for i in range(len(outputs)):
            total += 1
            if (outputs[i] == labels[i]):
                 correct += 1
        #
        #_, predicted = torch.max(outputs.data, 1)
        #total += labels.size(0)
        #correct += (predicted == labels).sum().item()

print(f'accuracy of the network: {100 * correct // total} %')

