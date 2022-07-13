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
from pandas import *
from torchsampler import ImbalancedDatasetSampler
#from torchmetrics import AUROC

print(torch.cuda.is_available())
print(torch.cuda.device_count())
#x_data = np.load('x_data.npy')
#y_data = np.load('y_data.npy')
#lig_data = np.load('lig_data.npy')
#total_data = []
x_data = []
lig_data = []
y_data = []
total_data = []

#read in data
for filename in tqdm(os.listdir('./model_files')):
    data = np.load(os.path.join('./model_files', filename), allow_pickle=True)
    for obj in tqdm(data):
            total_data = np.concatenate((total_data, obj), axis=None)

for i in tqdm(range(len(total_data))):
    x_data.append(total_data[i]['Protein Data'].numpy())
    lig_data.append(total_data[i]['Ligand Data'].numpy())
    y_data.append(total_data[i]['Binding'])

x_data = np.asarray(x_data)
lig_data = np.asarray(lig_data)
y_data = np.asarray(y_data)

#x_data = torch.from_numpy(x_data)
#lig_data = torch.from_numpy(lig_data)
#y_data = torch.from_numpy(y_data)

#print(x_data.size())
#print(lig_data.size())
#print(y_data.size())

train_x, test_x, train_lig, test_lig, train_y, test_y = train_test_split(x_data, lig_data, y_data,  test_size = 0.2, stratify=y_data)

final_train_x = []
final_train_lig = []
final_train_y = []

for i in tqdm(range(train_x.shape[0])):
    final_train_x.append(train_x[i])
    final_train_lig.append(train_lig[i])
    final_train_y.append(train_y[i])
    for j in range(0):
        final_train_x.append(rotate(train_x[i], angle=random.randrange(-15, 15), reshape=False))
        final_train_lig.append(rotate(train_lig[i], angle=random.randrange(-15, 15), reshape=False))
        final_train_y.append(train_y[i])

final_train_x = torch.from_numpy(np.asarray(final_train_x))
final_train_lig = torch.from_numpy(np.asarray(final_train_lig))
final_train_y = torch.from_numpy(np.asarray(final_train_y))

test_x = torch.from_numpy(np.asarray(test_x))
test_lig = torch.from_numpy(np.asarray(test_lig))
test_y = torch.from_numpy(np.asarray(test_y))

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
        #self.x = torch.from_numpy(np.asarray(final_train_x, dtype=np.float32))
        #self.lig = torch.from_numpy(np.asarray(final_train_lig, dtype=np.float32))
        #self.y = torch.from_numpy(np.asarray(final_train_y, dtype=np.float32))
        self.x = final_train_x
        self.lig = final_train_lig
        self.y = final_train_y
        self.n_samples = len(self.x)
    
    def __getitem__(self, index):
        #return self.x[index], self.lig[index], self.y[index]
        return self.x[index], self.lig[index], self.y[index]
    
    def __len__(self):
        return self.n_samples
    
class TestDataset(Dataset):
    def __init__(self):
        #self.x = torch.from_numpy(np.asarray(test_x, dtype=np.float32))
        #self.lig = torch.from_numpy(np.asarray(test_lig))
        #self.y = torch.from_numpy(np.asarray(test_y, dtype=np.float32))
        self.x = test_x
        self.lig = test_lig
        self.y = test_y
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
            nn.Conv3d(16, 32, 5),  # (in_channels, out_channels, kernel_size)
            nn.MaxPool3d(2),
            nn.ReLU(),
            #convolution 2
            nn.Dropout(p=0.2),
            nn.Conv3d(32, 64, 3),
            nn.MaxPool3d(2),
            nn.ReLU(),
            #convolution 3
            #nn.Conv3d(64, 128, 3, stride=3, padding=3),
            #nn.MaxPool3d(2),
            #nn.ReLU(),
            nn.Flatten()
            #nn.Softmax(1)
        )
        self.fc = nn.Sequential(
            nn.Linear(32768, 128),
            nn.Dropout(),
            nn.ReLU(),
            nn.Linear(128, 1)
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

#dataloader_train = DataLoader(dataset = dataset_train, batch_size = 4, shuffle = True, num_workers = 2)
dataloader_train = DataLoader(dataset = dataset_train, batch_size = 4, shuffle = True, sampler=ImbalancedDatasetSampler(dataset_train), num_workers = 2)

dataloader_test = DataLoader(dataset = dataset_test, batch_size = 4, shuffle = False, num_workers = 2)

print('training length: ', len(dataset_train))
#custom weighted bce for imbalanced dataset
#def weighted_binary_cross_entropy(output, target, weights=None):
    #output = torch.clamp(output,min=1e-8,max=1-1e-8)

    #if weights is not None:
    #   assert len(weights) == 2
    #    loss = weights[1] * (target * torch.log(output)) + \
    #           weights[0] * ((1 - target) * torch.log(1 - output))
    #else:
    #    loss = target * torch.log(output) + (1 - target) * torch.log(1 - output)

    #return torch.neg(torch.mean(loss))                 

criterion = nn.BCELoss()
#criterion = weighted_binary_cross_entropy
optimizer = optim.SGD(model.parameters(), lr=0.001, momentum=0.9)
model = model.float()

#print(labels)

total_epochs = 3
print('dataloader length: ')
print(len(dataloader_train))

for epoch in tqdm(range(total_epochs)):  # loop over the dataset multiple times
    running_loss = 0.0

    for i, data in tqdm(enumerate(dataloader_train, 0)):
        model.train()
        # get the inputs; data is a list of [prot, lig, labels]
        prot, lig, labels = data
        # zero the parameter gradients
        optimizer.zero_grad()

        # forward + backward + optimize
        outputs = model(prot.float(), lig.float())
        labels = labels.unsqueeze(1).float()
        loss = criterion(outputs, labels) #switch maybe?

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

#auroc = AUROC(pos_label=1)
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
        
        #print(auroc(outputs, labels))
        #
        #_, predicted = torch.max(outputs.data, 1)
        #total += labels.size(0)
        #correct += (predicted == labels).sum().item()

print(f'accuracy of the network: {100 * correct // total} %')

