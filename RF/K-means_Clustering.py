#This script follows the combination of k-means and random forest strategy applied in this paper: https://www.hindawi.com/journals/geofluids/2021/9321565/
 
#Imports
import matplotlib.pyplot as plt
from kneed import KneeLocator
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
import numpy as np
import CombineLigandsProteins
from sklearn.model_selection import train_test_split

#Import data
CombineLigandsProteins.import_final()
testX = CombineLigandsProteins.X
testY = CombineLigandsProteins.Y

#split into training and test set
X_train, X_test, y_train, y_test = train_test_split(testX, testY, stratify=testY, test_size=0.1) # 90% training and 10% test

#Cluster the data 
import faiss
import numpy as np

d = X_train.shape[1]
niter = 20
verbose = True
sse = [] #list holding SSE values
for k in range(2,2000,200): #Test different number of clusters
    ncentroids=k
    kmeans = faiss.Kmeans(d, ncentroids, niter = niter, verbose = verbose)
    kmeans.train(X_train.astype(np.float32))
    sse.append(kmeans.obj)
"""
sse = [] #list holding SSE values
for k in range(2,2000,200): #Test different number of clusters
    ncentroids=k
    kmeans = KMeans(n_clusters=k, verbose=1)
    kmeans.fit(X_train)
    sse.append(kmeans.inertia_)
"""

#See how k decreases sse
figure = plt.plot(range(2,2000,200), sse)
plt.savefig("ktest.png")