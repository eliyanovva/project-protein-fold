#This script follows the combination of k-means and random forest strategy applied in this paper: https://www.hindawi.com/journals/geofluids/2021/9321565/
 
#Imports
import matplotlib.pyplot as plt
from kneed import KneeLocator
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from sklearn.preprocessing import StandardScaler
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
kmeans_kwargs = {"verbose" : "1"} #Print SSE

sse = [] #list holding SSE values
for k in range(1,40): #Test different number of clusters
    kmeans = KMeans(n_clusters=k, **kmeans_kwargs)
    kmeans.fit(X_train)
    sse.append(kmeans.inertia_)

#See how k decreases sse
figure = plt.plot(range(1,40), sse)
plt.savefig("ktest.png")