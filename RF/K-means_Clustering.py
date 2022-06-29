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


#Faster k-means (from: https://towardsdatascience.com/k-means-8x-faster-27x-lower-error-than-scikit-learns-in-25-lines-eaedc7a3a0c8)
import faiss
import numpy as np

class FaissKMeans:
    def __init__(self, n_clusters, n_init=10, max_iter=300):
        self.n_clusters = n_clusters
        self.n_init = n_init
        self.max_iter = max_iter
        self.kmeans = None
        self.cluster_centers_ = None
        self.inertia_ = None

    def fit(self, X):
        self.kmeans = faiss.Kmeans(d=X.shape[1],
                                   k=self.n_clusters,
                                   niter=self.max_iter,
                                   nredo=self.n_init)
        self.kmeans.train(X.astype(np.float32))
        self.cluster_centers_ = self.kmeans.centroids
        self.inertia_ = self.kmeans.obj[-1]

    def predict(self, X):
        return self.kmeans.index.search(X.astype(np.float32), 1)[1]

#Cluster the data 
sse = [] #list holding SSE values
silhouette_coefficients = []
for k in range(2,41): #Test different number of clusters
    kmeans = FaissKMeans(n_clusters=k)
    kmeans.fit(X_train)
    sse.append(kmeans.inertia_)
    score = silhouette_score(X_train, kmeans.labels_)
    silhouette_coefficients.append(score)

#See how k decreases sse
figure = plt.plot(range(2,41), sse)
plt.savefig("ktest.png")