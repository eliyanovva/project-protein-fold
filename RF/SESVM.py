#based on algorithm description from https://dl.acm.org/doi/pdf/10.1145/3307339.3342141
#random documentation from https://www.geeksforgeeks.org/random-numbers-in-python/

from sklearn import svm
from sklearn.metrics import accuracy_score
import random
import numpy as np

#N+ = positive training set
#N- = negative training set
#N = training set = {N+, N-}
#P = test set
"""
For T iterations:
    randomly partition N- into M
    (assume N- can be paritioned into M sets without overlapping)
    M = |N-| / |N+| rounded up
    for m in [1, M]:
        S_m = {N+, N-_m}
        train SVM with Radial basis kernel (RBF) on S_m
        (-)m is the classifier that was trained
        a_m = prediction accuracy of (-)m on the entire set N
    Select the best classifer (-)*_t such that:
        (-)*_t = (-)k; k = arg max{a1, ..., aM}
    P* e {0,1} ^ |P| = prediction results
    P*_t = prediction vector of (-)*_t for the testing set P
P^* = MV(P*_1, P*2, ..., P*_T)
"""

#X = features
#Y = labels

X = [[1, 1, 1], [2, 1, 1], [3, 1, 1], [4, 1, 1], [5, 1, 1], [6, 1, 1], [7, 1, 1],
     [8, 1, 1], [9, 1, 1], [10, 1, 1], [11, 1, 1], [12, 1, 1], [13, 1, 1], [14, 1, 1]]
Y = [0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0]

def seperate_sets(N, Y):
    pos_set = []
    neg_set = []
    for i in range(len(Y)):
        if Y[i] == 0:
            neg_set.append(N[i])
        else:
            pos_set.append(N[i])
    return pos_set, neg_set

def create_partitions(pos_set, neg_set, M):
    partitions = []
    part_len = len(pos_set)
    overhang = part_len * M - len(neg_set)
    duplicate = []
    for item in neg_set:
        duplicate.append(item)

    for i in range(int(overhang)):
        dup_row = random.choice(duplicate)
        neg_set.append(dup_row)
        duplicate.remove(dup_row)

    for i in range(int(M)):
        n_set = []
        for j in range(part_len):
            row = random.choice(neg_set)
            n_set.append(row)
            neg_set.remove(row)
        partitions.append(n_set)
    return partitions


def SESVM(N, Y):
    pos_set, neg_set = seperate_sets(N, Y)
    M = np.ceil(float(len(neg_set)) / float(len(pos_set)))
    parts = create_partitions(pos_set, neg_set, M)
    labels = np.append(np.repeat(1, len(pos_set)), np.repeat(0, len(pos_set)))

    all_thetas = []
    accuracies = []

    for m in range(int(M)):
        #use GridSearchCV to optimize rbf hyperparameters???
        #C: decreasing C => more regulation, decrease C if there's a lot of noise
        #gamma
        features = np.concatenate((pos_set, parts[m]), axis=0)
        theta = svm.SVC(kernel='rbf')
        t = theta.fit(features, labels)
        all_thetas.append(t)
        a = accuracy_score(Y, theta.predict(N))
        accuracies.append(a)

    max_accuracy = 0
    max_index = 0

    for m in range(int(M)):
        if accuracies[m] > max_accuracy:
            max_accuracy = accuracies[m]
            max_index = m

    max_theta = all_thetas[max_index]

SESVM(X, Y)
