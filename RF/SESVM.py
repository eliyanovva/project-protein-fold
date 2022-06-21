#based on algorithm description from https://dl.acm.org/doi/pdf/10.1145/3307339.3342141
#random documentation from https://www.geeksforgeeks.org/random-numbers-in-python/

from sklearn import svm
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

def seperate_sets(features, labels):
    pos_set = []
    neg_set = []
    for i in range(len(labels)):
        if labels[i] == 0:
            neg_set.append(features[i])
        else:
            pos_set.append(features[i])
    return pos_set, neg_set

p, n = seperate_sets(X, Y)
print(p)
print(n)

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


def SESVM(pos_set, neg_set):
    M = np.ceil(float(len(neg_set)) / float(len(pos_set)))
    parts = create_partitions(pos_set, neg_set, M)
    for set in parts:
        print(set)
    #for m in range(M):
        #use GridSearchCV to optimize rbf hyperparameters???
        #C: decreasing C => more regulation, decrease C if there's a lot of noise
        #gamma
        #theta = svm.SVC(kernel='rbf')

SESVM(p, n)