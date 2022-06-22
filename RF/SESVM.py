#based on algorithm description from https://dl.acm.org/doi/pdf/10.1145/3307339.3342141
#random documentation from https://www.geeksforgeeks.org/random-numbers-in-python/

from sklearn import svm
from sklearn.metrics import accuracy_score
import random
import numpy as np
import CombineLigandsProteins

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

"""will be positive if:
    1st: geq than 6
    2nd: leq than 5
    3rd: quotient of 1st and 2nd
"""

X = [[8, 4, 2], [7, 6, 3], [8, 5, 3], [9, 3, 3], [5, 2, 3], [4, 3, 7], [8, 3, 1], [7, 2, 9], [3, 5, 6], [1, 4, 2], [8, 3, 2], [8, 6, 1]]
Y = [1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0]
#2-3 pos in X, 1 pos in P_X
P_X = [[5, 4, 1], [10, 5, 2], [8, 4, 1]]
P_Y = [0, 1, 0]

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


def SESVM(N, Y, T, P_X, P_Y):
    pos_set, neg_set = seperate_sets(N, Y)
    M = np.ceil(float(len(neg_set)) / float(len(pos_set)))
    predictions = []

    for i in range(T):
        neg_copy = []
        for item in neg_set:
            neg_copy.append(item)

        parts = create_partitions(pos_set, neg_copy, M)
        labels = np.append(np.repeat(1, len(pos_set)), np.repeat(0, len(pos_set)))

        all_features = []
        accuracies = []

        for m in range(int(M)):
            #use GridSearchCV to optimize rbf hyperparameters???
            #C: decreasing C => more regulation, decrease C if there's a lot of noise
            #gamma
            features = np.concatenate((pos_set, parts[m]), axis=0)
            theta = svm.SVC(kernel='rbf')
            theta.fit(features, labels)
            all_features.append(features)
            p = theta.predict(N)
            a = accuracy_score(Y, p)
            accuracies.append(a)

        max_accuracy = 0
        max_index = 0

        for m in range(int(M)):
            if accuracies[m] > max_accuracy:
                max_accuracy = accuracies[m]
                max_index = m

        opt_theta = svm.SVC(kernel='rbf')
        max_feat = all_features[max_index]
        opt_theta.fit(max_feat, labels)
        p = opt_theta.predict(P_X)
        predictions.append(p)
    return(predictions)

def MV(predictions):
    aggregate = []
    p = predictions[0]
    for i in range(len(p)):
        aggregate.append([0, 0])

    for i in range(len(predictions)):
        for j in range(len(p)):
            if predictions[i][j] == 0:
                aggregate[j][0] += 1
            else:
                aggregate[j][1] += 1

    mv_prediction = []
    for i in range(len(p)):
        if aggregate[i][0] >= aggregate[i][1]:
            mv_prediction.append(0)
        else:
            mv_prediction.append(1)
    return mv_prediction

mv_prediction = MV(SESVM(X, Y, 9, P_X, P_Y))
print("Accuracy: " + str(accuracy_score(P_Y, mv_prediction)))