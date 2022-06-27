#This script implements a Selective Ensemble Support Vector Machine algorithm
# based on descriptions from https://dl.acm.org/doi/pdf/10.1145/3307339.3342141
#random documentation from https://www.geeksforgeeks.org/random-numbers-in-python/

from sklearn import svm
from sklearn.metrics import accuracy_score
from sklearn.model_selection import train_test_split
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
    (assume N- can be partitioned into M sets without overlapping)
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

CombineLigandsProteins.import_final()

X = CombineLigandsProteins.final_matrix
Y = CombineLigandsProteins.logFCmat

print('Started the SESVM')

N, P, Y_n, Y_p = train_test_split(X, Y, test_size=.1)

# separate a test set into positive and negative observations
def seperate_sets(N, Y):
    pos_set = []
    neg_set = []

    for i in range(len(Y)):
        if Y[i] == 0:           #indicates a negative label
            #neg_set.append(list(N[i]))
            neg_set.append(i)
        else:
            #pos_set.append(list(N[i]))
            pos_set.append(i)
    return pos_set, neg_set

#partition the negative observations into M sets
#each paritioned set should be the same size of the set of positive observations
def create_partitions(pos_set, neg_set, M):
    partitions = []
    part_len = len(pos_set)
    print(len(pos_set))
    print(len(neg_set))
    print(M)
    overhang = part_len * M - len(neg_set)
    print(overhang)
    duplicate = []
    for item in neg_set:
        duplicate.append(item)
    print('made duplicate')
    for i in range(int(overhang)):
        dup_row = random.choice(duplicate)
        neg_set.append(dup_row)
        duplicate.remove(dup_row)
    print('completed the overhang')
    print(str(M))
    print(str(part_len))
    k = 0
    for i in range(int(M)):
        n_set = []
        for j in range(part_len):
            row = random.choice(neg_set)
            n_set.append(row)
            neg_set.remove(row)
            k += 1
            if k % 10 == 0:
                print('chose row ' + str(k))
        partitions.append(n_set)
        print('made a partition')
    return partitions

def SESVM(N, Y, T, P_X, P_Y):
    pos_set, neg_set = seperate_sets(N, Y)

    print('seperated the sets')

    M = np.ceil(float(len(neg_set)) / float(len(pos_set)))
    predictions = []

    for i in range(T):
        neg_copy = []
        for item in neg_set:
            neg_copy.append(item)
        print('made neg copy')
        parts = create_partitions(pos_set, neg_copy, M)
        print('made the partitions')
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

        #out of the M training sets, find the most accurate classifier
        max_accuracy = 0
        max_index = 0
        for m in range(int(M)):
            if accuracies[m] > max_accuracy:
                max_accuracy = accuracies[m]
                max_index = m
        #retrain a classifer based on the optimal training set, and use it to predict against the test set P
        opt_theta = svm.SVC(kernel='rbf')
        max_feat = all_features[max_index]
        opt_theta.fit(max_feat, labels)
        p = opt_theta.predict(P_X)
        predictions.append(p)
    return(predictions)

#Majority Voting algorithm
#Out of T predictions vectors, it will select the most frequent prediction for each sample
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

mv_prediction = MV(SESVM(N, Y_n, 3, P, Y_p))
print("Accuracy: " + str(accuracy_score(Y_p, mv_prediction)))

def MV_weighted(accuracies):
    print()