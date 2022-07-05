#This script extracts feature importance
#Code adapted from: https://www.datacamp.com/tutorial/random-forests-classifier-python 
#and https://scikit-learn.org/stable/auto_examples/ensemble/plot_forest_importances.html

#imports
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from imblearn.under_sampling import InstanceHardnessThreshold
import numpy as np

def train(features, labels, filter_feat):
    #define features and labels
    X = features #Kmers
    y = labels #Binds or not

    #split into training and test set
    X_train, X_test, y_train, y_test = train_test_split(X, y, stratify=y,test_size=0.1) # 90% training and 10% test

    #Oversampling was necessary, because most ligand/receptor pairs do not bind in our dataset
    ih = InstanceHardnessThreshold(n_jobs=-1, cv=3)

    X_res, y_res = ih.fit_resample(np.int_(X_train), np.int_(y_train))

    #Create a Gaussian Regression
    clf=RandomForestClassifier(n_estimators=100)

    #Train the model
    clf.fit(X_res,y_res)

    #Finding most important features
    importances = clf.feature_importances_

    #Writes the features of the model in order of importance to the important_features.txt file
    with open("important_features.txt", "w") as f:
        while importances.size > 0:
            index = np.argmax(importances)
            print(filter_feat[index], file=f)
            print(importances[index], file=f)
            importances = np.delete(importances, index)