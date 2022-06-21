#Adapted from: https://www.youtube.com/watch?v=SctFnD_puQI

#Imports
import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import HalvingRandomSearchCV
from imblearn.over_sampling import RandomOverSampler
from sklearn.model_selection import train_test_split
import CombineLigandsProteins

#Number of trees in the forest
n_estimators = [int(x) for x in np.linspace(start = 10, stop = 100, num = 10)]
#Features considered at each split
max_features = ['auto', 'sqrt']
#Number of levels in each tree
max_depth = [int(x) for x in np.linspace(start = 2, stop = 10, num = 5)]
#Samples required to split a node
min_samples_split = [int(x) for x in np.linspace(start = 10, stop = 100, num = 10)]
#Minimum samples at each leaf
min_samples_leaf = [1,2]
#Method of selecting samples for each tree
bootstrap = [True, False]

param_grid = {'n_estimators': n_estimators,
                'max_features' : max_features,
                'max_depth' : max_depth,
                'min_samples_split' : min_samples_split,
                'min_samples_leaf' : min_samples_leaf,
                'bootstrap' : bootstrap}

#Define model
model = RandomForestClassifier()
#Define random grid
grid = HalvingRandomSearchCV(estimator = model, param_distributions = param_grid, verbose = 2, n_jobs = -1)

#Import data
CombineLigandsProteins.import_final()
testX = CombineLigandsProteins.X
testY = CombineLigandsProteins.Y

#split into training and test set
X_train, X_test, y_train, y_test = train_test_split(testX, testY, test_size=0.1) # 90% training and 10% test

#Oversample
ros = RandomOverSampler(random_state = 42)
X_res, y_res = ros.fit_resample(X_train, y_train)

grid.fit(X_res, y_res)
print(grid.best_params_)

print(grid.score(X_test, y_test))