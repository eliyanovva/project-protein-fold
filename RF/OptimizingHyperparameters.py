#Adapted from: https://www.youtube.com/watch?v=SctFnD_puQI

#Imports
import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.experimental import enable_halving_search_cv
from sklearn.model_selection import HalvingRandomSearchCV
from imblearn.under_sampling import RandomUnderSampler
from sklearn.model_selection import train_test_split
import CombineLigandsProteins

#Number of trees in the forest
n_estimators = [np.linspace(60, 80, 1)]
#Criterion
criterion = ['gini', 'entropy', 'log_loss']
#Features considered at each split
max_features = ['log2', 'sqrt', None]
#Number of levels in each tree
max_depth = [np.linspace(10, 20, 5), None]
#Samples required to split a node
min_samples_split = [np.linspace(2, 70, 2)]
#Minimum samples at each leaf
min_samples_leaf = [np.linspace(1, 5, 1)]
#Method of selecting samples for each tree
bootstrap = [True]
#Best-first
max_leaf_nodes = [np.linspace(10, 100, 10), None]
#Should it split?
min_impurity_decrease = [np.linspace(0, 5, .5)]
#Estimate generalization score
oob_score = [True, False]
#Samples
max_samples = [np.linspace(.5, 5, .5), None]

param_grid = {'n_estimators': n_estimators,
                'criterion' : criterion,
                'max_features' : max_features,
                'max_depth' : max_depth,
                'min_samples_split' : min_samples_split,
                'min_samples_leaf' : min_samples_leaf,
                'bootstrap' : bootstrap,
                'max_leaf_nodes' : max_leaf_nodes,
                'min_impurity_decrease' : min_impurity_decrease,
                'oob_score' : oob_score,
                'max_samples' : max_samples
                }

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
ros = RandomUnderSampler(random_state = 42)
X_res, y_res = ros.fit_resample(X_train, y_train)

grid.fit(X_res, y_res)
print(grid.best_params_)

print(grid.score(X_test, y_test))