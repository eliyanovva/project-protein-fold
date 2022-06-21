#Adapted from: https://www.youtube.com/watch?v=SctFnD_puQI

#Imports
import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import RandomizedSearchCV

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

model = RandomForestClassifier()

grid = RandomizedSearchCV(estimator = model, param_distributions = param_grid, cv = 10, verbose = 2, n_jobs = -1)
