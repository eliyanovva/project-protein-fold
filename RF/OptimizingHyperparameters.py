#This script performs cross validation to optimize the hyperparameters of the Random Forest model. 

#Adapted from: https://www.youtube.com/watch?v=SctFnD_puQI

#Imports
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.experimental import enable_halving_search_cv
from sklearn.model_selection import HalvingRandomSearchCV
from imblearn.under_sampling import RandomUnderSampler
from sklearn.model_selection import train_test_split
import CombineLigandsProteins

#Number of trees in the forest
n_estimators = [int(x) for x in np.linspace(start = 10, stop = 1000, num = 100)]
#Features considered at each split
max_features = ['log2', 'sqrt', None]
#Number of levels in each tree
max_depth = [None]
#Samples required to split a node
min_samples_split = [int(x) for x in np.linspace(start = 10, stop = 1000, num = 100)]
#Minimum samples at each leaf
min_samples_leaf = [int(x) for x in np.linspace(start = 1, stop = 100, num = 100)]
#Method of selecting samples for each tree
bootstrap = [True]
#Best-first
max_leaf_nodes = [None]
#Should it split?
min_impurity_decrease = [int(x) for x in np.linspace(start = 0, stop = 5, num = 5)]
#Estimate generalization score
oob_score = [True, False]
#Samples
max_samples = [None]

param_grid = {'n_estimators': n_estimators,
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

#Train the model (performing cross validation)
grid.fit(X_res, y_res)
#Print optimized parameters
print(grid.best_params_)

#Test the model on the test set
print(grid.score(X_test, y_test))