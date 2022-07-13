#This script performs cross validation to optimize the hyperparameters of the Random Forest model. 

#Adapted from: https://www.youtube.com/watch?v=SctFnD_puQI and https://datascience.stackexchange.com/questions/44327/oversampling-before-cross-validation-is-it-a-problem


#Imports
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.experimental import enable_halving_search_cv
from sklearn.model_selection import HalvingRandomSearchCV
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import train_test_split
import CombineLigandsProteins
from sklearn import metrics
from sklearn.metrics import precision_recall_curve
from numpy import nanargmax

#Number of trees in the forest
n_estimators = [int(x) for x in np.linspace(start = 10, stop = 1000, num = 100)]
#Features considered at each split
max_features = ['log2', 'sqrt', None]
#Number of levels in each tree
max_depth = [None]
#Samples required to split a node
min_samples_split = [int(x) for x in np.linspace(start = 10, stop = 700, num = 100)]
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

#Import data
CombineLigandsProteins.import_final()
testX = CombineLigandsProteins.X
testY = CombineLigandsProteins.Y

#split into training and test set
X_train, X_test, y_train, y_test = train_test_split(testX, testY, stratify=testY, test_size=0.1) # 90% training and 10% test

#Define model
model = RandomForestClassifier(class_weight = "balanced")
#Define k-folds
kf = StratifiedKFold()

#Define random grid
grid = HalvingRandomSearchCV(estimator = model, param_distributions = param_grid, verbose = 2, n_jobs = -1,
                            cv=kf, return_train_score=True, scoring = "balanced_accuracy")

#Train the model (performing cross validation)
grid.fit(X_train, y_train)
#Print optimized parameters
print(grid.best_params_)

#Form predictions
y_pred=grid.predict_proba(X_test)[:,1]

precision, recall, thresholds = metrics.precision_recall_curve(y_test, y_pred)

acc = metrics.roc_auc_score(y_test, y_pred)
rec = metrics.auc(recall,precision)
#Print accuracy of the model
print("Accuracy:",acc)
print("Recall:",rec)


#Determine optimal threshold
precision, recall, thresholds = precision_recall_curve(y_test, y_pred)

fscore = (2 * precision * recall) / (precision + recall)

ix = nanargmax(fscore)
print('Best Threshold=%f, F-Score=%.3f' % (thresholds[ix], fscore[ix]))

#Run the model with the optimal threshold
y_pred = (grid.predict_proba(X_test)[:,1] >= thresholds[ix]).astype(bool)

bac = metrics.balanced_accuracy_score(y_test, y_pred, adjusted = "True")

print('Balanced Accuracy Score: ' + str(bac))
print(metrics.matthews_corrcoef(y_test, y_pred))