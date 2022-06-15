#Code informed by: https://machinelearningmastery.com/random-forest-ensemble-in-python/ 
#and https://scikit-learn.org/stable/modules/grid_search.html#grid-search
# and https://towardsdatascience.com/faster-hyperparameter-tuning-with-scikit-learn-71aa76d06f12 

#imports

from sklearn.ensemble import RandomForestRegressor
import CombineLigandsProteins
from sklearn.model_selection import train_test_split
from sklearn import metrics
from sklearn.experimental import enable_halving_search_cv
from sklearn.model_selection import HalvingGridSearchCV
import pandas as pd
from sklearn.metrics import mean_squared_log_error, make_scorer
import numpy as np

def train(features, labels):
    #define features and labels
    X = features #Globals.features (kmers)
    y = labels #logFC

    #split into training and test set
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3) # 70% training and 30% test

    model = RandomForestRegressor() #For regression RF

    print(model.get_params().keys())

    #Cross validate
    FACTOR = 5
    MAX_RESOURCE_DIVISOR = 4
    # Number of trees in random forest
    n_estimators = [int(x) for x in np.linspace(start = 200, stop = 2000, num = 10)]
    # Number of features to consider at every split
    max_features = ['auto', 'sqrt']
    # Maximum number of levels in tree
    max_depth = [int(x) for x in np.linspace(10, 110, num = 11)]
    max_depth.append(None)
    # Minimum number of samples required to split a node
    min_samples_split = [2, 5, 10]
    # Minimum number of samples required at each leaf node
    min_samples_leaf = [1, 2, 4]
    # Method of selecting samples for training each tree
    bootstrap = [True, False]
    # Create the random grid
    param_grid = {'n_estimators': n_estimators,
                'max_features': max_features,
                'max_depth': max_depth,
                'min_samples_split': min_samples_split,
                'min_samples_leaf': min_samples_leaf,
                'bootstrap': bootstrap}
 
    rmsle = lambda y_true, y_pred: metrics.r2_score(y_true, y_pred)
    scorer = make_scorer(rmsle, greater_is_better=False)
    
    grid_search_params = dict(estimator=model, param_grid=param_grid,
                          scoring=scorer,
                          cv=3,
                          n_jobs=-1,
                          verbose=2)

    n_samples = len(X_train)
    halving_results_n_samples = HalvingGridSearchCV(resource='n_samples', 
                                min_resources=n_samples//MAX_RESOURCE_DIVISOR,
                                factor=FACTOR,
                                **grid_search_params
                                ).fit(X_train, y_train)
    print(halving_results_n_samples.best_params_)

CombineLigandsProteins.import_final()
testX = CombineLigandsProteins.X
testY = CombineLigandsProteins.Y
train(testX, testY)





