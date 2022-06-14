#Code informed by: https://machinelearningmastery.com/random-forest-ensemble-in-python/ 

#imports
from numpy import mean
from numpy import std
from sklearn.datasets import make_regression
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import RepeatedKFold
from sklearn.ensemble import RandomForestRegressor
import CombineLigandsProteins
from sklearn.model_selection import train_test_split
from sklearn import metrics

def train(features, labels):
    #define features and labels
    X = features #Globals.features (kmers)
    y = labels #logFC

    #split into training and test set
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3) # 70% training and 30% test

    model = RandomForestRegressor() #For regression RF

    #Perform cross validation
    crossval = RepeatedKFold(n_splits=10, n_repeats=3, random_state=1) 
    n_scores = cross_val_score(model, X_train, y_train, scoring='neg_mean_absolute_error', cv=crossval, n_jobs=-1, error_score='raise')
    #Scoring the model
    print('MAE: %.3f (%.3f)' % (mean(n_scores), std(n_scores)))

    model.fit(X_train, y_train)
    #Form predictions
    y_pred=model.predict(X_test)
    print(y_pred)
    print(y_test)

    #Print accuracy of the model
    print("Accuracy:",metrics.r2_score(y_test, y_pred))

CombineLigandsProteins.import_final()
testX = CombineLigandsProteins.X
testY = CombineLigandsProteins.Y
train(testX, testY)





