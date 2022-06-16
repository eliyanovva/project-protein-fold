#Code adapted from: https://www.datacamp.com/tutorial/random-forests-classifier-python 

#imports
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn import metrics
import CombineLigandsProteins

def train(features, labels):
    #define features and labels
    X = features #Globals.features (kmers)
    y = labels #logFC

    #split into training and test set
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.1) # 90% training and 10% test

    #Create a Gaussian Regression
    clf=RandomForestRegressor(n_estimators=100)

    #Train the model
    clf.fit(X_train,y_train)

    #Form predictions
    y_pred=clf.predict(X_test)
    print(y_pred)
    print(y_test)

    #Print accuracy of the model
    print("Accuracy:",metrics.r2_score(y_test, y_pred))


CombineLigandsProteins.import_final()
testX = CombineLigandsProteins.X
testY = CombineLigandsProteins.Y
#train(testX, testY)