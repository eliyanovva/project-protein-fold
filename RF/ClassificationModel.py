#Code adapted from: https://www.datacamp.com/tutorial/random-forests-classifier-python 

#imports
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn import metrics
import Resample
Resample.resampled_matrices()
testX = Resample.input
testY = Resample.output

def train(features, labels):
    #define features and labels
    X = features #Globals.features (kmers)
    y = labels #logFC

    #split into training and test set
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.1) # 90% training and 10% test

    #Create a Gaussian Regression
    clf=RandomForestClassifier(n_estimators=100)

    #Train the model
    clf.fit(X_train,y_train)

    #Form predictions
    y_pred=clf.predict(X_test)

    #Print accuracy of the model
    #print("Accuracy:",metrics.roc_auc_score(y_test, y_pred))

train(testX, testY)