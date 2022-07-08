#This script determines the best probability threshold for determining positive vs negative pairings
#Adapted from: https://machinelearningmastery.com/threshold-moving-for-imbalanced-classification/

#imports
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn import metrics
from sklearn.metrics import precision_recall_curve
from numpy import nanargmax
from imblearn.over_sampling import RandomOverSampler

def train(features, labels):
    #define features and labels
    X = features #Kmers
    y = labels #Binds or not

    #split into training and test set
    X_train, X_test, y_train, y_test = train_test_split(X, y, stratify=y,test_size=0.1) # 90% training and 10% test
    print('split data')

    #compare to random undersampling

    #Create a Gaussian Regression
    clf=RandomForestClassifier(n_estimators=100, class_weight="balanced") #add in , class_weight="balanced"
    print('made classifier')
    #Train the model
    clf.fit(X_train,y_train)
    print('fit the data')

    #Train the model
    clf.fit(X_train,y_train)

    #Form predictions
    y_pred=clf.predict_proba(X_test)

    y_pred = y_pred[:,1]

    precision, recall, thresholds = precision_recall_curve(y_test, y_pred)

    fscore = (2 * precision * recall) / (precision + recall)

    ix = nanargmax(fscore)
    print('Best Threshold=%f, F-Score=%.3f' % (thresholds[ix], fscore[ix]))
