#This script is a Classification Random Forest Model with Undersampling
#Code adapted from: https://www.datacamp.com/tutorial/random-forests-classifier-python 

#imports
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn import metrics
from imblearn.under_sampling import InstanceHardnessThreshold
import numpy as np

def train(features, labels):
    #define features and labels
    X = features #Kmers
    y = labels #Binds or not

    #split into training and test set
    X_train, X_test, y_train, y_test = train_test_split(X, y, stratify=y,test_size=0.1) # 90% training and 10% test
  
    #Undersampling was necessary, because most ligand/receptor pairs do not bind in our dataset
    ih = InstanceHardnessThreshold(n_jobs=-1, cv=3)
    
    X_res, y_res = ih.fit_resample(np.int_(X_train), np.int_(y_train))
  

    #Create a Gaussian Regression
    clf=RandomForestClassifier(n_estimators=100)

    #Train the model
    clf.fit(X_res,y_res)

    #Form predictions
    y_pred=clf.predict_proba(X_test)[:,1]

    precision, recall, thresholds = metrics.precision_recall_curve(y_test, y_pred)

    acc = metrics.roc_auc_score(y_test, y_pred)
    rec = metrics.auc(recall,precision)
    #Print accuracy of the model
    print("Accuracy:",acc)
    print("Recall:",rec)

    return acc,rec

