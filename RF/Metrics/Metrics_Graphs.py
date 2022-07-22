#This script creates graphs of the Precision-Recall and Receiver Operating Characteristic curves
#Code adapted from: https://www.datacamp.com/tutorial/random-forests-classifier-python 

#imports
import statistics
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn import metrics
import matplotlib.pyplot as plt

def train(features, labels):
    #define features and labels
    X = features #Kmers
    y = labels #Binds or not
    
    rp = {}
    roc = {}
    for i in range(50):
        #split into training and test set
        X_train, X_test, y_train, y_test = train_test_split(X, y, stratify=y,test_size=0.1) # 90% training and 10% test

        #compare to random undersampling

        #Create a Gaussian Regression
        clf=RandomForestClassifier(n_estimators=100, class_weight="balanced") #add in , class_weight="balanced"

        #Train the model
        clf.fit(X_train,y_train)

        #Form predictions
        y_pred=clf.predict_proba(X_test)[:,1]
        precision, recall, thresholds = metrics.precision_recall_curve(y_test, y_pred)
        fpr, tpr, thresholds = metrics.roc_curve(y_test, y_pred)
        print(precision)
        #Collect values for graphs
        for j in range(len(recall)):
            if recall[j] not in rp:
                rp[recall[j]] = []
                rp[recall[j]].append(precision[j])
        for j in range(len(fpr)):
            if fpr[j] not in roc:
                roc[fpr[j]] = []
            roc[fpr[j]].append(tpr[j])

    #Take the average of collected values
    for num in list(rp.keys()):
        rp[num] = statistics.mean(rp[num])
    for num in list(roc.keys()):
        roc[num] = statistics.mean(roc[num])

    print(rp.values())
    
    #Precision Recall
    fig1 = plt.plot(list(rp.keys()), list(rp.values()))
    plt.xlabel("Recall")
    plt.ylabel("Precision")
    plt.title("Precision-Recall Curve")
    plt.savefig("Precision_Recall.png")

    #Receiver Operating Characteristic
    plt.clf()
    fig2 = plt.plot(list(roc.keys()), list(roc.values()))
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.title("ROC Curve")
    plt.savefig("ROC.png")

    