#This script is a Classification Random Forest Model with Cost Sensitive Learning
#Code adapted from: https://www.datacamp.com/tutorial/random-forests-classifier-python 

#imports
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn import metrics

def train(features, labels, BALANCE):
    #define features and labels
    X = features #Kmers
    y = labels #Binds or not

    #Create a Gaussian Regression
    #For balanced datasets
    if BALANCE == True:
        # split into training and test set
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.1)  # 90% training and 10% test
        clf=RandomForestClassifier(n_estimators=100)
        #for unbalanced datasets
    elif BALANCE == False:
        # split into training and test set
        X_train, X_test, y_train, y_test = train_test_split(X, y, stratify=y, test_size=0.1)  # 90% training and 10% test
        clf = RandomForestClassifier(n_estimators=100, class_weight="balanced")  # cost-sensitive learning
    #Train the model
    clf.fit(X_train,y_train)

    #Form predictions
    y_pred=clf.predict_proba(X_test)[:,1]
    precision, recall, thresholds = metrics.precision_recall_curve(y_test, y_pred)

    acc = metrics.roc_auc_score(y_test, y_pred)
    rec = metrics.auc(recall,precision)
    #Print accuracy of the model
    #print("Accuracy:",acc)
    #print("Recall:",rec)

    y_pred=clf.predict(X_test)

    mat = (metrics.matthews_corrcoef(y_test, y_pred))
    bac = metrics.balanced_accuracy_score(y_test, y_pred)
    #print('Balanced Accuracy Score: ' + str(bac))

    TN, FN, TP, FP = matthew_counts(y_test, y_pred)

    return acc,rec,bac,mat,TN, FN, TP, FP

#Examine TP and TN rates
def matthew_counts(y_test, y_pred):
    TN = 0
    FN = 0
    TP = 0
    FP = 0

    for i in range(len(y_test)):
        if (y_test[i] == 0) & (y_pred[i] == 0):
            TN += 1
        if (y_test[i] == 1) & (y_pred[i] == 0):
            FN += 1
        if (y_test[i] == 1) & (y_pred[i] == 1):
            TP += 1
        if (y_test[i] == 0) & (y_pred[i] == 1):
            FP += 1
    return TN, FN, TP, FP

#For training the model on new pairs
def train_new_pairs(features, labels, new_features, balance):
    X = features
    Y = labels
    X_n = new_features

    if balance == False:
        clf = RandomForestClassifier(n_estimators=100, class_weight="balanced")  # add in , class_weight="balanced"
        X_train, X_test, y_train, y_test = train_test_split(X, Y, stratify=Y, test_size=0.1)  # 90% training and 10% test
    else:
        clf = RandomForestClassifier(n_estimators=100)
        X_train, X_test, y_train, y_test = train_test_split(X, Y, test_size=0.1)  # 90% training and 10% test
    # Train the model
    clf.fit(X_train, y_train)

    n_pred = clf.predict(X_n)

    return n_pred




