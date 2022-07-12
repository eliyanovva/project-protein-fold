#This script is a Classification Random Forest Model with Undersampling
#Code adapted from: https://www.datacamp.com/tutorial/random-forests-classifier-python 

#imports
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn import metrics

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

    #Form predictions
    y_pred=clf.predict_proba(X_test)[:,1]
    print('made predictions')
    precision, recall, thresholds = metrics.precision_recall_curve(y_test, y_pred)

    acc = metrics.roc_auc_score(y_test, y_pred)
    rec = metrics.auc(recall,precision)
    #Print accuracy of the model
    print("Accuracy:",acc)
    print("Recall:",rec)

    y_pred=clf.predict(X_test)
    print(metrics.matthews_corrcoef(y_test, y_pred))
    bac = metrics.balanced_accuracy_score(y_test, y_pred)
    print('Balanced Accuracy Score: ' + str(bac))

    return acc,rec,bac

