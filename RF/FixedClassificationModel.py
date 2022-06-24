#This script is a Classification Random Forest Model with Oversampling
#Code adapted from: https://www.datacamp.com/tutorial/random-forests-classifier-python 

#imports
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn import metrics
from imblearn.over_sampling import RandomOverSampler

def train(features, labels):
    #define features and labels
    X = features #Kmers
    y = labels #Binds or not

    #split into training and test set
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.1) # 90% training and 10% test

    #Oversampling was necessary, because most ligand/receptor pairs do not bind in our dataset
    ros = RandomOverSampler(random_state = 42)

    X_res, y_res = ros.fit_resample(X_train, y_train)

    #Create a Gaussian Regression
    clf=RandomForestClassifier(n_estimators=100)

    #Train the model
    clf.fit(X_res,y_res)

    #Form predictions
    y_pred=clf.predict(X_test)

    #Print accuracy of the model
    print("Accuracy:",metrics.roc_auc_score(y_test, y_pred))

