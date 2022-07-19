import time
import CombineLigandsProteins
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn import metrics

#Time to process the data
start = time.time()

CombineLigandsProteins.import_final()
testX = CombineLigandsProteins.X
testY = CombineLigandsProteins.Y

X_train, X_test, y_train, y_test = train_test_split(testX, testY, stratify=testY, test_size=0.1) 

end = time.time()

print("The time execution to process the data is:", end-start)

#Time to train the model
start = time.time()

clf=RandomForestClassifier(n_estimators=100, class_weight="balanced") #add in , class_weight="balanced"
clf.fit(X_train,y_train)

end = time.time()

print("The time execution to train the model is:", end-start)

#Time to make predictions
start = time.time()

y_pred=clf.predict_proba(X_test)[:,1]

end = time.time()

print("The time execution to predict is:", end-start)