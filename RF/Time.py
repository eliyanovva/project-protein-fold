import time
import CombineLigandsProteins
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn import metrics

with open('Time.csv', 'w') as f:
    print('process,train,predict', file = f)
    totprocess = 0
    tottrain = 0
    totpredict = 0
    for i in range(1000):
        print(i)
        #Time to process the data
        start = time.time()

        CombineLigandsProteins.import_final()
        testX = CombineLigandsProteins.X
        testY = CombineLigandsProteins.Y

        X_train, X_test, y_train, y_test = train_test_split(testX, testY, stratify=testY, test_size=0.1) 

        end = time.time()

        process =  end-start

        #Time to train the model
        start = time.time()

        clf=RandomForestClassifier(n_estimators=100, class_weight="balanced") #add in , class_weight="balanced"
        clf.fit(X_train,y_train)

        end = time.time()

        train = end-start

        #Time to make predictions
        start = time.time()

        y_pred=clf.predict_proba(X_test)[:,1]

        end = time.time()

        predict = end-start

        print(str(process) + ',' + str(train) + ',' + str(predict), file = f)
        totprocess += process
        tottrain += train
        totpredict += predict

    avprocess = totprocess / 1000
    avtrain = tottrain / 1000
    avpredict = totpredict / 1000
    print(avprocess)
    print(avtrain)
    print(avpredict)
    print(str(avprocess)+','+str(avtrain)+','+str(avpredict), file = f)