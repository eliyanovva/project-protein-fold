#This script is a Classification Random Forest Model with Undersampling
#Code adapted from: https://www.datacamp.com/tutorial/random-forests-classifier-python 

#imports
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn import metrics
from imblearn.under_sampling import InstanceHardnessThreshold
import numpy as np

def train(features, labels, protein_freqs, ligand_freqs):
    #define features and labels
    X = features #Kmers
    y = labels #Binds or not

    #split into training and test set
    X_train, X_test, y_train, y_test = train_test_split(X, y, stratify=y,test_size=0.1) # 90% training and 10% test
    print('split data')
    #Undersampling was necessary, because most ligand/receptor pairs do not bind in our dataset
    ih = InstanceHardnessThreshold(n_jobs=3, cv=3)
    print('assigned iht')
    X_res, y_res = ih.fit_resample(np.int_(X_train), np.int_(y_train))
    print('iht sampling')
    print(len(X_res))
    #obs in X_train: 28296
    #obs in X_res: 23523 (?)
    """
    f = open('pl_pairs3.txt', "w")

    i = 0

    X_res_str = []
    X_train_str = []

    k = 0
    for obs in X_train:
        print(k)
        k += 1
        string = ""
        for freq in obs:
            string += str(freq)
        X_train_str.append(string)

    k = 0
    for obs in X_res:
        print(k)
        k += 1
        string = ""
        for freq in obs:
            string += str(freq)
        X_res_str.append(string)

    j = 0
    for i in range(len(X_train)):
        print(j)
        string = ""
        for freq in X_train[i]:
            string += str(freq)

        if X_res_str.count(string) == 0:

            TM3_AA = X_res[i][:904]
            TM5_AA = X_res[i][904:1906]
            TM6_AA = X_res[i][1906:2629]
            TM7_AA = X_res[i][2629:3607]
            TM3_Di = X_res[i][3607:3890]
            TM5_Di = X_res[i][3890:4471]
            TM6_Di = obs[4471:5300]
            TM7_Di = obs[5300:6014]
            ligand = X_res[i][6014:]

            print('NOT SAMPLED ' + str(i) + "\n")
            i+= 1
            p = ""
            l = ""

            #https://www.geeksforgeeks.org/numpy-array_equal-python/
            for prot in protein_freqs:
                if ((np.array_equal(np.array(TM3_AA), np.array(protein_freqs[prot][0]))) &
                    (np.array_equal(np.array(TM5_AA), np.array(protein_freqs[prot][1]))) &
                    (np.array_equal(np.array(TM6_AA), np.array(protein_freqs[prot][2]))) &
                    (np.array_equal(np.array(TM7_AA), np.array(protein_freqs[prot][3])))):
                        p = prot
                        break

            for lig in ligand_freqs:
                if (np.array_equal(np.array(ligand_freqs[lig]), np.array(ligand))):
                    l = lig
                    break

            f.write(p + " " + l + "\n")

    f.close()
    """

    #Create a Gaussian Regression
    clf=RandomForestClassifier(n_estimators=100)
    print('made classifier')
    #Train the model
    clf.fit(X_res,y_res)
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

    return acc,rec

