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
    ih = InstanceHardnessThreshold(n_jobs=4, cv=3)
    print('assigned iht')
    X_res, y_res = ih.fit_resample(np.int_(X_train), np.int_(y_train))
    print('iht sampling')

    #obs in X_train: 39402
    #obs in X_res: 9248

    #f = open('pl_pairs.txt', "w")

    for obs in X_res:
        TM3_AA = obs[:904]
        TM5_AA = obs[904:1906]
        TM6_AA = obs[1906:2629]
        TM7_AA = obs[2629:3607]
        TM3_Di = obs[3607:3890]
        TM5_Di = obs[3890:4471]
        TM6_Di = obs[4471:5300]
        TM7_Di = obs[5300:6014]
        ligand = obs[6014:]

        print(len(TM3_AA))
        print(len(TM5_AA))
        print(len(TM6_AA))
        print(len(TM7_AA))
        print(len(TM3_Di))
        print(len(TM5_Di))
        print(len(TM6_Di))
        print(len(TM7_Di))
        print(len(ligand))

        p = ""
        l = ""

        #https://www.geeksforgeeks.org/numpy-array_equal-python/
        for prot in protein_freqs:
            if ((np.array_equal(np.array(TM3_AA), np.array(protein_freqs[prot][0]))) &
                (np.array_equal(np.array(TM5_AA), np.array(protein_freqs[prot][1]))) &
                (np.array_equal(np.array(TM6_AA), np.array(protein_freqs[prot][2]))) &
                (np.array_equal(np.array(TM7_AA), np.array(protein_freqs[prot][3]))) &
                (np.array_equal(np.array(TM3_Di), np.array(protein_freqs[prot][4]))) &
                (np.array_equal(np.array(TM5_Di), np.array(protein_freqs[prot][5]))) &
                (np.array_equal(np.array(TM6_Di), np.array(protein_freqs[prot][6]))) &
                (np.array_equal(np.array(TM7_Di), np.array(protein_freqs[prot][7])))):
                    p = prot
                    break

        for lig in ligand_freqs:
            if (np.array_equal(np.array(ligand_freqs[lig]), np.array(ligand))):
                l = lig
                break

        #f.write(p + " " + l + "\n")

    #f.close()

    print(len(X_train))
    print(len(X_train[0]))
    print(len(X_res))
    print(len(X_res[0]))

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

