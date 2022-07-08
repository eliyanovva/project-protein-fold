from sklearn.model_selection import learning_curve
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn import metrics
import CombineLigandsProteins
import numpy as np
import matplotlib.pyplot as plt

CombineLigandsProteins.import_final()
X = CombineLigandsProteins.X #Kmers
y = CombineLigandsProteins.Y #Binds or not

#split into training and test set
X_train, X_test, y_train, y_test = train_test_split(X, y, stratify=y,test_size=0.1) # 90% training and 10% test

train_sizes, train_scores, valid_scores = learning_curve(RandomForestClassifier(n_estimators=100, class_weight="balanced"), X_train, y_train, train_sizes=np.arange(50, 20372, 100), cv=5)

plt.plot(train_sizes, train_scores)
plt.plot(train_sizes, valid_scores)
plt.savefig("Learning_Curve.png")

