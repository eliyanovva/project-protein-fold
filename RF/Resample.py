#Adapted from: https://www.geeksforgeeks.org/introduction-to-resampling-methods/
# importing libraries
import pandas as pd
import numpy as np
import seaborn as sns
from sklearn.preprocessing import StandardScaler
from imblearn.under_sampling import RandomUnderSampler, TomekLinks
from imblearn.over_sampling import RandomOverSampler, SMOTE
import CombineLigandsProteins

ros = RandomOverSampler(random_state = 42)

CombineLigandsProteins.import_final()
testX = CombineLigandsProteins.X
testY = CombineLigandsProteins.Y


X_res, y_res = ros.fit_resample(testX, testY)
   
X_res = pd.DataFrame(X_res)
Y_res = pd.DataFrame(y_res)
   
   
print("After Over Sampling Of Minor Class Total Samples are :", len(Y_res))
print('Class 0        :', round(Y_res[0].value_counts()[0]
                /len(Y_res) * 100, 2), '% of the dataset')
   
print('Class 1:', round(Y_res[0].value_counts()[1]
                /len(Y_res) * 100, 2), '% of the dataset')