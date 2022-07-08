#This script tests the Random Forest on the entire preprocessed data set

import CombineLigandsProteins
import FixedClassificationModel
import AdjustingThreshold
#import FeatureImportance

CombineLigandsProteins.import_final()
testX = CombineLigandsProteins.X
testY = CombineLigandsProteins.Y

AdjustingThreshold.train(testX, testY)

#FeatureImportance.train(testX, testY, CombineLigandsProteins.feats)
FixedClassificationModel.train(testX, testY)

"""
accuracy = 0
recall = 0
for i in range(50):
    print("run " + str(i))
    acc,rec = FixedClassificationModel.train(testX, testY)
    accuracy += acc
    recall += rec
print('Average Accuracy: ' + str(accuracy/50))
print('Average Recall: ' + str(recall/50))
"""