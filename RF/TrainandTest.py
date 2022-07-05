#This script tests the Random Forest on the entire preprocessed data set

import CombineLigandsProteins
import FixedClassificationModel
#import FeatureImportance

CombineLigandsProteins.import_final()
testX = CombineLigandsProteins.X
testY = CombineLigandsProteins.Y

#FeatureImportance.train(testX, testY, CombineLigandsProteins.feats) 
print("imported matrices")

FixedClassificationModel.train(testX, testY)

"""
accuracy = 0
recall = 0

with open("scores.txt", "w") as f:
    for i in range(10):
        print("run " + str(i), file=f)
        acc,rec = FixedClassificationModel.train(testX, testY)
        accuracy += acc
        recall += rec

    print('Average Accuracy: ' + str(accuracy/10), file=f)
    print('Average Recall: ' + str(recall/10), file=f)
"""