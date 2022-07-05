#This script tests the Random Forest on the entire preprocessed data set

import CombineLigandsProteins
import FixedClassificationModel

CombineLigandsProteins.import_final()
testX = CombineLigandsProteins.X
testY = CombineLigandsProteins.Y

accuracy = 0
recall = 0

#acc,rec = FixedClassificationModel.train(testX, testY)


for i in range(5):
    print("run " + str(i))
    acc,rec = FixedClassificationModel.train(testX, testY)
    accuracy += acc
    recall += rec

print('Average Accuracy: ' + str(accuracy/5))
print('Average Recall: ' + str(recall/5))
