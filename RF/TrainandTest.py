#This script tests the Random Forest on the entire preprocessed data set

import CombineLigandsProteins
import FixedClassificationModel
#import FeatureImportance

CombineLigandsProteins.import_final()
testX = CombineLigandsProteins.X
testY = CombineLigandsProteins.Y

print("imported matrices")

#FeatureImportance.train(testX, testY, CombineLigandsProteins.feats)
#FixedClassificationModel.train(testX, testY)


accuracy = []
recall = []

f = open('results_100.csv', "w")

for i in range(13):
    print("run " + str(i))
    acc,rec = FixedClassificationModel.train(testX, testY)
    accuracy.append(acc)
    recall.append(rec)
    f.write(str(i+1) + "," + str(acc) + "," + str(rec) + "\n")

print(accuracy)
print(recall)
#print('Average Accuracy: ' + str(accuracy/5))
#print('Average Recall: ' + str(recall/5))

f.close()