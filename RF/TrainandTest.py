#This script tests the Random Forest on the entire preprocessed data set

import CombineLigandsProteins
import FixedClassificationModel
#import AdjustingThreshold
#import Feature_Importance.FeatureImportance as fi

CombineLigandsProteins.import_final()
testX = CombineLigandsProteins.X
testY = CombineLigandsProteins.Y

#AdjustingThreshold.train(testX, testY)

#fi.train(testX, testY, CombineLigandsProteins.feats)
FixedClassificationModel.train(testX, testY)


accuracy = 0
recall = 0
BAC = 0
#f = open('results_filter10.csv', 'w')
#f.write('Run, Accuracy, Recall, Balanced Score' + "\n")

for i in range(50):
    print("run " + str(i))
    acc, rec, bac = FixedClassificationModel.train(testX, testY)
    accuracy += acc
    recall += rec
    BAC += bac
    #f.write(str(i+1) + ", " + str(acc) + ", " + str(rec) + ", " + str(bac) + "\n")

print('Average Accuracy: ' + str(accuracy/50))
print('Average Recall: ' + str(recall/50))
print('Average Balanced: ' + str(BAC/50))
#f.close()
