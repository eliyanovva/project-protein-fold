#This script tests the Random Forest on the entire preprocessed data set
import CombineLigandsProteins
#import Sequence_only
import FixedClassificationModel
#import Structure_only
#import AdjustingThreshold
#import Feature_Importance.FeatureImportance as fi
import Metrics_Graphs

CombineLigandsProteins.import_final()
testX = CombineLigandsProteins.X
testY = CombineLigandsProteins.Y

#Sequence_only.import_final()
#seqX = Sequence_only.X
#seqY = Sequence_only.Y

#Structure_only.import_final()
#structX = Structure_only.X
#structY = Structure_only.Y

#AdjustingThreshold.train(testX, testY)

#fi.importance_file(testX, testY, CombineLigandsProteins.feats)
#FixedClassificationModel.train(testX, testY)
Metrics_Graphs.train(testX, testY)

"""
accuracy = 0
recall = 0
BAC = 0
MAT = 0

f = open('results_true_false_vals_struct.csv', 'w')
f.write('Run,TN,FN,TP,FP' + "\n")

for i in range(50):
    print("run " + str(i))
    acc, rec, bac, mat, TN, FN, TP, FP = FixedClassificationModel.train(structX, structY)
    accuracy += acc
    recall += rec
    BAC += bac
    MAT += mat
    f.write(str(i+1) + ", " + str(TN) + ", " + str(FN) + ", " + str(TP) + ", " + str(FP) + "\n")

print('Average Accuracy: ' + str(accuracy/50))
print('Average Recall: ' + str(recall/50))
print('Average Balanced: ' + str(BAC/50))
print('Average Matthew: ' + str(mat/50))

f.close()

"""
