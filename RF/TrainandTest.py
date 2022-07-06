#This script tests the Random Forest on the entire preprocessed data set

import CombineLigandsProteins
import FixedClassificationModel
#import FeatureImportance

CombineLigandsProteins.import_final()
testX = CombineLigandsProteins.X
testY = CombineLigandsProteins.Y

print("imported matrices")

#FeatureImportance.train(testX, testY, CombineLigandsProteins.feats)
FixedClassificationModel.train(testX, testY)
