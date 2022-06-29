#This script tests the Random Forest on the entire preprocessed data set

import CombineLigandsProteins
import FeatureImportance

CombineLigandsProteins.import_final()
testX = CombineLigandsProteins.X
testY = CombineLigandsProteins.Y
FeatureImportance.train(testX, testY)