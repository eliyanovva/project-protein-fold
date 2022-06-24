#This script tests the Random Forest on the entire preprocessed data set

import CombineLigandsProteins
import FixedClassificationModel

CombineLigandsProteins.import_final()
testX = CombineLigandsProteins.X
testY = CombineLigandsProteins.Y
FixedClassificationModel.train(testX, testY)