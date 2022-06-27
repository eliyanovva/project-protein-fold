#This script tests the Random Forest on the entire preprocessed data set

import CombineLigandsProteins
import AdjustingThreshold

CombineLigandsProteins.import_final()
testX = CombineLigandsProteins.X
testY = CombineLigandsProteins.Y
AdjustingThreshold.train(testX, testY)