#This script tests the Random Forest on the entire preprocessed data set

import CombineLigandsProteins
import FeatureImportance

CombineLigandsProteins.import_final()
testX = CombineLigandsProteins.X
testY = CombineLigandsProteins.Y
AA_feat = CombineLigandsProteins.feat1
Di_feat = CombineLigandsProteins.feat2
Lig_feat = CombineLigandsProteins.feat3
FeatureImportance.train(testX, testY, AA_feat, Di_feat, Lig_feat)