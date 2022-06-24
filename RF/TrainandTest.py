import CombineLigandsProteins
import FixedClassificationModel

CombineLigandsProteins.import_final()
testX = CombineLigandsProteins.X
testY = CombineLigandsProteins.Y
FixedClassificationModel.train(testX, testY)