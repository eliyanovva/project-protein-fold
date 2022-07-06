#This script tests the Random Forest on the entire preprocessed data set

import CombineLigandsProteins
import FixedClassificationModel
#import FeatureImportance

CombineLigandsProteins.import_final()
testX = CombineLigandsProteins.X
testY = CombineLigandsProteins.Y
protein_freqs = CombineLigandsProteins.all_protein
ligand_freqs = CombineLigandsProteins.all_ligand
print("imported matrices")

#FeatureImportance.train(testX, testY, CombineLigandsProteins.feats)
FixedClassificationModel.train(testX, testY, protein_freqs, ligand_freqs)
