#This script tests the Random Forest on the entire preprocessed data set

import CombineLigandsProteins
import FixedClassificationModel
#import FeatureImportance

CombineLigandsProteins.import_final()
testX = CombineLigandsProteins.X
testY = CombineLigandsProteins.Y
protein_freqs = CombineLigandsProteins.all_protein
ligand_freqs = CombineLigandsProteins.all_ligand

"""
#FeatureImportance.train(testX, testY, CombineLigandsProteins.feats)
FixedClassificationModel.train(testX, testY, protein_freqs, ligand_freqs)
"""

accuracy = 0
recall = 0
for i in range(5):
    print("run " + str(i))
    acc,rec = FixedClassificationModel.train(testX, testY, protein_freqs, ligand_freqs)
    accuracy += acc
    recall += rec
print('Average Accuracy: ' + str(accuracy/5))
print('Average Recall: ' + str(recall/5))
