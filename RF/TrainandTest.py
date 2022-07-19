#This script tests the Random Forest on the entire preprocessed data set
import CombineLigandsProteins
#import Sequence_only
import FixedClassificationModel
#import Structure_only
#import AdjustingThreshold
#import Feature_Importance.FeatureImportance as fi
#import Metrics_Graphs

CombineLigandsProteins.import_final()
testX = CombineLigandsProteins.X
testY = CombineLigandsProteins.Y

n_testX = CombineLigandsProteins.nMat
n_proteins = CombineLigandsProteins.nproteins
n_ligands = CombineLigandsProteins.nligands

neutral_dict = {}
for protein in n_proteins:
    row = {}
    for lig in n_ligands:
        row[lig] = 0
    neutral_dict[protein] = row

num_ligands = len(n_ligands)

for i in range(50):
    print('run ' + str(i))
    n_pred = FixedClassificationModel.neutral_train(testX, testY, n_testX)

    j = 0
    for p in n_pred:
        p_ind = j // num_ligands
        if p_ind == 0:
            l_ind = j
        else:
            l_ind = j % num_ligands
        neutral_dict[n_proteins[p_ind]][n_ligands[l_ind]] += p
        j += 1

f = open('neutral_predict_none.csv', 'w')

f.write('Protein,')
for lig in n_ligands:
    f.write(lig + ",")
f.write("\n")

for id in n_proteins:
    f.write(id + ",")
    for lig in n_ligands:
        f.write(str(neutral_dict[id][lig]) + ',')
    f.write("\n")

#Sequence_only.import_final()
#seqX = Sequence_only.X
#seqY = Sequence_only.Y

#Structure_only.import_final()
#structX = Structure_only.X
#structY = Structure_only.Y

#AdjustingThreshold.train(testX, testY)

#fi.importance_file(testX, testY, CombineLigandsProteins.feats)
#FixedClassificationModel.train(testX, testY)
#Metrics_Graphs.train(testX, testY)

"""
accuracy = 0
recall = 0
BAC = 0
MAT = 0

f = open('results_true_false_vals_struct.csv', 'w')
f.write('Run,TN,FN,TP,FP' + "\n")

for i in range(50):
    print("run " + str(i))
    acc, rec, bac, mat, TN, FN, TP, FP = FixedClassificationModel.train(testX, testY)
    accuracy += acc
    recall += rec
    BAC += bac
    MAT += mat
    f.write(str(i+1) + ", " + str(TN) + ", " + str(FN) + ", " + str(TP) + ", " + str(FP) + "\n")

print('Average Accuracy: ' + str(accuracy/50))
print('Average Recall: ' + str(recall/50))
print('Average Balanced: ' + str(BAC/50))
#print('Average Matthew: ' + str(mat/50))
#f.close()
"""