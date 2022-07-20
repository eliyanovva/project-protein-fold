#This script tests the Random Forest on the entire preprocessed data set
import CombineLigandsProteins
#import Sequence_only
import FixedClassificationModel
import PredictNewCombos
#import Structure_only
#import AdjustingThreshold
#import Feature_Importance.FeatureImportance as fi
#import Metrics_Graphs

CombineLigandsProteins.import_final()
testX = CombineLigandsProteins.X
testY = CombineLigandsProteins.Y
logFC = CombineLigandsProteins.logFC_data
FDR = CombineLigandsProteins.FDR_data

"""
PredictNewCombos.import_final()
newX = PredictNewCombos.X
new_combos = PredictNewCombos.combos
all_ligs = set()

combo_list = []
for id in new_combos:
    for lig in new_combos[id]:
        combo_list.append([id, lig])
        all_ligs.add(lig)

combo_dict = {}
for id in new_combos:
    row = {}
    for lig in new_combos[id]:
        row[lig] = 0
    combo_dict[id] = row

for i in range(50):
    print('run ' + str(i))
    n_pred = FixedClassificationModel.neutral_train(testX, testY, newX)

    for j in range(len(n_pred)):
        combo = combo_list[j]
        combo_dict[combo[0]][combo[1]] += n_pred[j]

f = open('new_combo_prediction.csv', 'w')

f.write('Protein')
for lig in all_ligs:
    f.write("," + lig)
f.write("\n")

for id in new_combos:
    f.write(id)
    for lig in all_ligs:
        if new_combos[id].count(lig) == 0:
            f.write(",-")
        else:
            f.write("," + str(combo_dict[id][lig]))
    f.write("\n")
f.close()

f = open('high_FDR6_pairs.txt', 'w')
f.write("Protein,Ligand,logFC,FDR,Positive Obs,Classification" + "\n")

for id in new_combos:
    for csv in new_combos[id]:
        FC_val = logFC[id][csv]
        pos_obs = combo_dict[id][csv]
        if (FDR[id][csv] > .15) & (FDR[id][csv] <= .6):
            f.write(id + "," + csv + "," + str(logFC[id][csv]) +
                    "," + str(FDR[id][csv]) + "," + str(combo_dict[id][csv]))
            if (FC_val < 1) & (pos_obs <= 25):
                f.write(", TN" + "\n")
            elif (FC_val >= 1) & (pos_obs <= 25):
                f.write(", FN" + "\n")
            elif (FC_val < 1) & (pos_obs > 25):
                f.write(", FP" + "\n")
            elif (FC_val >= 1) & (pos_obs > 25):
                f.write(", TP" + "\n")
f.close()

"""

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