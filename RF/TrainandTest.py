#This script tests the Random Forest on the entire preprocessed data set
import CombineLigandsProteins
#import Sequence_only
#import FixedClassificationModel
#import PredictNewCombos
#import Structure_only
#import AdjustingThreshold
#import Feature_Importance.FeatureImportance as fi
import RF.Metrics.Metrics_Graphs as Metrics_Graphs

CombineLigandsProteins.import_final()
testX = CombineLigandsProteins.X
testY = CombineLigandsProteins.Y
logFC = CombineLigandsProteins.logFC_data
FDR = CombineLigandsProteins.FDR_data
BALANCE = CombineLigandsProteins.balance
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
    n_pred = FixedClassificationModel.train_new_pairs(testX, testY, newX, BALANCE)

    for j in range(len(n_pred)):
        combo = combo_list[j]
        combo_dict[combo[0]][combo[1]] += n_pred[j]

#f1 = open('results_true_false_newpairs_FDR1.csv', 'w')
#f1.write("Protein,Ligand,logFC,FDR,Positive Obs,Classification" + "\n")

i = 0
j = 0
new_pos = 0
new_neg = 0
for id in new_combos:
    for csv in new_combos[id]:
        i += 1
        FC_val = logFC[id][csv]
        pos_obs = combo_dict[id][csv]
        if (FDR[id][csv] > .15) & (FDR[id][csv] <= .4):
            j += 1
            if FC_val >= 1:
                new_pos += 1
            else:
                new_neg += 1

print('New Pos Pairs: ' + str(new_pos))
print('New Neg Pairs: ' + str(new_neg))
print(i)
print(j)

        if (FDR[id][csv] > .1) & (FDR[id][csv] <= .4):
            f1.write(id + "," + csv + "," + str(logFC[id][csv]) +
                    "," + str(FDR[id][csv]) + "," + str(combo_dict[id][csv]))
            if (FC_val < 1) & (pos_obs <= 25):
                f1.write(", TN" + "\n")
            elif (FC_val >= 1) & (pos_obs <= 25):
                f1.write(", FN" + "\n")
            elif (FC_val < 1) & (pos_obs > 25):
                f1.write(", FP" + "\n")
            elif (FC_val >= 1) & (pos_obs > 25):
                f1.write(", TP" + "\n")
            
#f1.close()

#Sequence_only.import_final()
#seqX = Sequence_only.X
#seqY = Sequence_only.Y

#Structure_only.import_final()
#structX = Structure_only.X
#structY = Structure_only.Y

#AdjustingThreshold.train(testX, testY)

#fi.importance_file(testX, testY, CombineLigandsProteins.feats)
#FixedClassificationModel.train(testX, testY, BALANCE)
"""
Metrics_Graphs.train(testX, testY)

"""
accuracy = 0
recall = 0
BAC = 0
MAT = 0

f1 = open('FiltAll_FDR1_ROC.csv', 'w')
f2 = open('FiltAll_FDR1_PreRec.csv', 'w')

f1.write('Run,TN,FP' + "\n")
f2.write('Run,Precision,Recall' + "\n")

for i in range(50):
    print("run " + str(i))
    acc, rec, bac, mat, TN, FN, TP, FP = FixedClassificationModel.train(testX, testY, BALANCE)
    accuracy += acc
    recall += rec
    BAC += bac
    MAT += mat
    #f2.write(str(i+1) + ", " + str(TN) + ", " + str(FN) + ", " + str(TP) + ", " + str(FP) + "\n")
    #f2.write(str(i+1)+","+str(acc) + "\n")
print('Average Accuracy: ' + str(accuracy/50))
print('Average Recall: ' + str(recall/50))
#print('Average Balanced: ' + str(BAC/50))
#print('Average Matthew: ' + str(mat/50))
#f2.close()
"""