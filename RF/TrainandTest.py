#This script trains and test the Random Forest 

#Imports
import CombineLigandsProteins
#import Sequence_only
import FixedClassificationModel
#import PredictNewCombos
#import Structure_only
#import AdjustingThreshold
#import Feature_Importance.FeatureImportance as fi
#import Metrics.Metrics_Graphs as Metrics_Graphs

#Initialize matrices
CLP_vars = CombineLigandsProteins.develop_matrices('../Ligands_withSMILE/ligand_SMILEs.csv', "../data_files/TMdomains/TM.csv",
                        "../data_files/3DiSequences/fullset_ss.fasta")
testX = CLP_vars['X']
testY = CLP_vars['Y']
logFC = CLP_vars['logFC_data']
FDR = CLP_vars['FDR_data']
BALANCE = CLP_vars['balance']

#Predicting with new combinations
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
"""

#Test with protein sequence only
#Sequence_only.import_final()
#seqX = Sequence_only.X
#seqY = Sequence_only.Y

#Test with protein structure only
#Structure_only.import_final()
#structX = Structure_only.X
#structY = Structure_only.Y

#Find the best threshold for classification
#AdjustingThreshold.train(testX, testY)

#Find the most important features
#fi.importance_file(testX, testY, CombineLigandsProteins.feats)

#Perform classic training and testing of the model
#FixedClassificationModel.train(testX, testY, BALANCE)

#Create graphs of Precision-Recall and Receiver Operating Characteristic curves
#Metrics_Graphs.train(testX, testY)

#Examine TP and FN rates

accuracy = 0
recall = 0
BAC = 0
MAT = 0

a_TN = 0
a_FN = 0
a_TP = 0
a_FP = 0

loss = 0
for i in range(50):
    print("run " + str(i))
    acc, rec, bac, mat, TN, FN, TP, FP, log_loss = FixedClassificationModel.train(testX, testY, BALANCE)
    accuracy += acc
    recall += rec
    BAC += bac
    MAT += mat
    loss += log_loss

    a_TN += TN
    a_FN += FN
    a_TP += TP
    a_FP += FP



print('Average Accuracy: ' + str(accuracy/50))
print('Average Recall: ' + str(recall/50))
print('Average Balanced: ' + str(BAC/50))

print('Average TN: ' + str(a_TN/50))
print('Average FN: ' + str(a_FN/50))
print('Average TP: ' + str(a_TP/50))
print('Average FP: ' + str(a_FP/50))

print('Average Log Loss: ' + str(loss/50))

