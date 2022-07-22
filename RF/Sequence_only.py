#This script creates matrices containing only protein sequence and ligand data

#Imports
import SmileKmer
import numpy as np
import ReadingFasta
import labels
import Globals
import Filtering
import Duplicates

#Create classification dictionary
acc_ids = Globals.initialize_protein_list()
logFC, FDR = labels.labels()
classified, pos_counts, neg_counts, pos_pairs, neg_pairs = labels.classified_logFC_FDR(logFC, FDR, acc_ids)

proteins_toconsider = set()     #has 392 proteins
ligands_toconsider = set()      #has 49 ligands

for pair in pos_pairs:
    proteins_toconsider.add(pair[0])
    ligands_toconsider.add(pair[1])

for pair in neg_pairs:
    proteins_toconsider.add(pair[0])
    ligands_toconsider.add(pair[1])

#pos_pairs = 565
#neg_pairs = 236
#proteins_toconsider = 392

pairs_by_prot = {}

for id in proteins_toconsider:
    pairs_by_prot[id] = 0

for pair in pos_pairs:
    id = pair[0]
    pairs_by_prot[id] += 1

for pair in neg_pairs:
    id = pair[0]
    pairs_by_prot[id] += 1

#Initialize Variables
#categorized variables
categorized_features_TM3 = set()
categorized_seqs_TM3 = {}
categorized_matrix_TM3 = []
categorized_features_TM5 = set()
categorized_seqs_TM5 = {}
categorized_matrix_TM5 = []
categorized_features_TM6 = set()
categorized_seqs_TM6 = {}
categorized_matrix_TM6 = []
categorized_features_TM7 = set()
categorized_seqs_TM7 = {}
categorized_matrix_TM7 = []



#Create AA output for TMs 3,5,6,7
AA_dict = Globals.initialize_AA_dict(list(proteins_toconsider))   #create dict with proteins from pos / neg pairs
AA_seqvar_TM3, AA_features_TM3 = ReadingFasta.make_seqvar_TMS(AA_dict, 0, 5, categorized_seqs_TM3, categorized_features_TM3)
AA_seqvar_TM5, AA_features_TM5 = ReadingFasta.make_seqvar_TMS(AA_dict, 1, 5, categorized_seqs_TM5, categorized_features_TM5)
AA_seqvar_TM6, AA_features_TM6 = ReadingFasta.make_seqvar_TMS(AA_dict, 2, 5, categorized_seqs_TM6, categorized_features_TM6)
AA_seqvar_TM7, AA_features_TM7 = ReadingFasta.make_seqvar_TMS(AA_dict, 3, 5, categorized_seqs_TM7, categorized_features_TM7)

AA_filter_TM3, feat1 = Filtering.richness_protein(AA_features_TM3, AA_seqvar_TM3, pos_counts, neg_counts, "TM3")
AA_filter_TM5, feat2 = Filtering.richness_protein(AA_features_TM5, AA_seqvar_TM5, pos_counts, neg_counts, "TM5")
AA_filter_TM6, feat3 = Filtering.richness_protein(AA_features_TM6, AA_seqvar_TM6, pos_counts, neg_counts, "TM6")
AA_filter_TM7, feat4 = Filtering.richness_protein(AA_features_TM7, AA_seqvar_TM7, pos_counts, neg_counts, "TM7")

AA_seqvar = [AA_seqvar_TM3, AA_seqvar_TM5, AA_seqvar_TM6, AA_seqvar_TM7]
AA_feat = [AA_filter_TM3, AA_filter_TM5, AA_filter_TM6, AA_filter_TM7]

unique_proteins = Duplicates.remove_proteins(AA_seqvar, AA_feat, [], [], pairs_by_prot)

pos_dict = {}
neg_dict = {}

for id in unique_proteins:
    pos_dict[id] = []
    neg_dict[id] = []

for pair in pos_pairs:
    id = pair[0]
    lig = pair[1]
    if (id in pos_dict) & (lig in ligands_toconsider):
        pos_dict[id].append(lig)

for pair in neg_pairs:
    id = pair[0]
    lig = pair[1]
    if (id in neg_dict) & (lig in ligands_toconsider):
        neg_dict[id].append(lig)

#Import dictionary matching ligands to SMILES String
ligand_dict = Globals.initialize_ligand_dict()
#Create ligands matrix
ligand_features, ligand_counts = SmileKmer.ligand_matrix(ligand_dict, 5, ligands_toconsider)

#ligand_counts only has ligands from ligands_toconsider as keys

pos_by_lig = {}
neg_by_lig = {}
total_by_lig = {}

for lig in ligands_toconsider:
    pos_by_lig[lig] = 0
    neg_by_lig[lig] = 0
    total_by_lig[lig] = 0

for id in unique_proteins:
    p_pairs = pos_dict[id]
    for lig in p_pairs:
        pos_by_lig[lig] += 1
        total_by_lig[lig] += 1

    n_pairs = neg_dict[id]
    for lig in n_pairs:
        neg_by_lig[lig] += 1
        total_by_lig[lig] += 1

lig_counts_filter = Filtering.richness_ligand(ligand_counts, pos_by_lig, neg_by_lig)
unique_ligands = Duplicates.remove_ligands(ligand_counts, total_by_lig)

pos_total = 0
neg_total = 0
lig_mat = []
for id in unique_proteins:
    for lig in pos_dict[id]:
        if (unique_ligands.count(lig) != 0):
            lig_mat.append(np.array(list(ligand_counts[lig].values())))
            pos_total += 1

for id in unique_proteins:
    for lig in neg_dict[id]:
        if (unique_ligands.count(lig) != 0):
            lig_mat.append(np.array(list(ligand_counts[lig].values())))
            neg_total += 1

pos_AA_mat_TM3 = ReadingFasta.makematrix(AA_seqvar_TM3, AA_filter_TM3, categorized_matrix_TM3, unique_ligands, pos_dict)
pos_AA_mat_TM5 = ReadingFasta.makematrix(AA_seqvar_TM5, AA_filter_TM5, categorized_matrix_TM5, unique_ligands, pos_dict)
pos_AA_mat_TM6 = ReadingFasta.makematrix(AA_seqvar_TM6, AA_filter_TM6, categorized_matrix_TM6, unique_ligands, pos_dict)
pos_AA_mat_TM7 = ReadingFasta.makematrix(AA_seqvar_TM7, AA_filter_TM7, categorized_matrix_TM7, unique_ligands, pos_dict)

AA_mat_TM3 = ReadingFasta.makematrix(AA_seqvar_TM3, AA_filter_TM3, pos_AA_mat_TM3, unique_ligands, neg_dict)
AA_mat_TM5 = ReadingFasta.makematrix(AA_seqvar_TM5, AA_filter_TM5, pos_AA_mat_TM5, unique_ligands, neg_dict)
AA_mat_TM6 = ReadingFasta.makematrix(AA_seqvar_TM6, AA_filter_TM6, pos_AA_mat_TM6, unique_ligands, neg_dict)
AA_mat_TM7 = ReadingFasta.makematrix(AA_seqvar_TM7, AA_filter_TM7, pos_AA_mat_TM7, unique_ligands, neg_dict)

AA_matrix = np.concatenate((np.array(AA_mat_TM3, dtype = np.uint8), np.array(AA_mat_TM5, dtype = np.uint8),
                            np.array(AA_mat_TM6, dtype = np.uint8), np.array(AA_mat_TM7, dtype = np.uint8)) , axis = 1)

#Concatenate protein and ligand matrices
final_matrix = np.concatenate((AA_matrix, np.array(lig_mat, dtype = np.uint8)), axis = 1)

#Create Classification Vector

pos_array = np.repeat(1, int(pos_total))
neg_array = np.repeat(0, int(neg_total))
logFCmat = np.concatenate((pos_array, neg_array), axis=0)

#Return the number of repeated entries. Adapted from: https://www.geeksforgeeks.org/print-unique-rows/
def uniquematrix(matrix):
    rowCount = len(matrix)
    if rowCount == 0:
        return
    columnCount = len(matrix[0])
    if columnCount == 0:
        return
    unique = {}
    ret = 0
    for row in matrix:
        rowkey =  " ".join(["%s"] * columnCount) % tuple(row)
        if rowkey not in unique:
            unique[rowkey] = True
        else:
            ret += 1
    return ret

#unique returned

#The following code checks whether all entries are unique
#print(uniquematrix(intermed_matrix))

def import_final():
    #For No3Di.py
    global AA
    #AA = AA_mat
    #For OneLigandRF.py and OneProteinRF.py
    global proteins
    #proteins = seqvar1
    global dictionary
    dictionary = classified
    global protmat
    protmat = AA_matrix
    #For TrainandTest.py
    global X
    X = final_matrix
    global Y
    Y = logFCmat
    global feats
    feat1.extend(feat2)
    feat1.extend(feat3)
    feat1.extend(feat4)
    feat1.extend(ligand_features)
    feats = feat1
