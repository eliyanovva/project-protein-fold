#This script creates the protein matrix and ligand matrix to train and test the Random Forest Algorithm

#Imports
import SmileKmer
import numpy as np
import ReadingFasta
import labels
import Globals
import Filtering
import Duplicates

#filter_strength = variable to set strength of the kmer filter
#use 'None' to select no filter
#use 'All' to only select kmers that occur exclusively in positive or negative pairs
#otherwise, input an integer value to select the filter strength

#prot_filter_strength = 'None'
prot_filter_strength = 6
#prot_filter_strength = 'All'

#lig_filter_strength = 'None'
lig_filter_strength = 6
#lig_filter_strength = 'All'

#Create classification dictionary
acc_ids = Globals.initialize_protein_list()
logFC, FDR = labels.labels()
classified, pos_counts, neg_counts, pos_pairs, neg_pairs, neutral_pairs = labels.classified_logFC_FDR(logFC, FDR, acc_ids)
#classified = key: protein id, value: (key = ligand, value = {1 if bind, 0 if not bind})
#pos_counts = key: protein id, value: number of positive protein interactions
#neg_counts = key: protein id, value: number of negative protein interactions
#pos_pairs = list of positive protein-ligand pairs; pos_pairs[i] = [protein id, ligand]
#neg_pairs = list of negative protein-ligand pairs; neg_pairs[i] = [protein id, ligand]
#neutral_pairs = list of neutral protein-ligand pairs; neutral_pairs[i] = [protein id, ligand]

proteins_toconsider = set()     #proteins that can form either a positive or negative pair with a ligand

#Extract the proteins and ligands that interact in positive or negative pairs;
#Only these proteins and ligands can be used in the feature matrix
for pair in pos_pairs:
    proteins_toconsider.add(pair[0])
for pair in neg_pairs:
    proteins_toconsider.add(pair[0])

pairs_by_prot = {}      #key = protein id, value = # of pairs involving the protein

#Iterate through all positive and negative prot-lig pairs to update pairs_by_prot
for pair in pos_pairs:
    id = pair[0]
    if id not in pairs_by_prot:
        pairs_by_prot[id] = 0
    pairs_by_prot[id] += 1

for pair in neg_pairs:
    id = pair[0]
    if id not in pairs_by_prot:
        pairs_by_prot[id] = 0
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
#3Di variables
di_features_TM3 = set()
di_seqs_TM3 = {}
di_matrix_TM3 = []
di_features_TM5 = set()
di_seqs_TM5 = {}
di_matrix_TM5 = []
di_features_TM6 = set()
di_seqs_TM6 = {}
di_matrix_TM6 = []
di_features_TM7 = set()
di_seqs_TM7 = {}
di_matrix_TM7 = []

#Create dict of AA sequences only with proteins from pos or neg pairs
AA_dict = Globals.initialize_AA_dict(list(proteins_toconsider))

#Create AA output for TMs 3,5,6,7
AA_seqvar_TM3, AA_features_TM3 = ReadingFasta.make_seqvar_TMS(AA_dict, 0, 5, categorized_seqs_TM3, categorized_features_TM3)
AA_seqvar_TM5, AA_features_TM5 = ReadingFasta.make_seqvar_TMS(AA_dict, 1, 5, categorized_seqs_TM5, categorized_features_TM5)
AA_seqvar_TM6, AA_features_TM6 = ReadingFasta.make_seqvar_TMS(AA_dict, 2, 5, categorized_seqs_TM6, categorized_features_TM6)
AA_seqvar_TM7, AA_features_TM7 = ReadingFasta.make_seqvar_TMS(AA_dict, 3, 5, categorized_seqs_TM7, categorized_features_TM7)

#Filter the AA kmers for TMs 3,5,6,7
AA_filter_TM3, feat1 = Filtering.richness_protein(AA_features_TM3, AA_seqvar_TM3, pos_counts, neg_counts, "TM3", prot_filter_strength)
AA_filter_TM5, feat2 = Filtering.richness_protein(AA_features_TM5, AA_seqvar_TM5, pos_counts, neg_counts, "TM5", prot_filter_strength)
AA_filter_TM6, feat3 = Filtering.richness_protein(AA_features_TM6, AA_seqvar_TM6, pos_counts, neg_counts, "TM6", prot_filter_strength)
AA_filter_TM7, feat4 = Filtering.richness_protein(AA_features_TM7, AA_seqvar_TM7, pos_counts, neg_counts, "TM7", prot_filter_strength)

#Create dict of 3Di sequences only with proteins from pos or neg pairs
Di_dict = Globals.initialize_3Di_dict(list(proteins_toconsider))
#Create 3Di output for TMs 3,5,6,7
Di_seqvar_TM3, Di_features_TM3 = ReadingFasta.make_seqvar_TMS(Di_dict, 0, 5, di_seqs_TM3, di_features_TM3)
Di_seqvar_TM5, Di_features_TM5 = ReadingFasta.make_seqvar_TMS(Di_dict, 1, 5, di_seqs_TM5, di_features_TM5)
Di_seqvar_TM6, Di_features_TM6 = ReadingFasta.make_seqvar_TMS(Di_dict, 2, 5, di_seqs_TM6, di_features_TM6)
Di_seqvar_TM7, Di_features_TM7 = ReadingFasta.make_seqvar_TMS(Di_dict, 3, 5, di_seqs_TM7, di_features_TM7)

#Filter the 3Di kmers for TMs 3,5,6,7
Di_filter_TM3, feat5 = Filtering.richness_protein(Di_features_TM3, Di_seqvar_TM3, pos_counts, neg_counts, "TM3", prot_filter_strength)
Di_filter_TM5, feat6 = Filtering.richness_protein(Di_features_TM5, Di_seqvar_TM5, pos_counts, neg_counts, "TM5", prot_filter_strength)
Di_filter_TM6, feat7 = Filtering.richness_protein(Di_features_TM6, Di_seqvar_TM6, pos_counts, neg_counts, "TM6", prot_filter_strength)
Di_filter_TM7, feat8 = Filtering.richness_protein(Di_features_TM7, Di_seqvar_TM7, pos_counts, neg_counts, "TM7", prot_filter_strength)

AA_seqvar = [AA_seqvar_TM3, AA_seqvar_TM5, AA_seqvar_TM6, AA_seqvar_TM7]
AA_feat = [AA_filter_TM3, AA_filter_TM5, AA_filter_TM6, AA_filter_TM7]
Di_seqvar = [Di_seqvar_TM3, Di_seqvar_TM5, Di_seqvar_TM6, Di_seqvar_TM7]
Di_feat = [Di_filter_TM3, Di_filter_TM5, Di_filter_TM6, Di_filter_TM7]

#Extract proteins with unique AA and 3di kmer frequencies
unique_proteins = Duplicates.remove_proteins(AA_seqvar, AA_feat, Di_seqvar, Di_feat, pairs_by_prot, list(proteins_toconsider))

pos_dict = {}       #key = protein id, value = list of ligands that the protein binds with
neg_dict = {}       #key = protein id, value = list of ligands that the protein does not bind with

for id in unique_proteins:
    pos_dict[id] = []
    neg_dict[id] = []

ligands_from_unip = set()       #set of ligands that form pos or neg pairs with the set of unique proteins

#Iterate through all positive prot-lig pairs to update pos_dict
for pair in pos_pairs:
    id = pair[0]
    lig = pair[1]
    #Check that id is from the set of unique proteins
    if (id in pos_dict):
        pos_dict[id].append(lig)
        ligands_from_unip.add(lig)
#Iterate through all negative prot-lig pairs to update neg_dict
for pair in neg_pairs:
    id = pair[0]
    lig = pair[1]
    #Check that id is from the set of unique proteins
    if (id in neg_dict):
        neg_dict[id].append(lig)
        ligands_from_unip.add(lig)
#Import dictionary matching ligands to SMILES String
ligand_dict = Globals.initialize_ligand_dict()

#Create ligands matrix
ligand_features, ligand_counts = SmileKmer.ligand_matrix(ligand_dict, 5, ligands_from_unip)
#ligand_features = list of kmers found in the ligands with pos / neg pairs
#ligand_counts = key: ligand, value: dict (key: kmer, value: freq. of kmer in the ligand)
#   ligand_counts only uses ligands from ligands_from_unip as keys

pos_by_lig = {}         #key: ligand, value: # of pos. pairs with the ligand
neg_by_lig = {}         #key: ligand, value: # of neg. pairs with the ligand
total_by_lig = {}       #key: ligand, value: # of all pairs the involve the ligand

for lig in ligands_from_unip:
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

#Update ligand_counts to only use filtered kmers
lig_counts_filter = Filtering.richness_ligand(ligand_counts, pos_by_lig, neg_by_lig, lig_filter_strength)
#Extract ligands with unique kmer frequencies
unique_ligands = Duplicates.remove_ligands(lig_counts_filter, total_by_lig)

pos_total = 0           #num. of positive pairs with ligands from unique_ligands
neg_total = 0           #num. of negative pairs with ligands from unique_ligands
lig_mat = []            #matrix to store ligand features
for id in unique_proteins:
    for lig in pos_dict[id]:
        if (unique_ligands.count(lig) != 0):
            lig_mat.append(np.array(list(lig_counts_filter[lig].values())))
            pos_total += 1

for id in unique_proteins:
    for lig in neg_dict[id]:
        if (unique_ligands.count(lig) != 0):
            lig_mat.append(np.array(list(lig_counts_filter[lig].values())))
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

pos_Di_mat_TM3 = ReadingFasta.makematrix(Di_seqvar_TM3, Di_filter_TM3, di_matrix_TM3, unique_ligands, pos_dict)
pos_Di_mat_TM5 = ReadingFasta.makematrix(Di_seqvar_TM5, Di_filter_TM5, di_matrix_TM5, unique_ligands, pos_dict)
pos_Di_mat_TM6 = ReadingFasta.makematrix(Di_seqvar_TM6, Di_filter_TM6, di_matrix_TM6, unique_ligands, pos_dict)
pos_Di_mat_TM7 = ReadingFasta.makematrix(Di_seqvar_TM7, Di_filter_TM7, di_matrix_TM7, unique_ligands, pos_dict)

Di_mat_TM3 = ReadingFasta.makematrix(Di_seqvar_TM3, Di_filter_TM3, pos_Di_mat_TM3, unique_ligands, neg_dict)
Di_mat_TM5 = ReadingFasta.makematrix(Di_seqvar_TM5, Di_filter_TM5, pos_Di_mat_TM5, unique_ligands, neg_dict)
Di_mat_TM6 = ReadingFasta.makematrix(Di_seqvar_TM6, Di_filter_TM6, pos_Di_mat_TM6, unique_ligands, neg_dict)
Di_mat_TM7 = ReadingFasta.makematrix(Di_seqvar_TM7, Di_filter_TM7, pos_Di_mat_TM7, unique_ligands, neg_dict)

Di_matrix = np.concatenate((np.array(Di_mat_TM3, dtype = np.uint8), np.array(Di_mat_TM5, dtype = np.uint8),
                            np.array(Di_mat_TM6, dtype = np.uint8), np.array(Di_mat_TM7, dtype = np.uint8)) , axis = 1)

#Concatenate AA and 3Di matrices
#intermed_matrix = np.concatenate((np.array(AA_mat, dtype = np.uint8), np.array(Di_mat, dtype = np.uint8)) , axis = 1)
intermed_matrix = np.concatenate((np.array(AA_matrix, dtype = np.uint8), np.array(Di_matrix, dtype = np.uint8)) , axis = 1)

#Concatenate protein and ligand matrices
final_matrix = np.concatenate((intermed_matrix, np.array(lig_mat, dtype = np.uint8)), axis = 1)

#Create Classification Vector

pos_array = np.repeat(1, int(pos_total))
neg_array = np.repeat(0, int(neg_total))
logFCmat = np.concatenate((pos_array, neg_array), axis=0)

print('Finished Part 1')

acc_set = set(acc_ids)
neutral_proteins = acc_set.difference(proteins_toconsider) #https://www.w3schools.com/python/ref_set_difference.asp

n_protein_list = list(neutral_proteins)
n_protein_list.sort(reverse=True)

neutral_ligands = Globals.initialize_ligand_list()
neutral_ligands.sort()

nAA_dict = Globals.initialize_AA_dict(n_protein_list)
nDi_dict = Globals.initialize_3Di_dict(n_protein_list)

nAA_seqvar_TM3, ignoreAA3 = ReadingFasta.make_seqvar_TMS(nAA_dict, 0, 5, {}, set())
nAA_seqvar_TM5, ignoreAA5 = ReadingFasta.make_seqvar_TMS(nAA_dict, 1, 5, {}, set())
nAA_seqvar_TM6, ignoreAA6 = ReadingFasta.make_seqvar_TMS(nAA_dict, 2, 5, {}, set())
nAA_seqvar_TM7, ignoreAA7 = ReadingFasta.make_seqvar_TMS(nAA_dict, 3, 5, {}, set())

nDi_seqvar_TM3, ignoreDi3 = ReadingFasta.make_seqvar_TMS(nDi_dict, 0, 5, {}, set())
nDi_seqvar_TM5, ignoreDi5 = ReadingFasta.make_seqvar_TMS(nDi_dict, 1, 5, {}, set())
nDi_seqvar_TM6, ignoreDi6 = ReadingFasta.make_seqvar_TMS(nDi_dict, 2, 5, {}, set())
nDi_seqvar_TM7, ignoreDi7 = ReadingFasta.make_seqvar_TMS(nDi_dict, 3, 5, {}, set())

nAA_seqvar = [nAA_seqvar_TM3, nAA_seqvar_TM5, nAA_seqvar_TM6, nAA_seqvar_TM7]
nDi_seqvar = [nDi_seqvar_TM3, nDi_seqvar_TM5, nDi_seqvar_TM6, nDi_seqvar_TM7]
n_unip = Duplicates.n_remove_proteins(nAA_seqvar, AA_feat, nDi_seqvar, Di_feat, n_protein_list)

filter_kmers = list(lig_counts_filter['pS6_DE_1p_dimethyltrisulfide.csv'].keys())
n_lig_counts = SmileKmer.n_ligand_matrix(ligand_dict, 5, neutral_ligands, filter_kmers)
n_uni_lig = Duplicates.n_remove_ligands(n_lig_counts)

num_ligands = len(n_uni_lig)

nlig_mat = []
for lig in n_uni_lig:
    nlig_mat.append(np.array(list(n_lig_counts[lig].values())))

nAA_mat_TM3 = ReadingFasta.make_nmatrix(nAA_seqvar_TM3, AA_filter_TM3, [], n_unip, num_ligands)
nAA_mat_TM5 = ReadingFasta.make_nmatrix(nAA_seqvar_TM5, AA_filter_TM5, [], n_unip, num_ligands)
nAA_mat_TM6 = ReadingFasta.make_nmatrix(nAA_seqvar_TM6, AA_filter_TM6, [], n_unip, num_ligands)
nAA_mat_TM7 = ReadingFasta.make_nmatrix(nAA_seqvar_TM7, AA_filter_TM7, [], n_unip, num_ligands)

nDi_mat_TM3 = ReadingFasta.make_nmatrix(nDi_seqvar_TM3, Di_filter_TM3, [], n_unip, num_ligands)
nDi_mat_TM5 = ReadingFasta.make_nmatrix(nDi_seqvar_TM5, Di_filter_TM5, [], n_unip, num_ligands)
nDi_mat_TM6 = ReadingFasta.make_nmatrix(nDi_seqvar_TM6, Di_filter_TM6, [], n_unip, num_ligands)
nDi_mat_TM7 = ReadingFasta.make_nmatrix(nDi_seqvar_TM7, Di_filter_TM7, [], n_unip, num_ligands)

nAA_mat = np.concatenate((np.array(nAA_mat_TM3, dtype= np.uint8), np.array(nAA_mat_TM5, dtype= np.uint8),
                          np.array(nAA_mat_TM6, dtype= np.uint8), np.array(nAA_mat_TM7, dtype= np.uint8)), axis = 1)

nDi_mat = np.concatenate((np.array(nDi_mat_TM3, dtype= np.uint8), np.array(nDi_mat_TM5, dtype= np.uint8),
                          np.array(nDi_mat_TM6, dtype= np.uint8), np.array(nDi_mat_TM7, dtype= np.uint8)), axis = 1)

n_intermed = np.concatenate((np.array(nAA_mat, dtype = np.uint8), np.array(nDi_mat, dtype = np.uint8)) , axis = 1)
n_final_lig = np.repeat(nlig_mat, len(n_unip), axis = 0)

n_final_mat = np.concatenate((n_intermed, n_final_lig), axis=1)
"""
print(len(n_unip))                  #404
print(len(neutral_ligands))         #55
print(len(filter_kmers))            #182
print(len(n_intermed))              #21412
print(len(n_intermed[0]))           #10547
print(len(n_uni_lig))               #53
print(len(n_final_lig))             #21412
print(len(n_final_lig[0]))          #182
"""

for lig in n_uni_lig:
    print(lig)

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
    protmat = intermed_matrix
    #For TrainandTest.py
    global X
    X = final_matrix
    global Y
    Y = logFCmat
    global nMat
    nMat = n_final_mat
    global nproteins
    nproteins = n_unip
    global nligands
    nligands = n_uni_lig
    global feats
    feat1.extend(feat2)
    feat1.extend(feat3)
    feat1.extend(feat4)
    feat1.extend(feat5)
    feat1.extend(feat6)
    feat1.extend(feat7)
    feat1.extend(feat8)
    feat1.extend(ligand_features)
    feats = feat1