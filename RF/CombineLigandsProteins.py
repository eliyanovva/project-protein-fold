#This script creates the protein matrix and ligand matrix to train and test the Random Forest Algorithm

#Imports
import SmileKmer
import numpy as np
import ReadingFasta
import labels
import Globals
import Filtering

#Create classification dictionary
logFC, FDR = labels.labels()
classified, pos_counts, neg_counts = labels.classified_logFC_FDR(logFC, FDR)

#Initialize Variables
#categorized variables
categorized_features_TM3 = set()
categorized_seqs_TM3 = []
categorized_matrix_TM3 = []
categorized_features_TM5 = set()
categorized_seqs_TM5 = []
categorized_matrix_TM5 = []
categorized_features_TM6 = set()
categorized_seqs_TM6 = []
categorized_matrix_TM6 = []
categorized_features_TM7 = set()
categorized_seqs_TM7 = []
categorized_matrix_TM7 = []
#3Di variables
di_features_TM3 = set()
di_seqs_TM3 = []
di_matrix_TM3 = []
di_features_TM5 = set()
di_seqs_TM5 = []
di_matrix_TM5 = []
di_features_TM6 = set()
di_seqs_TM6 = []
di_matrix_TM6 = []
di_features_TM7 = set()
di_seqs_TM7 = []
di_matrix_TM7 = []

"""
#Creating output for categorized amino acids
#Read fasta file
fasta1 = open("../data_files/AminoAcidSequences/fully_categorized.fasta")
#Create kmer frequency dictionary
seqvar1, features1 = ReadingFasta.make_seqvar(fasta1, categorized_seqs, categorized_features)
#Remove insignificant kmers
filter_feat = Filtering.richness_protein(features1, seqvar1, pos_counts, neg_counts)
# Make the matrix
AA_mat = ReadingFasta.makematrix(seqvar1, filter_feat, categorized_matrix)
"""

#Create AA output for TMs 3,5,6,7
AA_dict = Globals.initialize_AA_dict()
AA_seqvar_TM3, AA_features_TM3 = ReadingFasta.make_seqvar_TMS(AA_dict, 0, 5, categorized_seqs_TM3, categorized_features_TM3)
AA_seqvar_TM5, AA_features_TM5 = ReadingFasta.make_seqvar_TMS(AA_dict, 1, 5, categorized_seqs_TM5, categorized_features_TM5)
AA_seqvar_TM6, AA_features_TM6 = ReadingFasta.make_seqvar_TMS(AA_dict, 2, 5, categorized_seqs_TM6, categorized_features_TM6)
AA_seqvar_TM7, AA_features_TM7 = ReadingFasta.make_seqvar_TMS(AA_dict, 3, 5, categorized_seqs_TM7, categorized_features_TM7)


print(AA_features_TM3)
print(len(AA_features_TM3))
print(AA_features_TM5)
print(len(AA_features_TM5))
print(AA_features_TM6)
print(len(AA_features_TM6))
print(AA_features_TM7)
print(len(AA_features_TM7))

AA_mat_TM3 = ReadingFasta.makematrix(AA_seqvar_TM3, AA_features_TM3, categorized_matrix_TM3)
print(len(AA_mat_TM3[0]))
AA_mat_TM5 = ReadingFasta.makematrix(AA_seqvar_TM5, AA_features_TM5, categorized_matrix_TM5)
print(len(AA_mat_TM5[0]))
AA_mat_TM6 = ReadingFasta.makematrix(AA_seqvar_TM6, AA_features_TM6, categorized_matrix_TM6)
print(len(AA_mat_TM6[0]))
AA_mat_TM7 = ReadingFasta.makematrix(AA_seqvar_TM7, AA_features_TM7, categorized_matrix_TM7)
print(len(AA_mat_TM7[0]))

"""
#Creating output for 3Di sequences
# Read fasta file
fasta2 = open("../data_files/3DiSequences/fullset_ss.fasta")
#Create kmer frequency dictionary
seqvar2, features2 = ReadingFasta.make_seqvar(fasta2, di_seqs, di_features)
#Remove insignificant kmers
filter_feat2 = Filtering.richness_protein(features2, seqvar2, pos_counts, neg_counts)
# Make the matrix
Di_mat = ReadingFasta.makematrix(seqvar2, filter_feat2, di_matrix)
"""

#Create 3Di output for Tms 3,5,6,7
Di_dict = Globals.initialize_3Di_dict()
Di_seqvar_TM3, Di_features_TM3 = ReadingFasta.make_seqvar_TMS(Di_dict, 0, 5, di_seqs_TM3, di_features_TM3)
Di_seqvar_TM5, Di_features_TM5 = ReadingFasta.make_seqvar_TMS(Di_dict, 1, 5, di_seqs_TM5, di_features_TM5)
Di_seqvar_TM6, Di_features_TM6 = ReadingFasta.make_seqvar_TMS(Di_dict, 2, 5, di_seqs_TM6, di_features_TM6)
Di_seqvar_TM7, Di_features_TM7 = ReadingFasta.make_seqvar_TMS(Di_dict, 3, 5, di_seqs_TM7, di_features_TM7)

print(Di_features_TM3)
print(len(Di_features_TM3))
print(Di_features_TM5)
print(len(Di_features_TM5))
print(Di_features_TM6)
print(len(Di_features_TM6))
print(Di_features_TM7)
print(len(Di_features_TM7))

Di_mat_TM3 = ReadingFasta.makematrix(Di_seqvar_TM3, Di_features_TM3, di_matrix_TM3)
print(len(Di_mat_TM3[0]))
Di_mat_TM5 = ReadingFasta.makematrix(Di_seqvar_TM5, Di_features_TM5, di_matrix_TM5)
print(len(Di_mat_TM5[0]))
Di_mat_TM6 = ReadingFasta.makematrix(Di_seqvar_TM6, Di_features_TM6, di_matrix_TM6)
print(len(Di_mat_TM6[0]))
Di_mat_TM7 = ReadingFasta.makematrix(Di_seqvar_TM7, Di_features_TM7, di_matrix_TM7)
print(len(Di_mat_TM7[0]))

"""
#Concatenate AA and 3Di matrices
intermed_matrix = np.concatenate((np.array(AA_mat, dtype = np.uint8), np.array(Di_mat, dtype = np.uint8)) , axis = 1)
#Expand the protein matrix to account for the ligands
ligand_count = 55
proteins_matrix = np.repeat(intermed_matrix, repeats = ligand_count, axis = 0)

#Import dictionary matching ligands to SMILES String
ligand_dict = Globals.initialize_ligand_dict()
#Create ligands matrix
#ligand_matrix, ligand_features = SmileKmer.ligand_matrix(ligand_dict, 5, 1084)
ligand_matrix, ligand_features = SmileKmer.ligand_matrix(ligand_dict, 5, 100)

#Concatenate protein and ligand matrices
final_matrix = np.concatenate((proteins_matrix, np.array(ligand_matrix, dtype = np.uint8)), axis = 1)

#Create Classification Vector
proteins = seqvar1
logFCmat = []
for protein in proteins:
    for ligand in list(ligand_dict.keys()):
        logFCmat.append(float(classified[str(protein.name)][ligand]))

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

#The following code checks whether all entries are unique
#print(uniquematrix(final_matrix))

def import_final():
    #For No3Di.py
    global AA
    AA = AA_mat
    #For OneLigandRF.py and OneProteinRF.py
    global proteins
    proteins = seqvar1
    global dictionary
    dictionary = classified
    global protmat
    protmat = intermed_matrix
    #For TrainandTest.py
    global X
    X = final_matrix
    global Y
    Y = logFCmat
    global feat1
    feat1 = filter_feat
    global feat2
    feat2 = filter_feat2
    global feat3
    feat3 = ligand_features
"""
