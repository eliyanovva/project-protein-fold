import SmileKmer
import numpy as np
import ReadingFasta
import labels
import Globals
import Filtering

ligand_dict = Globals.initialize_ligand_dict()

logFC, pVal = labels.labels()
classified, pos_counts, neg_counts = labels.classified_logFC_pVal(logFC, pVal)

def import_labels():
    global positives
    positives = pos_counts
    global negatives
    negatives = neg_counts

def exportdicts():
    global logdic
    logdic = logFC
    global pdic
    pdic = pVal
    global class_dict
    class_dict = classified

#Import proteins matrix
#PreparingMatrix.access_matrix()
#proteins_matrix = PreparingMatrix.intermediate_matrix

#Initialize Set of Features
#categorized variables
categorized_features = set()
categorized_seqs = []
categorized_matrix = []
#3Di variables
di_features = set()
di_seqs = []
di_matrix = []

#Creating output for categorized amino acids
#Read fasta file
fasta1 = open("../AminoAcidSequences/fully_categorized.fasta")
#fasta1 = open("../AminoAcidSequences/categorized.fasta")
#Create kmer frequency dictionary
seqvar1, features1 = ReadingFasta.make_seqvar(fasta1, categorized_seqs, categorized_features)
#Remove insignificant kmers
filter_feat = Filtering.richness_protein(features1, seqvar1, pos_counts, neg_counts)
# Make the matrix
AA_mat = ReadingFasta.makematrix(seqvar1, filter_feat, categorized_matrix)
#AA_mat = ReadingFasta.makematrix(seqvar1, features1, Globals.categorized_matrix)

#Creating output for 3Di sequences
# Read fasta file
fasta2 = open("../3DiSequences/fullset_ss.fasta")
#fasta2 = open("../3DiSequences/outputDb_ss.fasta")
#Create kmer frequency dictionary
seqvar2, features2 = ReadingFasta.make_seqvar(fasta2, di_seqs, di_features)
#Remove insignificant kmers
filter_feat2 = Filtering.richness_protein(features2, seqvar2, pos_counts, neg_counts)
# Make the matrix
Di_mat = ReadingFasta.makematrix(seqvar2, filter_feat2, di_matrix)
#Di_mat = ReadingFasta.makematrix(seqvar2, features2, Globals.di_matrix)

intermed_matrix = np.concatenate((np.array(AA_mat, dtype = np.uint8), np.array(Di_mat, dtype = np.uint8)) , axis = 1)
ligand_count = 38
proteins_matrix = np.repeat(intermed_matrix, repeats = ligand_count, axis = 0)

print(len(AA_mat))
print(len(AA_mat[0]))
print(len(Di_mat))
print(len(Di_mat[0]))



#Import ligands matrix
SmileKmer.importmatrix(ligand_dict, 5, 1084)
#SmileKmer.importmatrix(ligand_dict, 5, 230)
ligand_matrix = SmileKmer.ligmat

print("Num. Ligand kmers: " + str(len(ligand_matrix[0])))

#Concatenate protein and ligand matrices
final_matrix = np.concatenate((proteins_matrix, np.array(ligand_matrix, dtype = np.uint8)), axis = 1)

#Create logFC vector
proteins = seqvar1
logFCmat = []
for protein in proteins:
    for ligand in list(ligand_dict.keys()):
        logFCmat.append(float(classified[str(protein.name)][ligand]))

print(len(final_matrix))

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

print(uniquematrix(final_matrix))

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
    #For Regular Trials
    global X
    X = final_matrix
    global Y
    Y = logFCmat