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

for i in range(10):
    print(intermed_matrix[i])

#Import ligands matrix
SmileKmer.importmatrix(ligand_dict, 5, 1084)
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

def import_final():
    global X
    X = final_matrix
    global Y
    Y = logFCmat