from tkinter import Y

from cv2 import log
import PreparingMatrix
import SmileKmer
import numpy as np
import ReadingFasta
import labels
import SMILE
import Globals
import Filtering

ligand_dict = SMILE.create_ligand_dict()
cit_logFC, cit_pval = labels.cit_labels()

logFC, pVal = labels.labels()
classified, pos_counts, neg_counts = labels.classified_logFC_pVal(logFC, pVal)

def import_labels():
    global positives
    positives = pos_counts
    global negatives
    negatives = neg_counts

def exportdicts():
    global citlog
    citlog = cit_logFC
    global citp
    citp = cit_pval
    #global citcor
    #citcor = cit_corrected
    global logdic
    logdic = logFC
    global pdic
    pdic = pVal
    global class_dict
    class_dict = classified

#Import proteins matrix
#PreparingMatrix.access_matrix()
#proteins_matrix = PreparingMatrix.intermediate_matrix

#Creating output for categorized amino acids
#Read fasta file
fasta1 = open("/home/users/sml96/bin/project-protein-fold/AminoAcidSequences/categorized.fasta")
#Create kmer frequency dictionary
seqvar1, features1 = ReadingFasta.make_seqvar(fasta1, Globals.categorized_seqs, Globals.categorized_features)
#Remove insignificant kmers
filter_feat = Filtering.richness_protein(features1, seqvar1, pos_counts, neg_counts)
# Make the matrix
AA_mat = ReadingFasta.makematrix(seqvar1, filter_feat, Globals.categorized_matrix)

#Creating output for 3Di sequences
# Read fasta file
fasta2 = open("/home/users/sml96/bin/project-protein-fold/foldseek-master/foldseek/foldseek/outputDb_ss.fasta")
#Create kmer frequency dictionary
seqvar2, features2 = ReadingFasta.make_seqvar(fasta2, Globals.di_seqs, Globals.di_features)
#Remove insignificant kmers
filter_feat2 = Filtering.richness_protein(features2, seqvar2, pos_counts, neg_counts)
# Make the matrix
Di_mat = ReadingFasta.makematrix(seqvar2, filter_feat, Globals.di_matrix)

intermed_matrix = np.concatenate((np.array(AA_mat), np.array(Di_mat)) , axis = 1)
ligand_count = 38
proteins_matrix = np.repeat(intermed_matrix, repeats = ligand_count, axis = 0)

#Import ligands matrix
SmileKmer.importmatrix(ligand_dict, 5, 230)
ligand_matrix = SmileKmer.ligmat

#Concatenate protein and ligand matrices
final_matrix = np.concatenate((proteins_matrix, ligand_matrix), axis = 1)

#Create logFC vector
ReadingFasta.import_variables()
proteins = ReadingFasta.sequence_seqs
logFCmat = []
for protein in proteins:
    for ligand in list(ligand_dict.keys()):
        logFCmat.append(float(classified[str(protein.name)][ligand]))

def import_final():
    global X
    X = final_matrix
    global Y
    Y = logFCmat