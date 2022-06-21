from tkinter import Y

from cv2 import log
import PreparingMatrix
import SmileKmer
import numpy as np
import ReadingFasta
import labels
import SMILE
import Filtering

ligand_dict = SMILE.create_ligand_dict()
cit_logFC, cit_pval = labels.cit_labels()
logFC, pVal = labels.labels()
classified, pos_counts, neg_counts = labels.classified_logFC_pVal()

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
PreparingMatrix.access_matrix()
proteins_matrix = PreparingMatrix.intermediate_matrix

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

freq_dict = {}
kmers = []
Filtering.richness_ligand(kmers, freq_dict, pos_counts, neg_counts)

def import_final():
    global X
    X = final_matrix
    global Y
    Y = logFCmat