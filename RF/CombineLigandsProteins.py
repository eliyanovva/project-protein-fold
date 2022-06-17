from tkinter import Y

from cv2 import log
import PreparingMatrix
import SmileKmer
import numpy as np
import ReadingFasta
import labels
import SMILE

ligand_dict = SMILE.create_ligand_dict()
cit_logFC, cit_pval, cit_corrected = labels.cit_labels()
logFC, pVal = labels.labels()
classified = labels.classified_logFC_pVal()

def exportdicts():
    global citlog  
    citlog = cit_logFC
    global citp
    citp = cit_pval
    global citcor
    citcor = cit_corrected
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
i = 0
for protein in proteins:
    for ligand in list(ligand_dict.keys()):
        logFCmat.append(float(classified[str(protein.name)][ligand]))
        i+=1
        if (float(classified[str(protein.name)][ligand]) == 1.0):
            np.concatenate((final_matrix, np.reshape(final_matrix[i], (1,-1))), axis = 0)
            proteins.append(protein)
print(logFCmat[230:])


def import_final():
    global X
    X = final_matrix
    global Y
    Y = logFCmat
