from tkinter import Y
import PreparingMatrix
import SmileKmer
import numpy as np
import ReadingFasta
import labels
import SMILE

ligand_dict = SMILE.create_ligand_dict()
cit_logFC, cit_pval = labels.cit_labels()
logFC, pVal = labels.labels()

#Import proteins matrix
PreparingMatrix.access_matrix()
proteins_matrix = PreparingMatrix.intermediate_matrix

#Import ligands matrix
SmileKmer.importmatrix(ligand_dict, 5, 29)
ligand_matrix = SmileKmer.ligmat
#print(ligand_matrix)

#Import dictionary
PreparingMatrix.access_dictionary()
logFC_dict = PreparingMatrix.dictionary
#print(logFC_dict)

#Concatenate protein and ligand matrices
#print(len(proteins_matrix))
#print(len(proteins_matrix[0]))
#print(len(ligand_matrix))
#print(len(ligand_matrix[0]))
final_matrix = np.concatenate((proteins_matrix, ligand_matrix), axis = 1)
#print(final_matrix)

#Create logFC vector
ReadingFasta.import_variables()
proteins = ReadingFasta.sequence_seqs
logFC = []
for protein in proteins:
    if protein.name == 'A0A1L1SRP1':
        break
    for ligand in list(ligand_dict.keys()):
        logFC.append(float(logFC_dict[str(protein.name)][ligand]))
#print(len(final_matrix))
#print(len(logFC))

"""#logFC vector
ReadingFasta.import_variables()
proteins = ReadingFasta.sequence_seqs
logFC = []
for protein in proteins:
    for ligand in list(ligand_dict.keys()):
        logFC.append(logFC_dict[frozenset([protein.name,ligand])])"""

def import_final():
    global X
    X = final_matrix
    global Y
    Y = logFC
