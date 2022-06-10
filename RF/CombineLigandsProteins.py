import PreparingMatrix
import SmileKmer


#Import proteins matrix
PreparingMatrix.access_matrix()
proteins_matrix = PreparingMatrix.intermediate_matrix

#Import ligands matrix
SmileKmer.importmatrix()
ligand_matrix = SmileKmer.ligmat

#Import dictionary
PreparingMatrix.access_dictionary()
logFC_dict = PreparingMatrix.dictionary
