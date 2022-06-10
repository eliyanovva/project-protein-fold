import PreparingMatrix
import SmileKmer


#Import proteins matrix
PreparingMatrix.access_matrix()
proteins_matrix = PreparingMatrix.intermediate_matrix

SmileKmer.importmatrix()
print(SmileKmer.ligmat)