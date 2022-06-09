import ReadingFasta
import pandas as pd
import numpy as np

#Load in the matrices
ReadingFasta.import_variables()

#Convert csv to dictionary
dict = {}
df = pd.read_csv('/home/users/sml96/bin/project-protein-fold/olfr_de/uniprot_ligand_logfc.csv', header = None)

for i in range(0, 35477, 2):
    dict[frozenset(df.loc[i])] = df.loc[i+1][1]

#Concatenate AA and 3Di sequences
protein_matrix = np.concatenate((np.array(ReadingFasta.sequence_matrix), np.array(ReadingFasta.structure_matrix)) , axis = 1)


def expand(matrix, ligand_count):
    return np.repeat(matrix, repeats = ligand_count, axis = 0)

expanded_matrix = expand(protein_matrix, 39)

#Check for proper expansion
#print(len(protein_matrix))
#print(len(expanded_matrix))

#Allow access of the expanded matrix in another script
def access_matrix():
    global intermediate_matrix
    intermediate_matrix = expanded_matrix