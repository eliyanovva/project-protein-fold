from cmath import nan
from sklearn.decomposition import DictionaryLearning
import ReadingFasta
import pandas as pd
import numpy as np

#Load in the matrices
ReadingFasta.import_variables()

#Convert csv to dictionary of proteins with nested dictionary of ligands to logFC
protdict = {}
df = pd.read_csv('/home/users/sml96/bin/project-protein-fold/olfr_de/uniprot_ligand_logfc.csv', header = None)

for i in range (0, 35476, 2):
    if str(df.loc[i][0]) not in protdict:
        protdict[str(df.loc[i][0])] = {}    
    ligdict = protdict[str(df.loc[i][0])]
    ligdict[str(df.loc[i][1])] = str(df.loc[i+1][1])
    """if (str(df.loc[i][0]) == 'Q8VEZ3'):
        print(ligdict)
        print(protdict[str(df.loc[i][0])])
        print(str(df.loc[i][0]))"""
    protdict[str(df.loc[i][0])] = ligdict
    """if (str(df.loc[i][0]) == 'Q8VEZ3'):
        print(protdict[str(df.loc[i][0])])"""

"""i = 0
while i < 35477:
    if (df.loc[i][0] == 'Q8VEZ3'):
        print(i)
    i += 2"""


#See dictionary
"""print(protdict.keys())
print(protdict['Q8VEZ3'])"""

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

def access_dictionary():
    global dictionary
    dictionary = protdict