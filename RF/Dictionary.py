import pandas as pd

file = open('/home/users/sml96/bin/project-protein-fold/olfr_de/uniprot_ligand_logfc.csv')
reader = pd.read_csv(file, header = None)
dictionary = {}
print(reader)
for i in range(0, 35477, 2):
    dictionary[frozenset(reader.loc[i])] = reader.loc[i+1][1]

