import ReadingFasta
import pandas as pd

#Load in the matrices
ReadingFasta.import_variables()

#Convert csv to dictionary
dict = {}
df = pd.read_csv('/home/users/sml96/bin/project-protein-fold/olfr_de/uniprot_ligand_logfc.csv', header = None)

for i in range(0, 35477, 2):
    dict[frozenset(df.loc[i])] = df.loc[i+1][1]

print(dict)