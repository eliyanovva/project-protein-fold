import ReadingFasta
import pandas as pd

#Load in the matrices
ReadingFasta.import_variables()

#Convert csv to dictionary
dict = {}
data = pd.read_csv('/home/users/sml96/bin/project-protein-fold/olfr_de/uniprot_ligand_logfc.csv')
