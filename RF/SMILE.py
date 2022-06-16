import pandas as pd

ligand_dict = {}
df = pd.read_csv('ligand_SMILEs.csv')

def create_ligand_dict():
    files = df['ligand file'].tolist()
    smiles = df['SMILE'].tolist()
    for i in range(len(files)):
        ligand_dict[files[i]] = smiles[i]

    return ligand_dict
