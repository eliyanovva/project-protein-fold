#This script creates functions to import lists of ligands and proteins and match ligands with SMILES strings

#Imports
import pandas as pd

#Additional coding help from:
#https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.shape.html
#https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.at.html
#https://stackoverflow.com/questions/23748995/pandas-dataframe-column-to-list
#https://www.geeksforgeeks.org/writing-to-file-in-python/
#https://www.geeksforgeeks.org/how-to-read-from-a-file-in-python/
#https://www.w3schools.com/python/ref_string_replace.asp


#Function to create a dictionary of ligands matched to SMILES strings
def initialize_ligand_dict(smile_location):
    ligand_dict = {}
    df = pd.read_csv(smile_location)
    files = initialize_ligand_list(smile_location)
    smiles = df['SMILE'].tolist()
    for i in range(len(files)):
        ligand_dict[files[i]] = smiles[i]

    return ligand_dict

#List of filenames for ligands that we have matched SMILES strings to
#Cannot distinguish chirality
def initialize_ligand_list(smile_location):
    df = pd.read_csv(smile_location)
    ligands = df['Ligands'].tolist()

    return ligands

#Function to create list of protein accessions
def initialize_protein_list(TM_location):
    df = pd.read_csv(TM_location)
    protein_list = list(df.iloc[:, 0])

    return protein_list


def initialize_AA_dict(proteins, TM_csv):
    df = pd.read_csv(TM_csv)
    print(df)
    TMs_by_id = {}
    num_rows = df.shape[0]
    for i in range(num_rows):
        id = df.at[i, 'protein']
        if id in proteins:
            TMs = [str(df.at[i, 'TM3']), str(df.at[i, 'TM5']), str(df.at[i, 'TM6']), str(df.at[i, 'TM7'])]
            TMs_by_id[id] = TMs

    return categorize(TMs_by_id)

def initialize_indices(proteins, TM_csv):
    df = pd.read_csv(TM_csv)

    TM_indices = {}
    num_rows = df.shape[0]
    for i in range(num_rows):
        id = df.at[i, 'protein']
        if id in proteins:
            indices = [int(df.at[i,'s3']), int(df.at[i,'e3']), int(df.at[i,'s5']), int(df.at[i,'e5']),
                       int(df.at[i,'s6']), int(df.at[i,'e6']), int(df.at[i,'s7']), int(df.at[i,'e7'])]
            TM_indices[id] = indices

    return TM_indices

def initialize_3Di_dict(p_list, TM_location, Di_location):
    TM_indices = initialize_indices(p_list, TM_location)
    Di_dict = {}
    Di = open(Di_location, "r")
    lines = Di.readlines()

    for i in range(len(lines)):
        if i % 2 == 0:
            id = lines[i][1:-1]
            seq = lines[i + 1][:-1]
            Di_dict[id] = seq
    Di.close()

    Di_TMs = {}
    for id in TM_indices:
        TMs = []
        seq = Di_dict[id]

        for i in range(4):
            start = TM_indices[id][2 * i]
            end = TM_indices[id][(2 * i) + 1]
            TMs.append(seq[start - 1:end])

        Di_TMs[id] = TMs

    return Di_TMs

def categorize(TM_dict):
    categorize_dict = {}
    for id in TM_dict:
        categorize_TMs = []
        for TM in TM_dict[id]:
            TM = TM.replace('A', 'a').replace('G', 'a').replace('V', 'a')
            TM = TM.replace('I', 'b').replace('L', 'b').replace('F', 'b').replace('P', 'b')
            TM = TM.replace('Y', 'c').replace('M', 'c').replace('T', 'c').replace('S', 'c')
            TM = TM.replace('H', 'd').replace('N', 'd').replace('Q', 'd').replace('W', 'd')
            TM = TM.replace('R', 'e').replace('K', 'e')
            TM = TM.replace('D', 'f').replace('E', 'f')
            TM = TM.replace('C', 'g')
            categorize_TMs.append(TM)
        categorize_dict[id] = categorize_TMs
    return categorize_dict