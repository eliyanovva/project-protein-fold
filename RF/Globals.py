#This script creates functions to import lists of ligands and proteins and match ligands with SMILES strings

#Methodology for categorization is based on 'String-Based Models for Predicting RNA-Protein Interaction',
#https://dl.acm.org/doi/10.1145/3107411.3107508

#Imports
import pandas as pd

#Additional coding help from:
#https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.shape.html
#https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.at.html
#https://stackoverflow.com/questions/23748995/pandas-dataframe-column-to-list
#https://www.geeksforgeeks.org/writing-to-file-in-python/
#https://www.geeksforgeeks.org/how-to-read-from-a-file-in-python/
#https://www.w3schools.com/python/ref_string_replace.asp

def initialize_ligand_dict(smile_location):
    """
    This function creates a dictionary mapping ligands to their SMILE formulas.

    Args:
        smile_location (str): file path to data table of ligands and their SMILE formulas

    Returns:
        ligand_dict (dict): dictionary mapping a ligand to its SMILE formula
    """
    ligand_dict = {}
    df = pd.read_csv(smile_location)
    files = initialize_ligand_list(smile_location)
    smiles = df['SMILE'].tolist()
    for i in range(len(files)):
        ligand_dict[files[i]] = smiles[i]

    return ligand_dict

def initialize_ligand_list(smile_location):
    """
    This function creates a list of every ligand that has a SMILE formula.
    The current set of SMILE formulas cannot distinguish chirality

    Args:
        smile_location (str): file path to data table of ligands and their SMILE formulas

    Returns:
        ligands (list): list of every ligand that has a SMILE formula.
    """
    df = pd.read_csv(smile_location)
    ligands = df['Ligands'].tolist()

    return ligands

def initialize_protein_list(TM_location):
    """
    This function returns a list of proteins that have sequence data for TMs 3,5,6, and 7.

    Args:
        TM_location (str): file path to data table of proteins and their sequences for TMs 3,5,6 and 7

    Returns:
        protein_list (list): list of proteins that have sequence data for TMs 3,5,6, and 7
    """
    df = pd.read_csv(TM_location)
    proteins = list(df.iloc[:, 0])

    return proteins


def initialize_AA_dict(proteins, TM_location):
    """
    This function returns a dictionary mapping protein ids to the amino acid sequences for TMs 3,5,6, and 7.

    Args:
        proteins (list): list of proteins to extract TM data for
        TM_location (str): file path to data table of proteins and their data for TMs 3,5,6 and 7

    Returns:
        AA_dict (dict): dictionary mapping a protein id to a list of amino acid sequences for TMs 3,5,6, and 7
            ex: AA_dict[id] = [AA seq for id's TM3, AA seq for id's TM5, AA seq for id's TM6, AA seq for id's TM7]
    """
    df = pd.read_csv(TM_location)

    TMs_by_id = {}
    num_rows = df.shape[0]
    for i in range(num_rows):
        id = df.at[i, 'protein']
        if id in proteins:
            TMs = [str(df.at[i, 'TM3']), str(df.at[i, 'TM5']), str(df.at[i, 'TM6']), str(df.at[i, 'TM7'])]
            TMs_by_id[id] = TMs

    AA_dict = categorize(TMs_by_id)

    return AA_dict

def initialize_indices(proteins, TM_location):
    """
    This function returns a dictionary mapping protein ids to the start/stop indices for TMs 3,5,6, and 7.

    Args:
        proteins (list): list of proteins to extract TM data for
        TM_location (str): file path to data table of proteins and their data for TMs 3,5,6 and 7

    Returns:
        TM_indices (dict): dictionary mapping a protein id to the start/stop indices for TMs 3,5,6, and 7
            ex: TM_indices[id] = [start index of id's TM3, stop index of id's TM3, start index of id's TM5,
                                    stop index of id's TM5, start index of id's TM6, stop index of id's TM6,
                                    start index of id's TM7, stop index of id's TM7]
    """
    df = pd.read_csv(TM_location)

    TM_indices = {}
    num_rows = df.shape[0]
    for i in range(num_rows):
        id = df.at[i, 'protein']
        if id in proteins:
            indices = [int(df.at[i,'s3']), int(df.at[i,'e3']), int(df.at[i,'s5']), int(df.at[i,'e5']),
                       int(df.at[i,'s6']), int(df.at[i,'e6']), int(df.at[i,'s7']), int(df.at[i,'e7'])]
            TM_indices[id] = indices

    return TM_indices

def initialize_3Di_dict(proteins, TM_location, Di_location):
    """
    This function returns a dictionary mapping protein ids to the 3di sequences for TMs 3,5,6, and 7.

    Args:
        proteins (list): list of proteins to extract TM data for
        TM_location (str): file path to data table of proteins and their data for TMs 3,5,6 and 7
        Di_location (str): file path to fasta of proteins and their 3di sequences

    Returns:
        Di_TMs (dict): dictionary mapping protein ids to the 3di sequences for TMs 3,5,6, and 7.
            ex: Di_TMs[id] = [3di seq for id's TM3, 3di seq for id's TM5, 3di seq for id's TM6, 3di seq for id's TM7]
    """
    TM_indices = initialize_indices(proteins, TM_location)
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

def categorize(AA_dict):
    """
    This function categorizes the amino acid sequences in AA_dict.

    Args:
        AA_dict (dict): dictionary mapping a protein id to a list of amino acid sequences for TMs 3,5,6, and 7
            ex: AA_dict[id] = [AA seq for id's TM3, AA seq for id's TM5, AA seq for id's TM6, AA seq for id's TM7]

    Returns:
        Updates the AA sequences stored in AA_dict.
        For each sequence in AA_dict[id], the amino acids are categorized based on their
            dipole moments and side chain volumes.
    """
    categorize_dict = {}
    for id in AA_dict:
        categorize_TMs = []
        for TM in AA_dict[id]:
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