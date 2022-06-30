#This script creates functions to import lists of ligands and proteins and match ligands with SMILES strings

#Imports
import pandas as pd

#Function to create a dictionary of ligands matched to SMILES strings
def initialize_ligand_dict():
    ligand_dict = {}
    df = pd.read_csv('../Ligands_withSMILE/ligand_SMILEs.csv')
    files = df['ligand file'].tolist()
    smiles = df['SMILE'].tolist()
    for i in range(len(files)):
        ligand_dict[files[i]] = smiles[i]

    return ligand_dict

#List of filenames for ligands that we have matched SMILES strings to
#Cannot distinguish chirality
def initialize_ligand_list():
    ligands = ['pS6_DE_1p_citronellol.csv', 'pS6_DE_1p_isoamylAcetate.csv', 'pS6_DE_1p_ethylTiglate.csv',
               'pS6_DE_1p_bIonone.csv', 'pS6_DE_1p_butyricAcid.csv',
               'pS6_DE_1p_paraCresol.csv', 'pS6_DE_1p_bCaryophyllene.csv', 'pS6_DE_p1_isovalericAcid.csv',
               'pS6_DE_1p_Octanal.csv', 'pS6_DE_1p_heptanal.csv',
               'pS6_DE_1p_tbm.csv', 'pS6_DE_1p_bDamascone.csv', 'pS6_DE_1p_pyridine.csv', 'pS6_DE_1p_propionicAcid.csv',
               'pS6_DE_1p_methylSalicylate.csv',
               'pS6_DE_p01_e2butene1thiol.csv', 'pS6_DE_1p_3methyl1butanethiol.csv', 'pS6_DE_1p_ethylButyrate.csv',
               'pS6_DE_1p_hexylTiglate.csv', 'pS6_DE_1p_indole.csv',
               'pS6_DE_500mM_2propylthietane.csv', 'pS6_DE_1p_2heptanone.csv', 'pS6_DE_p01_cyclopentanethiol.csv',
               'pS6_DE_1p_dimethyltrisulfide.csv',
               'pS6_DE_1p_guaiacol.csv', 'pS6_DE_1p_Benzaldehyde.csv', 'pS6_DE_p01_citral.csv',
               'pS6_DE_3mM_androstenone.csv', 'pS6_DE_100p_ebFarnesene.csv',
               'pS6_DE_1p_acetophenone.csv', 'pS6_DE_1p_transCinnamaldehyde.csv', 'pS6_DE_1p_linalool.csv',
               'pS6_DE_1p_2hexanone.csv', 'pS6_DE_1p_isopropylTiglate.csv',
               'pS6_DE_1p_aPinene.csv', 'pS6_DE_1p_diacetyl.csv', 'pS6_DE_1p_geranoil.csv',
               'pS6_DE_1p_heptanoicAcid.csv', 'pS6_DE_1p_2e3mp.csv', 'pS6_DE_1p_2hac.csv', 'pS6_DE_1p_2m2p.csv', 
               'pS6_DE_1p_2m2t.csv', 'pS6_DE_1p_2phenylAlcohol.csv', 'pS6_DE_1p_4methylAC.csv', 'pS6_DE_1p_25dmp.csv',
               'pS6_DE_1p_nCarvone.csv', 'pS6_DE_1p_nDihydrocarveol.csv', 'pS6_DE_1p_nMenthol.csv', 'pS6_DE_1p_ntmt.csv', 
               'pS6_DE_1p_pCarvone.csv', 'pS6_DE_1p_pDihydrocarveol.csv', 'pS6_DE_1p_pLimonene.csv', 'pS6_DE_1p_pMenthol.csv', 'pS6_DE_1p_sbt.csv', 'pS6_DE_1p_tmt.csv']
    return ligands

#Function to create list of protein accessions
def initialize_protein_list():
    df = pd.read_csv("TMs.csv")
    protein_list = list(df.iloc[:, 0])
    """
    acc_ids = []
    fr = open("../data_files/AminoAcidSequences/allsequences.fasta", "r")
    lines = fr.readlines()
    for line in lines:
        if line[0] == ">":
            acc_ids.append(line[1:-1])
    fr.close()
    """
    return protein_list


def initialize_AA_dict():
    df = pd.read_csv("TMs.csv")
    protein_list = initialize_protein_list()

    TMs_by_id = {}

    for i in range(len(protein_list)):
        TMs = [str(df.iloc[i, 1]), str(df.iloc[i, 4]), str(df.iloc[i, 7]), str(df.iloc[i, 10])]
        TMs_by_id[protein_list[i]] = TMs

    return categorize(TMs_by_id)

def initialize_indices():
    df = pd.read_csv("TMs.csv")
    protein_list = initialize_protein_list()

    TM_indices = {}
    for i in range(len(protein_list)):
        indices = [int(df.iloc[i, 2]), int(df.iloc[i, 3]), int(df.iloc[i, 5]), int(df.iloc[i, 6]), int(df.iloc[i, 8]),
                   int(df.iloc[i, 9]), int(df.iloc[i, 11]), int(df.iloc[i, 12]), ]
        TM_indices[protein_list[i]] = indices

    return TM_indices

def initialize_3Di_dict():
    TM_indices = initialize_indices()
    Di_dict = {}
    Di = open("../data_files/3DiSequences/fullset_ss.fasta", "r")
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
