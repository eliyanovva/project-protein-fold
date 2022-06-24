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
               'pS6_DE_1p_heptanoicAcid.csv']
    return ligands

#Function to create list of protein accessions
def initialize_protein_list():
    acc_ids = []
    fr = open("../AminoAcidSequences/allsequences.fasta", "r")
    lines = fr.readlines()
    for line in lines:
        if line[0] == ">":
            acc_ids.append(line[1:-1])
    fr.close()

    return acc_ids

