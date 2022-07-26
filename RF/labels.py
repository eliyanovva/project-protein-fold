# This script generates dictionaries of logFC, FDR, and classification with protein/ligand pairs as keys

import pandas as pd
import Globals

# Generate lists of proteins and ligands
TM_location = "../data_files/TMdomains/TM.csv"
smile_location = "../Ligands_withSMILE/ligand_SMILES.csv"
acc_ids = Globals.initialize_protein_list(TM_location)
ligands = Globals.initialize_ligand_list(smile_location)

def labels(ligand_folder):
    """
    This function extracts the experimental logFC and FDR values for the protein-ligand pairs.

    Args:
        ligand_folder (str): path file to folder of binding experiment datafiles

    Returns:
        logFC_byID (dict): dictionary mapping a protein-ligand pair to its logFC value
            ex: logFC_byID[id][lig] = logFC of the protein-ligand pair of id-lig
        FDR_byID (dict): dictionary mapping a protein-ligand pair to its FDR value
            ex: logFC_byID[id][lig] = FDR of the protein-ligand pair of id-lig
    """
    # Initialize variables
    logFC_byID = {}
    FDR_byID = {}
    # id = accession number of a given protein
    for id in acc_ids:
        logFC_byID[id] = {}
        FDR_byID[id] = {}

    fas_df = pd.read_csv('uniprot_ensemble.csv', index_col='accession number')

    # Read each csv file for the corresponding ligand
    for lig in ligands:
        file_name = ligand_folder + lig
        curr_df = pd.read_csv(file_name, index_col='ensembl_gene_id')

        for id in acc_ids:
            ensem_id = fas_df.loc[id]['ensembl_gene_id']  # The ENSEMBLE id corresponding to the accession number

            logFC_byID[id][lig] = (curr_df.loc[ensem_id]['logFC'])  # Find logFC for the ligand-protein pair
            FDR_byID[id][lig] = (curr_df.loc[ensem_id]['FDR'])  # Find FDR for the ligand-protein pair

    # Return dictionaries with protein-ligand pair keys and logFC and FDR values
    return logFC_byID, FDR_byID

def extract_new_combos(FDR_byID, proteins, ligands):
    new_combos = {}
    i = 0
    for id in proteins:
        new_combos[id] = []
        for lig in ligands:
            FDR = FDR_byID[id][lig]
            if (FDR > .1) & (FDR <= .4):
                new_combos[id].append(lig)
                i += 1
    return new_combos

# Create a classification dictionary with protein-ligand pair keys and bind (1) or not bind (0) as values
def classified_logFC_FDR(logFC_byID, FDR_byID, protein_list):
    """
    This function classifies protein-ligand pairs as to whether or not they bind with each other.

    Args:
        logFC_byID (dict): dictionary mapping a protein-ligand pair to its logFC value
            ex: logFC_byID[id][lig] = logFC of the protein-ligand pair of id-lig
        FDR_byID (dict): dictionary mapping a protein-ligand pair to its FDR value
            ex: logFC_byID[id][lig] = FDR of the protein-ligand pair of id-lig
        protein_list (list): proteins from all experimental protein-ligand pairs

    Returns:
        classified (dict): dictionary mapping a protein-ligand to their classification
            ex: classified[id][lig] = {1 if id and lig bind, 0 if they do not}
        pos_counts (dict): dictionary mapping a protein id to the # of positive protein-ligand pairs with id as the protein
            ex: If the protein id only binds with the ligands L1, L2, and L3, then pos_counts[id] = 3
        neg_counts (dict): dictionary mapping a protein id to the # of negative protein-ligand pairs with id as the protein
            ex: If the protein id doesn't bind with the ligands L1 and L4, then neg_counts[id] = 2
        pos_dict (dict): dictionary mapping a protein id to the list of ligands that id binds with
            ex: If the protein id only binds with the ligands L1, L2, and L3, then pos_counts[id] = [L1, L2, L3]
        neg_dict (dict): dictionary mapping a protein id to the list of ligands that id does not bind with
            ex: If the protein id doesn't bind with the ligands L1 and L4, then neg_counts[id] = [L1, L4]
        proteins_toconsider (list): sorted list of proteins that have at least 1 positive or negative interaction with a ligand
    """
    classified = {}
    pos_counts = {}
    neg_counts = {}
    pos_dict = {}
    neg_dict = {}
    proteins_toconsider = set()

    for id in protein_list:
        pos = 0         #running count of pos pairs with id
        neg = 0         #running count of neg pairs with id
        classified[id] = {}
        for lig in ligands:
            if FDR_byID[id][lig] <= .1:
                if logFC_byID[id][lig] >= 1:    # The protein and ligand bind
                    classified[id][lig] = 1
                    pos += 1                    #only update pos count if pair isn't removed

                    if id not in pos_dict:
                        pos_dict[id] = []
                    pos_dict[id].append(lig)

                    proteins_toconsider.add(id)

                elif logFC_byID[id][lig] < 1:   # The protein and ligand do not bind
                    classified[id][lig] = 0
                    neg += 1                    #only update neg count if pair isn't removed

                    if id not in neg_dict:
                        neg_dict[id] = []
                    neg_dict[id].append(lig)

                    proteins_toconsider.add(id)

                pos_counts[id] = pos
                neg_counts[id] = neg

    proteins_toconsider = list(proteins_toconsider)
    proteins_toconsider.sort()

    return classified, pos_counts, neg_counts, pos_dict, neg_dict, proteins_toconsider

#labels('../olfr_de/')