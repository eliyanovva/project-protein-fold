# This script generates dictionaries of logFC, p-value, and classification with protein/ligand pairs as keys

import pandas as pd
import Globals

# Generate lists of proteins and ligands
acc_ids = Globals.initialize_protein_list()
csvs = Globals.initialize_ligand_list()

# Initialize variables
num_proteins = len(acc_ids)
num_ligands = len(csvs)
logFC_byID = {}
FDR_byID = {}
cit_logFC = {}
cit_FDR = {}


# returns dictionaries of the logFC and p-values
# key: protein id, value: dict (key: ligand file name, value: data label)
def labels():
    # id = accession number of a given protein
    for id in acc_ids:
        logFC_byID[id] = {}
        FDR_byID[id] = {}

    fas_df = pd.read_csv('uniprot_ensemble.csv', index_col='accession number')

    # Read each csv file for the corresponding ligand
    for csv in csvs:
        file_name = '../olfr_de/' + csv
        curr_df = pd.read_csv(file_name, index_col='ensembl_gene_id')

        for id in acc_ids:
            ensem_id = fas_df.loc[id]['ensembl_gene_id']  # The ENSEMBLE id corresponding to the accession number
            logFC_byID[id][csv] = (curr_df.loc[ensem_id]['logFC'])  # Find logFC for the ligand-protein pair
            FDR_byID[id][csv] = (curr_df.loc[ensem_id]['FDR'])  # Find p-value for the ligand-protein pair

    # Return dictionaries with protein-ligand pair keys and logFC and p-value values
    return logFC_byID, FDR_byID


# Create a classification dictionary with protein-ligand pair keys and bind (1) or not bind (0) as values
def classified_logFC_FDR(logFC_byID, FDR_byID, protein_list):
    classified = {}
    pos_counts = {}  # key: protein id, value: number of positive protein interactions
    neg_counts = {}  # key: protein id, value: number of negative protein interactions

    pos_pairs = []
    neg_pairs = []

    class_by_CSV = {}

    for csv in csvs:
        class_by_CSV[csv] = 0

    count = 0

    for id in protein_list:
        pos = 0
        neg = 0
        classified[id] = {}
        for csv in csvs:
            if FDR_byID[id][csv] <= .1:
                count += 1

                if logFC_byID[id][csv] >= 1:  # The protein and ligand bind
                    classified[id][csv] = 1
                    pos += 1    #only update pos count if pair isn't removed
                    class_by_CSV[csv] += 1
                    pos_pairs.append([id, csv])
                elif logFC_byID[id][csv] < 1:  # The protein and ligand do not bind
                    classified[id][csv] = 0
                    neg += 1    #only update neg count if pair isn't removed
                    neg_pairs.append([id, csv])

                pos_counts[id] = pos
                neg_counts[id] = neg

    return classified, pos_counts, neg_counts, pos_pairs, neg_pairs