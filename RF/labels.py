import pandas as pd
import numpy as np
import Globals

acc_ids = Globals.initialize_protein_list()
csvs = Globals.initialize_ligand_list()

fas_df = pd.read_csv('uniprot_ensemble.csv', index_col='accession number')
#fas_df = pd.read_csv('fasta_list.csv', index_col='accession number')

num_proteins = len(acc_ids)
num_ligands = len(csvs)
logFC_byID = {}
pVal_byID = {}
cit_logFC = {}
cit_pVal = {}

#id = accession number of a given protein
for id in acc_ids:
           logFC_byID[id] = {}
           pVal_byID[id] = {}

#returns vectors of the logFC and p-values
#key: protein id, value: dict (key: ligand file name, value: data label)
def labels():
    for csv in csvs:
        file_name = 'olfr_de_copy1/olfr_de/'+csv
        curr_df = pd.read_csv(file_name, index_col='ensembl_gene_id')
        #curr_df = pd.read_csv(file_name, index_col='name')

        for id in acc_ids:

            ensem_id = fas_df.loc[id]['ensembl_gene_id']
            logFC_byID[id][csv] = (curr_df.loc[ensem_id]['logFC'])
            pVal_byID[id][csv] = (curr_df.loc[ensem_id]['PValue'])
            """
            name = fas_df.loc[id]['receptor']
            logFC_byID[id][csv] = (curr_df.loc[name]['logFC'])
            pVal_byID[id][csv] = (curr_df.loc[name]['PValue'])
            """
    return logFC_byID, pVal_byID

def classified_logFC_pVal(logFC_byID, pVal_byID):
    classified = {}
    pos_counts = {} #key: protein id, value: number of positive protein interactions
    neg_counts = {} #key: protein id, value: number of negative protein interactions

    for id in pVal_byID:
        pos = 0
        neg = 0
        classified[id] = {}
        for csv in csvs:
            if (logFC_byID[id][csv] >= .5) & (pVal_byID[id][csv] <= .05):
                classified[id][csv] = 1
                pos += 1
            else:
                classified[id][csv] = 0
                neg += 1
        pos_counts[id] = pos
        neg_counts[id] = neg

    return classified, pos_counts, neg_counts

#pos_counts[id] = num. of pos interactions that protein id has