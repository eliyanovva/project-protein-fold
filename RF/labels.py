import pandas as pd
import numpy as np
import SMILE

acc_ids = SMILE.create_protein_list()
csvs = SMILE.create_ligand_list()

#acc_ids = ["P1", "P2", "P3"]
#csvs = ["L1", "L2", "L3"]

fas_df = pd.read_csv('fasta_list.csv', index_col='accession number')

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
           
def cit_labels():
           cit_df = pd.read_csv('olfr_de_copy1/olfr_de/pS6_DE_1p_citronellol.csv', index_col='name')
           
           for id in acc_ids:
                      name = fas_df.loc[id]['receptor']
                      cit_logFC[id] = (cit_df.loc[name]['logFC'])
                      cit_pVal[id] = (cit_df.loc[name]['PValue'])
           return cit_logFC, cit_pVal

#returns vectors of the logFC and p-values
#key: protein id, value: dict (key: ligand file name, value: data label)
def labels():
           for csv in csvs:
                      file_name = 'olfr_de_copy1/olfr_de/'+csv
                      curr_df = pd.read_csv(file_name, index_col='name')

                      for id in acc_ids:
                                 name = fas_df.loc[id]['receptor']
                                 logFC_byID[id][csv] = (curr_df.loc[name]['logFC'])
                                 pVal_byID[id][csv] = (curr_df.loc[name]['PValue'])
           return logFC_byID, pVal_byID

#cutoff for p-Val: > .05
#cutoff for logFC: < .5
#1 = yes, 0 = no

def classified_logFC(logFC_byID):
    classified = np.zeros_like(np.arange(num_proteins * num_ligands))
    i = 0
    for id in logFC_byID:
        for csv in csvs:
            if logFC_byID[id][csv] < .5:
                classified[i] = 0
            else:
                classified[i] = 1
            i += 1

    return classified

def classified_pVal(pVal_byID):
    classified = np.zeros_like(np.arange(num_proteins * num_ligands))
    i = 0
    for id in pVal_byID:
        for csv in csvs:
            if pVal_byID[id][csv] > .05:
                classified[i] = 0
            else:
                classified[i] = 1
            i += 1

    return classified

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

#logFC_byID = {"P1": {"L1": .8, "L2": .7, "L3": .6}, "P2": {"L1": .7, "L2": .2, "L3": .6}, "P3": {"L1": .3, "L2": .4, "L3": .6}}
#pVal_byID = {"P1": {"L1": .01, "L2": .04, "L3": .02}, "P2": {"L1": .01, "L2": .01, "L3": .01}, "P3": {"L1": .02, "L2": .1, "L3": .02}}
