import pandas as pd

acc_ids = SMILE.create_protein_list()
csvs = SMILE.create_ligand_list()

fas_df = pd.read_csv('fasta_list.csv', index_col='accession number')

logFC_byID = {}
pVal_byID = {}
cit_logFC = {}
cit_pVal = {}
cit_corrected = {}

for id in acc_ids:
           logFC_byID[id] = {}
           pVal_byID[id] = {}
           
def cit_labels():
           cit_df = pd.read_csv('olfr_de_copy1/olfr_de/pS6_DE_1p_citronellol.csv', index_col='name')
           
           for id in acc_ids:
                      name = fas_df.loc[id]['receptor']
                      cit_logFC[id] = (cit_df.loc[name]['logFC'])
                      cit_pVal[id] = (cit_df.loc[name]['PValue'])
                      cit_corrected[id] = (cit_df.loc[name]['logFC']*cit_df.loc[name]['No1']*cit_df.loc[name]['No2']*cit_df.loc[name]['No3']/3)
           return cit_logFC, cit_pVal, cit_corrected

def labels():
           for csv in csvs:
                      file_name = 'olfr_de_copy1/olfr_de/'+csv
                      curr_df = pd.read_csv(file_name, index_col='name')
                      for id in acc_ids:
                                 name = fas_df.loc[id]['receptor']
                                 logFC_byID[id][csv] = (curr_df.loc[name]['logFC'])
                                 pVal_byID[id][csv] = (curr_df.loc[name]['PValue'])
           return logFC_byID, pVal_byID

def classified_logFC():
    logFC_byID, pVal_byID = labels()

    classified = {}
    for id in logFC_byID:
        classified[id] = {}
        for csv in csvs:
            if logFC_byID[id][csv] < .5:
                classified[id][csv] = 0
            else:
                classified[id][csv] = 1

    return classified

def classified_pVal():
    logFC_byID, pVal_byID = labels()

    classified = {}
    for id in pVal_byID:
        classified[id] = {}
        for csv in csvs:
            if pVal_byID[id][csv] > .05:
                classified[id][csv] = 0
            else:
                classified[id][csv] = 1

    return classified

def classified_logFC_pVal(logFC_byID, pVal_byID):
    classified = {}
    i = 0
    pos_counts_id = {}
    neg_counts_id = {}

    """
    for csv in csvs:
        pos_counts_lig[csv] = 0
        neg_counts_lig[csv] = 0
    """

    for id in pVal_byID:
        id_counts = 0
        classified[id] = {}
        for csv in csvs:
            if (logFC_byID[id][csv] >= .5) & (pVal_byID[id][csv] <= .05):
                classified[id][csv] = 1
                id_counts += 1
                #pos_counts_lig[csv] += 1
            else:
                classified[id][csv] = 0
                #neg_counts_lig[csv] += 1
            i += 1
        pos_counts_id[id] = id_counts

    return classified, pos_counts_id, neg_counts_id
