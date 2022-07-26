#Script creates new pairings out of the unused pairings (logFC > 0.1) to test the model on

import numpy as np
import CombineLigandsProteins
import labels
import ReadingFasta


#Additional coding help from:
#https://www.w3schools.com/python/ref_list_sort.asp
#https://www.geeksforgeeks.org/python-convert-set-into-a-list/
#https://www.geeksforgeeks.org/python-dictionary-values/
#https://numpy.org/doc/stable/reference/generated/numpy.concatenate.html

def develop_new_combo_matrices(smile_location, TM_location, Di_location):

    CLP_vars = CombineLigandsProteins.develop_matrices(smile_location, TM_location, Di_location)

    logFC = CLP_vars['logFC_data']
    FDR = CLP_vars['FDR_data']
    AA_seqvar_TM3 = CLP_vars['AA3_seqs']
    AA_seqvar_TM5 = CLP_vars['AA5_seqs']
    AA_seqvar_TM6 = CLP_vars['AA6_seqs']
    AA_seqvar_TM7 = CLP_vars['AA7_seqs']
    Di_seqvar_TM3 = CLP_vars['Di3_seqs']
    Di_seqvar_TM5 = CLP_vars['Di5_seqs']
    Di_seqvar_TM6 = CLP_vars['Di6_seqs']
    Di_seqvar_TM7 = CLP_vars['Di7_seqs']
    AA_filter_TM3 = CLP_vars['AA3_kmers']
    AA_filter_TM5 = CLP_vars['AA5_kmers']
    AA_filter_TM6 = CLP_vars['AA6_kmers']
    AA_filter_TM7 = CLP_vars['AA7_kmers']
    Di_filter_TM3 = CLP_vars['Di3_kmers']
    Di_filter_TM5 = CLP_vars['Di5_kmers']
    Di_filter_TM6 = CLP_vars['Di6_kmers']
    Di_filter_TM7 = CLP_vars['Di7_kmers']
    unique_proteins = CLP_vars['uni_prot']
    unique_ligands = CLP_vars['uni_lig']
    lig_counts_filter = CLP_vars['lig_counts']

    #https://www.geeksforgeeks.org/randomly-select-n-elements-from-list-in-python/

    new_combos = labels.extract_new_combos(FDR, unique_proteins, unique_ligands)

    AA_mat_TM3 = ReadingFasta.make_combomatrix(AA_seqvar_TM3, AA_filter_TM3, [], new_combos)
    AA_mat_TM5 = ReadingFasta.make_combomatrix(AA_seqvar_TM5, AA_filter_TM5, [], new_combos)
    AA_mat_TM6 = ReadingFasta.make_combomatrix(AA_seqvar_TM6, AA_filter_TM6, [], new_combos)
    AA_mat_TM7 = ReadingFasta.make_combomatrix(AA_seqvar_TM7, AA_filter_TM7, [], new_combos)

    Di_mat_TM3 = ReadingFasta.make_combomatrix(Di_seqvar_TM3, Di_filter_TM3, [], new_combos)
    Di_mat_TM5 = ReadingFasta.make_combomatrix(Di_seqvar_TM5, Di_filter_TM5, [], new_combos)
    Di_mat_TM6 = ReadingFasta.make_combomatrix(Di_seqvar_TM6, Di_filter_TM6, [], new_combos)
    Di_mat_TM7 = ReadingFasta.make_combomatrix(Di_seqvar_TM7, Di_filter_TM7, [], new_combos)

    AA_matrix = np.concatenate((np.array(AA_mat_TM3, dtype = np.uint8), np.array(AA_mat_TM5, dtype = np.uint8),
                                np.array(AA_mat_TM6, dtype = np.uint8), np.array(AA_mat_TM7, dtype = np.uint8)) , axis = 1)

    Di_matrix = np.concatenate((np.array(Di_mat_TM3, dtype = np.uint8), np.array(Di_mat_TM5, dtype = np.uint8),
                                np.array(Di_mat_TM6, dtype = np.uint8), np.array(Di_mat_TM7, dtype = np.uint8)) , axis = 1)

    lig_mat = []            #matrix to store ligand features
    for id in new_combos:
        for lig in new_combos[id]:
                lig_mat.append(np.array(list(lig_counts_filter[lig].values())))

    intermed_matrix = np.concatenate((np.array(AA_matrix, dtype = np.uint8), np.array(Di_matrix, dtype = np.uint8)) , axis = 1)
    final_matrix = np.concatenate((intermed_matrix, np.array(lig_mat, dtype = np.uint8)), axis = 1)

    return {'X': final_matrix, 'combos': new_combos}

develop_new_combo_matrices('../Ligands_withSMILE/ligand_SMILEs.csv', "../data_files/TMdomains/TM.csv",
                            "../data_files/3DiSequences/fullset_ss.fasta")
