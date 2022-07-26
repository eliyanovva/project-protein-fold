#This script creates the protein matrix and ligand matrix to train and test the Random Forest Algorithm

#Imports
import RF.SmileKmer as SmileKmer
import numpy as np
import RF.ReadingFasta as ReadingFasta
import RF.Globals as Globals
import RF.labels as labels
import RF.Filtering as Filtering
import RF.Duplicates as Duplicates

#Additional coding help from:
#https://numpy.org/doc/stable/reference/generated/numpy.repeat.html
#https://www.w3schools.com/python/ref_list_sort.asp
#https://www.geeksforgeeks.org/python-convert-set-into-a-list/
#https://www.geeksforgeeks.org/python-dictionary-values/
#https://numpy.org/doc/stable/reference/generated/numpy.concatenate.html

"""TM_location = "../data_files/TMdomains/TM.csv"
Di_location = "../data_files/3DiSequences/fullset_ss.fasta"
smile_location = "../Ligands_withSMILE/ligand_SMILES.csv"""""

def develop_matrices(smile_location, TM_location, Di_location):
    """
    This function creates a feature matrix of kmer frequencies for the experimental protein-ligand pairs,
    and their corresponding vector of classification labels.

    Args:
        smile_location (str): file path to data table of ligands and their SMILE formulas
        TM_location (str): file path to data table of proteins and their sequences for TMs 3,5,6 and 7
        Di_location (str): file path to fasta of proteins and their 3di sequences

    Returns:

    """
    #filter_strength = variable to set strength of the kmer filter
    #use 'None' to select no filter
    #use 'All' to only select kmers that occur exclusively in positive or negative pairs
    #otherwise, input an integer value to select the filter strength

    #prot_filter_strength = 'None'
    #prot_filter_strength = 8
    prot_filter_strength = 'All'

    #lig_filter_strength = 'None'
    #lig_filter_strength = 8
    lig_filter_strength = 'All'

    #Create classification dictionary
    acc_ids = Globals.initialize_protein_list(TM_location)
    logFC, FDR = labels.labels('../olfr_de/')
    classified, pos_counts, neg_counts, pos_dict, neg_dict, proteins_toconsider = labels.classified_logFC_FDR(logFC, FDR, acc_ids)

    total_pos = 0
    total_neg = 0
    for id in pos_counts:
        total_pos += pos_counts[id]
    for id in neg_counts:
        total_neg += neg_counts[id]

    #boolean value that describes whether or not the input dataset is balanced;
    #in other words, if there are far more / fewer positive observations than there are negative
    #BALANCED = True => balanced dataset
    #BALANCED = False => imbalanced dataset
    #affects the training algorithm and filtering protocol used
    """
    if (total_pos / total_neg > 1.5) | (total_pos / total_neg > 1.5):
        BALANCED = False
    else:
        BALANCED = True
    """

    BALANCED = False        #hardcoded for now; can uncomment the lines above

    pairs_by_prot = {}      #key = protein id, value = # of pairs involving the protein

    for id in proteins_toconsider:
        pairs_by_prot[id] = 0
    for id in pos_counts:
        pairs_by_prot[id] += pos_counts[id]
    for id in neg_counts:
        pairs_by_prot[id] += neg_counts[id]

    #Create dict of AA sequences only with proteins from pos or neg pairs
    AA_dict = Globals.initialize_AA_dict(proteins_toconsider, TM_location)

    #Create AA output for TMs 3,5,6,7
    AA_seqvar_TM3, AA_features_TM3 = ReadingFasta.make_seqvar_TMS(AA_dict, 0, 5)
    AA_seqvar_TM5, AA_features_TM5 = ReadingFasta.make_seqvar_TMS(AA_dict, 1, 5)
    AA_seqvar_TM6, AA_features_TM6 = ReadingFasta.make_seqvar_TMS(AA_dict, 2, 5)
    AA_seqvar_TM7, AA_features_TM7 = ReadingFasta.make_seqvar_TMS(AA_dict, 3, 5)

    #Filter the AA kmers for TMs 3,5,6,7
    if BALANCED == True:
        AA_filter_TM3, feat1 = Filtering.richness_prot_balance(AA_features_TM3, AA_seqvar_TM3, pos_counts, neg_counts, "TM3", prot_filter_strength)
        AA_filter_TM5, feat2 = Filtering.richness_prot_balance(AA_features_TM5, AA_seqvar_TM5, pos_counts, neg_counts, "TM5", prot_filter_strength)
        AA_filter_TM6, feat3 = Filtering.richness_prot_balance(AA_features_TM6, AA_seqvar_TM6, pos_counts, neg_counts, "TM6", prot_filter_strength)
        AA_filter_TM7, feat4 = Filtering.richness_prot_balance(AA_features_TM7, AA_seqvar_TM7, pos_counts, neg_counts, "TM7", prot_filter_strength)
    elif BALANCED == False:
        AA_filter_TM3, feat1 = Filtering.richness_prot_imbalance(AA_features_TM3, AA_seqvar_TM3, pos_counts, neg_counts, "TM3", prot_filter_strength)
        AA_filter_TM5, feat2 = Filtering.richness_prot_imbalance(AA_features_TM5, AA_seqvar_TM5, pos_counts, neg_counts, "TM5", prot_filter_strength)
        AA_filter_TM6, feat3 = Filtering.richness_prot_imbalance(AA_features_TM6, AA_seqvar_TM6, pos_counts, neg_counts, "TM6", prot_filter_strength)
        AA_filter_TM7, feat4 = Filtering.richness_prot_imbalance(AA_features_TM7, AA_seqvar_TM7, pos_counts, neg_counts, "TM7", prot_filter_strength)

    #Create dict of 3Di sequences only with proteins from pos or neg pairs
    Di_dict = Globals.initialize_3Di_dict(proteins_toconsider, TM_location, Di_location)
    #Create 3Di output for TMs 3,5,6,7
    Di_seqvar_TM3, Di_features_TM3 = ReadingFasta.make_seqvar_TMS(Di_dict, 0, 5)
    Di_seqvar_TM5, Di_features_TM5 = ReadingFasta.make_seqvar_TMS(Di_dict, 1, 5)
    Di_seqvar_TM6, Di_features_TM6 = ReadingFasta.make_seqvar_TMS(Di_dict, 2, 5)
    Di_seqvar_TM7, Di_features_TM7 = ReadingFasta.make_seqvar_TMS(Di_dict, 3, 5)

    #Filter the 3Di kmers for TMs 3,5,6,7
    if BALANCED == True:
        Di_filter_TM3, feat5 = Filtering.richness_prot_balance(Di_features_TM3, Di_seqvar_TM3, pos_counts, neg_counts, "TM3", prot_filter_strength)
        Di_filter_TM5, feat6 = Filtering.richness_prot_balance(Di_features_TM5, Di_seqvar_TM5, pos_counts, neg_counts, "TM5", prot_filter_strength)
        Di_filter_TM6, feat7 = Filtering.richness_prot_balance(Di_features_TM6, Di_seqvar_TM6, pos_counts, neg_counts, "TM6", prot_filter_strength)
        Di_filter_TM7, feat8 = Filtering.richness_prot_balance(Di_features_TM7, Di_seqvar_TM7, pos_counts, neg_counts, "TM7", prot_filter_strength)
    elif BALANCED == False:
        Di_filter_TM3, feat5 = Filtering.richness_prot_imbalance(Di_features_TM3, Di_seqvar_TM3, pos_counts, neg_counts, "TM3", prot_filter_strength)
        Di_filter_TM5, feat6 = Filtering.richness_prot_imbalance(Di_features_TM5, Di_seqvar_TM5, pos_counts, neg_counts, "TM5", prot_filter_strength)
        Di_filter_TM6, feat7 = Filtering.richness_prot_imbalance(Di_features_TM6, Di_seqvar_TM6, pos_counts, neg_counts, "TM6", prot_filter_strength)
        Di_filter_TM7, feat8 = Filtering.richness_prot_imbalance(Di_features_TM7, Di_seqvar_TM7, pos_counts, neg_counts, "TM7", prot_filter_strength)

    AA_seqvar = [AA_seqvar_TM3, AA_seqvar_TM5, AA_seqvar_TM6, AA_seqvar_TM7]
    AA_feat = [AA_filter_TM3, AA_filter_TM5, AA_filter_TM6, AA_filter_TM7]
    Di_seqvar = [Di_seqvar_TM3, Di_seqvar_TM5, Di_seqvar_TM6, Di_seqvar_TM7]
    Di_feat = [Di_filter_TM3, Di_filter_TM5, Di_filter_TM6, Di_filter_TM7]

    #Extract proteins with unique AA and 3di kmer frequencies
    unique_proteins = Duplicates.remove_proteins(AA_seqvar, AA_feat, Di_seqvar, Di_feat, pairs_by_prot, proteins_toconsider)
    unique_proteins.sort()

    pos_keys = list(pos_dict.keys())
    neg_keys = list(neg_dict.keys())
    ligands_from_unip = set()  # set of ligands that form pos or neg pairs with the set of unique proteins

    pos_by_lig = {}  # key: ligand, value: # of pos. pairs with the ligand
    neg_by_lig = {}  # key: ligand, value: # of neg. pairs with the ligand
    total_by_lig = {}  # key: ligand, value: # of all pairs the involve the ligand

    for i in range(len(pos_keys)):
        key = pos_keys[i]
        if key not in unique_proteins:
            pos_dict.pop(key)
        else:
            for lig in pos_dict[key]:
                ligands_from_unip.add(lig)

    for i in range(len(neg_keys)):
        key = neg_keys[i]
        if key not in unique_proteins:
            neg_dict.pop(key)
        else:
            for lig in neg_dict[key]:
                ligands_from_unip.add(lig)

    ligands_from_unip = list(ligands_from_unip)
    ligands_from_unip.sort()

    #Import dictionary matching ligands to SMILES String
    ligand_dict = Globals.initialize_ligand_dict(smile_location)
    #Create ligands matrix
    ligand_features, ligand_counts = SmileKmer.ligand_kmer_count(ligand_dict, 5, ligands_from_unip)

    for lig in ligands_from_unip:
        pos_by_lig[lig] = 0
        neg_by_lig[lig] = 0
        total_by_lig[lig] = 0

    for id in pos_dict:
        for lig in pos_dict[id]:
            pos_by_lig[lig] += 1
            total_by_lig[lig] += 1

    for id in neg_dict:
        for lig in neg_dict[id]:
            neg_by_lig[lig] += 1
            total_by_lig[lig] += 1

    #Update ligand_counts to only use filtered kmers
    if BALANCED == True:
        lig_counts_filter, filter_kmers = Filtering.richness_lig_balance(ligand_counts, pos_by_lig, neg_by_lig, lig_filter_strength, ligand_features)
    if BALANCED == False:
        lig_counts_filter, filter_kmers = Filtering.richness_lig_imbalance(ligand_counts, pos_by_lig, neg_by_lig, lig_filter_strength, ligand_features)

    #Extract ligands with unique kmer frequencies
    unique_ligands = Duplicates.remove_ligands(lig_counts_filter, total_by_lig)

    pos_total = 0           #num. of positive pairs with ligands from unique_ligands
    neg_total = 0           #num. of negative pairs with ligands from unique_ligands
    lig_mat = []            #matrix to store ligand features
    for id in pos_dict:
        for lig in pos_dict[id]:
            if (unique_ligands.count(lig) != 0):
                lig_mat.append(np.array(list(lig_counts_filter[lig].values())))
                pos_total += 1

    for id in neg_dict:
        for lig in neg_dict[id]:
            if (unique_ligands.count(lig) != 0):
                lig_mat.append(np.array(list(lig_counts_filter[lig].values())))
                neg_total += 1

    pos_AA_mat_TM3 = ReadingFasta.makematrix(AA_seqvar_TM3, AA_filter_TM3, [], unique_ligands, pos_dict)
    pos_AA_mat_TM5 = ReadingFasta.makematrix(AA_seqvar_TM5, AA_filter_TM5, [], unique_ligands, pos_dict)
    pos_AA_mat_TM6 = ReadingFasta.makematrix(AA_seqvar_TM6, AA_filter_TM6, [], unique_ligands, pos_dict)
    pos_AA_mat_TM7 = ReadingFasta.makematrix(AA_seqvar_TM7, AA_filter_TM7, [], unique_ligands, pos_dict)

    AA_mat_TM3 = ReadingFasta.makematrix(AA_seqvar_TM3, AA_filter_TM3, pos_AA_mat_TM3, unique_ligands, neg_dict)
    AA_mat_TM5 = ReadingFasta.makematrix(AA_seqvar_TM5, AA_filter_TM5, pos_AA_mat_TM5, unique_ligands, neg_dict)
    AA_mat_TM6 = ReadingFasta.makematrix(AA_seqvar_TM6, AA_filter_TM6, pos_AA_mat_TM6, unique_ligands, neg_dict)
    AA_mat_TM7 = ReadingFasta.makematrix(AA_seqvar_TM7, AA_filter_TM7, pos_AA_mat_TM7, unique_ligands, neg_dict)

    AA_matrix = np.concatenate((np.array(AA_mat_TM3, dtype = np.uint8), np.array(AA_mat_TM5, dtype = np.uint8),
                                np.array(AA_mat_TM6, dtype = np.uint8), np.array(AA_mat_TM7, dtype = np.uint8)) , axis = 1)

    pos_Di_mat_TM3 = ReadingFasta.makematrix(Di_seqvar_TM3, Di_filter_TM3, [], unique_ligands, pos_dict)
    pos_Di_mat_TM5 = ReadingFasta.makematrix(Di_seqvar_TM5, Di_filter_TM5, [], unique_ligands, pos_dict)
    pos_Di_mat_TM6 = ReadingFasta.makematrix(Di_seqvar_TM6, Di_filter_TM6, [], unique_ligands, pos_dict)
    pos_Di_mat_TM7 = ReadingFasta.makematrix(Di_seqvar_TM7, Di_filter_TM7, [], unique_ligands, pos_dict)

    Di_mat_TM3 = ReadingFasta.makematrix(Di_seqvar_TM3, Di_filter_TM3, pos_Di_mat_TM3, unique_ligands, neg_dict)
    Di_mat_TM5 = ReadingFasta.makematrix(Di_seqvar_TM5, Di_filter_TM5, pos_Di_mat_TM5, unique_ligands, neg_dict)
    Di_mat_TM6 = ReadingFasta.makematrix(Di_seqvar_TM6, Di_filter_TM6, pos_Di_mat_TM6, unique_ligands, neg_dict)
    Di_mat_TM7 = ReadingFasta.makematrix(Di_seqvar_TM7, Di_filter_TM7, pos_Di_mat_TM7, unique_ligands, neg_dict)

    Di_matrix = np.concatenate((np.array(Di_mat_TM3, dtype = np.uint8), np.array(Di_mat_TM5, dtype = np.uint8),
                                np.array(Di_mat_TM6, dtype = np.uint8), np.array(Di_mat_TM7, dtype = np.uint8)) , axis = 1)

    #Concatenate AA and 3Di matrices
    intermed_matrix = np.concatenate((np.array(AA_matrix, dtype = np.uint8), np.array(Di_matrix, dtype = np.uint8)) , axis = 1)

    #Concatenate protein and ligand matrices
    final_matrix = np.concatenate((intermed_matrix, np.array(lig_mat, dtype = np.uint8)), axis = 1)

    #Create Classification Vector
    pos_array = np.repeat(1, int(pos_total))
    neg_array = np.repeat(0, int(neg_total))
    logFCmat = np.concatenate((pos_array, neg_array), axis=0)
                                                            #All                            #None
    print('Pos Observations: ' + str(pos_total))    #FDR<.1: 221, FDR<.15: 113 | FDR<.1: 545, FDR<.15: 579
    print('Neg Observations: ' + str(neg_total))    #FDR<.1: 52, FDR<.15: 140  | FDR<.1: 236, FDR<.15: 490

    feat1.extend(feat2)
    feat1.extend(feat3)
    feat1.extend(feat4)
    feat1.extend(feat5)
    feat1.extend(feat6)
    feat1.extend(feat7)
    feat1.extend(feat8)
    feat1.extend(ligand_features)

    return {'X':final_matrix, 'Y':logFCmat, 'feats':feat1, 'balance':BALANCED, 'logFC_data':logFC, 
    'FDR_data':FDR, 'AA3_kmers':AA_filter_TM3, 'AA5_kmers':AA_filter_TM5, 'AA6_kmers':AA_filter_TM6, 
    'AA7_kmers':AA_filter_TM7, 'Di3_kmers':Di_filter_TM3, 'Di5_kmers':Di_filter_TM5, 'Di6_kmers':Di_filter_TM6, 
    'Di7_kmers':Di_filter_TM7, 'kmers':filter_kmers, 'uni_prot':unique_proteins, 'uni_lig':unique_ligands, 
    'AA3_seqs':AA_seqvar_TM3, 'AA5_seqs':AA_seqvar_TM5, 'AA6_seqs':AA_seqvar_TM6, 'AA7_seqs':AA_seqvar_TM7, 
    'Di3_seqs':Di_seqvar_TM3, 'Di5_seqs':Di_seqvar_TM5, 'Di6_seqs':Di_seqvar_TM6, 'Di7_seqs':Di_seqvar_TM7,
    'lig_counts':lig_counts_filter}

#result = develop_matrices('../Ligands_withSMILE/ligand_SMILEs.csv', "../data_files/TMdomains/TM.csv",
                       # "../data_files/3DiSequences/fullset_ss.fasta")