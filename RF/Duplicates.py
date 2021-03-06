#This script reduces the list of proteins and ligands used to make the final matrix.
#The proteins / ligands left will be those with 'unique frequencies'.
#For example, if both proteins P1 and P2 have 'unique frequencies',
#then P1 and P2 must have different frequencies of 1 or more kmers

#Additional coding help from:
#https://www.w3schools.com/python/ref_list_sort.asp
#https://www.geeksforgeeks.org/python-convert-set-into-a-list/
#https://www.geeksforgeeks.org/python-dictionary-values/
#https://www.w3schools.com/python/ref_dictionary_keys.asp

def remove_proteins(AA_seqvar, AA_feat, Di_seqvar, Di_feat, pairs_by_prot, all_ids):
    """
    This function returns a list of proteins with unique kmer frequencies.

    Args:
        AA_seqvar (list): list of Amino Acid kmer frequency dictionaries for TMs 3, 5, 6, and 7
            TM[0] = TM3, TM[1] = TM5, TM[2] = TM6, TM[3] = TM7
            ex: AA_seqvar[i][id][kmer] = (int) freq. of kmer in TM[i] of protein id
        AA_feat (list): list of Amino Acid kmer lists for TMs 3, 5, 6, and 7
            ex: AA_feat[i] (list): list of AA kmers in TM[i], post-filtering
        Di_seqvar (list): list of 3di kmer frequency dictionaries for TMs 3, 5, 6, and 7
            TM[0] = TM3, TM[1] = TM5, TM[2] = TM6, TM[3] = TM7
            ex: Di_seqvar[i][id][kmer] = (int) freq. of kmer in TM[i] of protein id
        Di_feat (list): list of 3di kmer lists for TMs 3, 5, 6, and 7
            ex: Di_feat[i]: list of 3di kmers in TM[i], post-filtering
        pairs_by_prot (dict): dictionary mapping a protein id to the # of protein-ligand pairs with id as the protein
            Refers to both positive and negative pairs.
            ex: If protein id binds with the ligands L1 and L2, and doesn't bind with the ligand L3, then
            pairs_by_prot[id] = 3
        all_ids (list): list of all protein ids used as keys in AA_seqvar and Di_seqvar

    Returns:
        unique_proteins (list): list of proteins that are associated with unique kmer frequencies
            Let K be a set of kmers. Proteins Pi and Pj have unique kmer frequences if, for at least
            one kmer k in K, k does not have the same frequency in Pi and Pj.
    """
    unique_seqs = {}
    # key: AA and 3di kmer frequencies (formatted as string), value: protein with those frequencies
    # will only store unique kmer frequencies as keys

    if (len(AA_seqvar) == 0) | (len(Di_seqvar) == 0):
        print('Invalid Sequence dictionaries')
        return []

    else:
        #check that AA_seqvar[i] has a frequency value for all kmers in AA_feat[i], for all i
        #ensures standardization when the freq_str are created later on
        for i in range(4):
            for id in AA_seqvar[i]:
                for kmer in AA_feat[i]:
                    #if the kmer isn't stored in AA_seqvar[i], initialize the kmer frequency as 0
                    if kmer not in AA_seqvar[i][id]:
                        AA_seqvar[i][id][kmer] = 0

        # check that Di_seqvar[i] has a frequency value for all kmers in Di_feat[i], for all i
        for i in range(4):
            for id in Di_seqvar[i]:
                for kmer in Di_feat[i]:
                    # if the kmer isn't stored in Di_seqvar[i], initialize the kmer frequency as 0
                    if kmer not in Di_seqvar[i][id]:
                        Di_seqvar[i][id][kmer] = 0

        for k in range(len(all_ids)):
            id = all_ids[k]
            #format the AA and 3di kmer frequencies of id as a single string: freq_str
            #freq_str will include the frequencies from all the TMs: TM3, TM5, TM6, and TM7
            freq_str = ""
            for i in range(4):
                for kmer in AA_feat[i]:
                    freq_str += str(AA_seqvar[i][id][kmer])
            for i in range(4):
                for kmer in Di_feat[i]:
                    freq_str += str(Di_seqvar[i][id][kmer])

            # if freq_str is unique, automatically store it in unique_seqs
            if freq_str not in unique_seqs:
                unique_seqs[freq_str] = id
            # otherwise, compare the 2 proteins with the kmer frequencies of freq_str
            # the proteins that has more pairs will be used in the final matrix
            else:
                old_id = unique_seqs[freq_str]
                if pairs_by_prot[id] > pairs_by_prot[old_id]:
                    unique_seqs[freq_str] = id

        unique_proteins = list(unique_seqs.values())

        return unique_proteins

def remove_ligands(ligand_counts, total_by_lig):
    """
    This function returns a list of ligands with unique frequencies.

    Args:
        ligand_counts (dict): dictionary mapping a (str) ligand to a frequency dictionary
            ex: ligand_counts[lig][kmer] = (int) freq. of kmer in lig
        total_by_lig (dict): dictionary mapping a ligand 'lig' to the # of protein-ligand pairs with lig as the ligand
            Refers to both positive and negative pairs.
            ex: If a ligand 'lig' binds with protein P1 and does not bind with proteins P2, P3, and P4,
            then total_by_lig[lig] = 4

    Returns:
        unique_ligands: list of ligands that are associated with unique kmer frequencies
            Let K be a set of kmers. Ligands Li and Lj have unique kmer frequences if, for at least
            one kmer k in K, k does not have the same frequency in Li and Lj.
    """
    unique_seqs = {}
    #key: kmer frequencies (formatted as string), value: ligand with those frequencies
    #will only store unique kmer frequencies as keys

    ligands = list(ligand_counts.keys())
    ligands.sort()

    for lig in ligands:
        #format the kmer frequencies of lig as a string: freq_str
        freq_str = ""
        for kmer in ligand_counts[lig]:
            freq_str += str(ligand_counts[lig][kmer])

        #if freq_str is unique, automatically store it in unique_seqs
        if freq_str not in unique_seqs:
            unique_seqs[freq_str] = lig
        #otherwise, compare the 2 ligands with the kmer frequencies of freq_str
        #the ligand that has more pairs will be used in the final matrix
        else:
            old_lig = unique_seqs[freq_str]
            if total_by_lig[lig] > total_by_lig[old_lig]:
                unique_seqs[freq_str] = lig

    unique_ligands = list(unique_seqs.values())

    return unique_ligands