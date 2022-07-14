#This script reduces the list of proteins and ligands used to make the final matrix.
#The proteins / ligands left will be those with 'unique frequencies'.
#For example, if both proteins P1 and P2 have 'unique frequencies',
#then P1 and P2 must have different frequencies of 1 or more kmers


#remove_protein: returns a list of proteins with unique frequencies

#AA_seqvar: list of Amino Acid (AA) kmer freq. dictionaries for TMs 3, 5, 6, and 7
#       AA_seqvar[i]: key = protein id, value = dict (key: kmer, value: freq. of AA kmer in protein for TM[i])
#AA_feat: list of lists of all AA kmers for TMs 3, 5, 6, and 7, respectively
#       AA_feat[i]: list of AA kmers in TM[i], post-filtering
#Di_seqvar: list of 3di kmer freq. dictionaries for TMs 3, 5, 6, and 7
#       Di_seqvar[i]: key = protein id, value = dict (key: kmer, value: freq. of 3di kmer in protein for TM[i])
#Di_feat: list of lists of all 3di kmers for TMs 3, 5, 6, and 7, respectively
#       Di_feat[i]: list of 3di kmers in TM[i], post-filtering
#pairs_by_prot: key = protein id, value = # of pairs involving the protein
#all_ids: list of all protein ids used as keys in AA_seqvar and Di_seqvar
def remove_proteins(AA_seqvar, AA_feat, Di_seqvar, Di_feat, pairs_by_prot, all_ids):
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
        # list of proteins that are associated with the unique kmer frequencies
        return unique_proteins

#remove_ligands: returns a list of ligands with unique frequencies

#ligand_counts: key: ligand, value: dict (key: kmer, value: freq. of kmer in the ligand)
#total_by_lig: key: ligand, value: # of pairs the involve the ligand
def remove_ligands(ligand_counts, total_by_lig):
    unique_seqs = {}
    #key: kmer frequencies (formatted as string), value: ligand with those frequencies
    #will only store unique kmer frequencies as keys

    for lig in ligand_counts:
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
    #set of proteins that are associated with the unique kmer frequencies

    return unique_ligands