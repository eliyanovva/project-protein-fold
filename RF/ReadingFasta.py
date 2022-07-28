#This script has functions to extract kmers and k-mer frequencies,
# and create a matrix of protein k-mer frequencies.

#Imports
import numpy as np

def make_seqvar_TMS(TM_dict, TM_num, k):
    """
    This function extracts a kmer frequency dictionary for a single transmembrane (TM) domain.

    Args:
        TM_dict (dict): a dictionary mapping a protein id to a list of input sequences for each TM
            ex: TM_dict[id] = [sequence for TM3, sequence for TM5, sequence for TM6, sequence for TM7]
        TM_num (int): index to extract desired TM
        k (int): k (int): desired kmer length

    Returns:
        seqvar (dict): dictionary mapping a protein id to a frequency dictionary
            ex: seqvar[id]: key = (str) kmer, value = (int) freq. of kmer in id
    """

    seqvar = {}
    feat = set()

    # We create a kmer frequency dictionary for all protein ids.
    # A kmer is a substring of length k.
    # Each frequency dictionary is based on the input sequences
    # for a single type of transmembrane domain (TM).

    for id in TM_dict:
        seq = TM_dict[id][TM_num]
        seqvar[id] = featurize(seq, k, feat)

    return seqvar, feat

def featurize(seq,k,feat):
    """
    This function extracts kmers and kmer frequencies from seq.

    Args:
        seq (str): sequence of amino acids or 3di input
        k (int): desired kmer length
        feat (set): set of kmers found in previous calls to featurize

    Returns:
        freq_dict (dict): dictionary mapping a (str) kmer to the (int) freq. of kmer in seq
    """

    freq_dict = {}

    # We iterate through all of the characters in seq to form substrings.
    # Each substring of length k is designated as a kmer.
    # New kmers are added to feat.
    # For each kmer k that is formed, the frequency of k in freq_dict is updated.

    for i in range(0, len(seq) - k + 1):
        kmer = ""
        for j in range(k):
            kmer += seq[i + j]
        feat.add(kmer)
        if kmer not in freq_dict:
            freq_dict[kmer] = 0
        freq_dict[kmer] += 1

    return freq_dict

def makematrix(seqvar, feat, mat, unique_l, counts):
    """
    This function creates a 2-D matrix of frequency values.

    Args:
        seqvar (dict): dictionary mapping a protein id to a frequency dictionary
            ex: seqvar[id]: key = (str) kmer, value = (int) freq. of kmer in id
        feat (list): list of kmers for all the protein ids in seqvar
        mat (list): matrix to add the frequency values to; can be empty or already have values
        unique_l (list): list of unique ligands to use in the final matrix
        counts (dict): dictionary mapping a protein id to the list of ligands that id can interact with
            counts will refer exclusively either to positive or negative pairs

    Returns:
        mat (list) = matrix of frequency values
            ex: mat[i][j] = freq. of kmer j in protein i
    """

    # We create a row of frequency values for every protein from the dict counts.
    # To do so, for each protein, we extract the protein's frequency of every kmer from feat.
    # If a feat kmer isn't found in the protein, the frequency value is initialized as 0.
    # We add the frequency row to the matrix. The row is replicated once for each unique ligand
    # that the protein interacts with.

    for id in counts:
        newseq = []
        for kmer in feat:
            if kmer not in seqvar[id]:          # Check if kmer can't be found in the sequence
                seqvar[id][kmer] = 0
            newseq.append(seqvar[id].get(kmer))
        for lig in counts[id]:
            if unique_l.count(lig) != 0:        #Check that the ligand is unique
                mat.append(np.array(newseq))

    return mat

def make_combomatrix(seqvar, feat, mat, new_combos):
    for id in new_combos:
        newseq = []
        for kmer in feat:
            #For kmers not found in the protein, populate the matrix with zeros
            if kmer not in seqvar[id]:
                seqvar[id][kmer] = 0
            #Add the frequency value of the kmer
            newseq.append(seqvar[id].get(kmer))

        # Add a frequency array for each protein
        for i in range(len(new_combos[id])):
            # Add an array for each unique ligand the protein can pair with
            mat.append(np.array(newseq))
    return mat

