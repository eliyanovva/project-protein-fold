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

    for i in range(0, len(seq) - k + 1):
        kmer = ""
        #form the kmer
        for j in range(k):
            kmer += seq[i + j]
        #update frequency of the kmer in dict
        if kmer not in freq_dict:
            freq_dict[kmer] = 0
            #add new kmers to feat
            feat.add(kmer)
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
    for id in counts:
        newseq = []
        for kmer in feat:
            #For kmers not found in the protein, populate the matrix with zeros
            if kmer not in seqvar[id]:
                seqvar[id][kmer] = 0
            #Add the frequency value of the kmer
            newseq.append(seqvar[id].get(kmer))

        # Add a frequency array for each protein
        for lig in counts[id]:
            # Add an array for each unique ligand the protein can pair with
            if unique_l.count(lig) != 0:
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

def make_unfiltered_matrix(seqvar, feat, num_ligands):
    mat = []
    for id in list(seqvar.keys()):
        newseq = []
        for kmer in feat:
            if kmer not in seqvar[id]:
                seqvar[id][kmer]= 0
            newseq.append(seqvar[id].get(kmer))

        for i in range(num_ligands):
            mat.append(np.array(newseq))
    
    return mat