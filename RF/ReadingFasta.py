#This script has functions to extract kmers and k-mer frequencies,
# and create a matrix of protein k-mer frequencies.

#Imports
import numpy as np

#make_seqvar_TMS: extracts a kmer frequency dictionary for a single transmembrane domain (TM)
#Return:
#   seqvar: key = protein id, value = dict (key: kmer, value: freq. of kmer in protein)
#   feat: list of all kmers for a given TM

#TM_dict: key = protein id, value = list(seq for TM3, seq for TM5, seq for TM6, seq for TM7)
#TM_num: int, index to extract desired TM
#k: int, kmer length
#seqvar: initialized empty dictionary
#feat: initialized empty set
def make_seqvar_TMS(TM_dict, TM_num, k, seqvar, feat):
    for id in TM_dict:
        seq = TM_dict[id][TM_num]
        seqvar[id] = featurize(seq, k, feat)

    return seqvar, feat

#featurize: function to extract kmers and kmer frequencies from a sequence
#Return:
#   dict: key = kmer, value = freq. of kmer in sequence

#seq: string of either Amino Acid or 3di input
#k: length of kmer
#feat: set of kmers found from all seq passed to featurize
def featurize(seq,k,feat):
    dict = {}

    for i in range(0, len(seq) - k + 1):
        kmer = ""
        #form the kmer
        for j in range(k):
            kmer += seq[i + j]
        #update frequency of the kmer in dict
        if kmer not in dict:
            dict[kmer] = 0
            #add new kmers to feat
            feat.add(kmer)
        dict[kmer] += 1
    return dict

#Return: 2 dimensional matrix (accession number by k-mer) of frequency values
#seqvar: key = protein id, value = dict (key: kmer, value: freq. of kmer in protein)
#feat: list of k-mers
#mat: empty matrix to be populated with frequency values
#unique_l: unique ligands to use in the final matrix
#counts: key = protein id, value = list of ligands that the protein can pair with
#   Depending on the variable passed, counts can refer to positive or negative pairs
def makematrix(seqvar, feat, mat, unique_l, counts):
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

def make_nmatrix(seqvar, feat, mat, unique, num_ligs):
    for id in unique:
        newseq = []
        for kmer in feat:
            #For kmers not found in the protein, populate the matrix with zeros
            if kmer not in seqvar[id]:
                seqvar[id][kmer] = 0
            #Add the frequency value of the kmer
            newseq.append(seqvar[id].get(kmer))
        # Add a frequency array for each protein
        for i in range(num_ligs):
                mat.append(np.array(newseq))
    return mat

