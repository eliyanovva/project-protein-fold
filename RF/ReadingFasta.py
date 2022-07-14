#This script has functions to extract k-mer frequency, 
# create objects for each protein storing accession, sequence, and k-mer frequency, 
# and create a matrix of protein k-mer frequency.

#Imports
import numpy as np

#Creating protein sequence class
class Seq:
    def __init__(self, name, sequence, dictionary):
        self.name = name
        self.sequence = sequence
        self.dictionary = dictionary #key: k-mer; value: frequency
    def __repr__(self):
        return self.name
    def __eq__(self, other):
        return self.name == other.name

#function to extract kmers and kmer-frequency from protein sequences
def featurize(seq,k,feat):
    dict = {}
    for i in range(0, len(seq) - k + 1):
        kmer = ""
        for j in range(k):
            kmer += seq[i + j]

        if kmer not in dict:
            dict[kmer] = 0
            feat.add(kmer)
        dict[kmer] += 1
    return dict

def make_seqvar_TMS(TM_dict, TM_num, k, seqvar, feat):
    for id in TM_dict:
        seq = TM_dict[id][TM_num]
        seqvar[id] = featurize(seq, k, feat)

    return seqvar, feat

#Return: List of sequence objects representing each protein; List of k-mers found in the protein
#fasta = file to read
#seqvar = list to be populated with sequence objects
#feat = list to be populated with k-mers

#Return: 2 dimensional matrix (accession number by k-mer) of frequency values
#seqvar = list of sequence objects
#feat = list of k-mers
#mat = empty matrix to be populated with frequency values

def makematrix(seqvar, feat, mat, unique_l, counts):
    for id in counts:
        newseq = []
        for kmer in feat:
                #For kmers not found in the protein, populate the matrix with zeros
            if kmer not in seqvar[id]:
                seqvar[id][kmer] = 0
                #Add the frequency value of the kmer
            newseq.append(seqvar[id].get(kmer))
            #Add a frequency array for each protein
        for lig in counts[id]:
            if unique_l.count(lig) != 0:
                mat.append(np.array(newseq))

    return mat

