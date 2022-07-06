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

def remove_duplicates(AA_seqvar, AA_feat, Di_seqvar, Di_feat):
    unique_seqs = set()
    unique_proteins = []

    print(AA_seqvar[0])
    print(Di_seqvar[0])
    print(AA_feat[0])
    print(Di_feat[0])

    for i in range(4):
        for seq in AA_seqvar[i]:
            for kmer in AA_feat[i]:
                if kmer not in seq.dictionary:
                    seq.dictionary[kmer] = 0
    for i in range(4):
        for seq in Di_seqvar[i]:
            for kmer in Di_feat[i]:
                if kmer not in seq.dictionary:
                    seq.dictionary[kmer] = 0

    all_ids = []
    for seq in AA_seqvar[0]:
        all_ids.append(seq.name)

    for id in all_ids:
        freq_str = ""
        for i in range(4):
            for seq in AA_seqvar[i]:
                if seq.name == id:
                    for kmer in AA_feat[i]:
                        freq_str += str(seq.dictionary[kmer])
        for i in range(4):
            for seq in Di_seqvar[i]:
                if seq.name == id:
                    for kmer in Di_feat[i]:
                        freq_str += str(seq.dictionary[kmer])
        print(len(freq_str))
        if freq_str not in unique_seqs:
            unique_seqs.add(freq_str)
            unique_proteins.append(id)

    """
    for seq in seqvar_a3:
        for kmer in feat_a3:
            if kmer not in seq.dictionary:
                seq.dictionary[kmer] = 0

    for seq in seqvar_a5:
        for kmer in feat_a5:
            if kmer not in seq.dictionary:
                seq.dictionary[kmer] = 0

    for seq in seqvar_a6:
        for kmer in feat_a6:
            if kmer not in seq.dictionary:
                seq.dictionary[kmer] = 0

    for seq in seqvar_a7:
        for kmer in feat_a7:
            if kmer not in seq.dictionary:
                seq.dictionary[kmer] = 0

    for seq in seqvar_d3:
        for kmer in feat_d3:
            if kmer not in seq.dictionary:
                seq.dictionary[kmer] = 0

    for seq in seqvar_d5:
        for kmer in feat_d5:
            if kmer not in seq.dictionary:
                seq.dictionary[kmer] = 0

    for seq in seqvar_d6:
        for kmer in feat_d6:
            if kmer not in seq.dictionary:
                seq.dictionary[kmer] = 0

    for seq in seqvar_d7:
        for kmer in feat_d7:
            if kmer not in seq.dictionary:
                seq.dictionary[kmer] = 0

    all_ids = []
    for seq in seqvar_a3:
        all_ids.append(seq.name)

    for id in all_ids:
        freq_str = ""
        for seq in seqvar_a3:
            if seq.name == id:
                for kmer in feat_a3:
                    freq_str += str(seq.dictionary[kmer])
        for seq in seqvar_a5:
            if seq.name == id:
                for kmer in feat_a5:
                    freq_str += str(seq.dictionary[kmer])
        for seq in seqvar_a6:
            if seq.name == id:
                for kmer in feat_a6:
                    freq_str += str(seq.dictionary[kmer])
        for seq in seqvar_a7:
            if seq.name == id:
                for kmer in feat_a7:
                    freq_str += str(seq.dictionary[kmer])
        for seq in seqvar_d3:
            if seq.name == id:
                for kmer in feat_d3:
                    freq_str += str(seq.dictionary[kmer])
        for seq in seqvar_d5:
            if seq.name == id:
                for kmer in feat_d5:
                    freq_str += str(seq.dictionary[kmer])
        for seq in seqvar_d6:
            if seq.name == id:
                for kmer in feat_d6:
                    freq_str += str(seq.dictionary[kmer])
        for seq in seqvar_d7:
            if seq.name == id:
                for kmer in feat_d7:
                    freq_str += str(seq.dictionary[kmer])
        """

    print(len(unique_proteins))


def make_seqvar_TMS(TM_dict, TM_num, k, seqvar, feat):
    for id in TM_dict:
        name = id
        seq = TM_dict[name][TM_num]
        seqvar.append(Seq(name, seq, featurize(seq, k, feat)))
    return seqvar, feat

#Return: List of sequence objects representing each protein; List of k-mers found in the protein
#fasta = file to read
#seqvar = list to be populated with sequence objects
#feat = list to be populated with k-mers
def make_seqvar(fasta, seqvar, feat):
    i = 0
    j = 0
    for line in fasta:
        #Identify and add new sequence objects
        if line[0] == '>':
            name = line.replace('\n','')
            name = name.replace('>', '')
            seqvar.append(Seq(name, '', featurize('',7,feat)))
            i += 1
        #Update sequence objects to store sequences and frequency dictionaries
        else:
            sequence = line.replace('\n','')
            seqvar[i-1].sequence = seqvar[i-1].sequence + sequence
            seqvar[i-1].dictionary = featurize(seqvar[i-1].sequence,7,feat)
        j += 1

    return seqvar, feat

#Return: 2 dimensional matrix (accession number by k-mer) of frequency values
#seqvar = list of sequence objects
#feat = list of k-mers
#mat = empty matrix to be populated with frequency values
def makematrix(seqvar, feat, mat, all_freqs):
    for seq in seqvar:
        newseq = []

        id = seq.name
        if id not in all_freqs:
            all_freqs[id] = []

        for kmer in feat:
            #For kmers not found in the protein, populate the matrix with zeros
            if kmer not in seq.dictionary:
                seq.dictionary[kmer] = 0
            #Add the frequency value of the kmer
            newseq.append(seq.dictionary.get(kmer))
        #Add a frequency array for each protein
        mat.append(np.array(newseq))
        all_freqs[id].append(newseq)
    return mat
