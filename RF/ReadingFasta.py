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
    unique_proteins = set()

    for i in range(4):
        for id in AA_seqvar[i]:
            for kmer in AA_feat[i]:
                if kmer not in AA_seqvar[i][id]:
                    AA_seqvar[i][id][kmer] = 0
    for i in range(4):
        for id in Di_seqvar[i]:
            for kmer in Di_feat[i]:
                if kmer not in Di_seqvar[i][id]:
                    Di_seqvar[i][id][kmer] = 0

    all_ids = []
    print(len(AA_seqvar[0]))
    for id in AA_seqvar[0]:
        all_ids.append(id)

    for k in range(len(all_ids)):
        id = all_ids[k]
        freq_str = ""
        for i in range(4):
            for kmer in AA_feat[i]:
                freq_str += str(AA_seqvar[i][id][kmer])
        for i in range(4):
            for kmer in Di_feat[i]:
                freq_str += str(Di_seqvar[i][id][kmer])
        if freq_str not in unique_seqs:
            unique_seqs.add(freq_str)
            unique_proteins.add(id)
    return unique_proteins

def make_seqvar_TMS(TM_dict, TM_num, k, seqvar, feat):
    for id in TM_dict:
        #name = id
        seq = TM_dict[id][TM_num]
        #seqvar.append(Seq(name, seq, featurize(seq, k, feat)))
        seqvar[id] = featurize(seq, k, feat)
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
def makematrix(seqvar, feat, mat, unique):
    for seq in seqvar:
        newseq = []
        id = seq.name
        if id in unique:

            for kmer in feat:
                #For kmers not found in the protein, populate the matrix with zeros
                if kmer not in seq.dictionary:
                    seq.dictionary[kmer] = 0
                #Add the frequency value of the kmer
                newseq.append(seq.dictionary.get(kmer))
            #Add a frequency array for each protein
            mat.append(np.array(newseq))
    return mat

def makematrix2(seqvar, feat, mat, unique, counts):
    for id in unique:
        newseq = []
        for kmer in feat:
            #For kmers not found in the protein, populate the matrix with zeros
            if kmer not in seqvar[id]:
                seqvar[id][kmer] = 0
            #Add the frequency value of the kmer
            newseq.append(seqvar[id].get(kmer))
        #Add a frequency array for each protein
        for i in range(counts[id]):
            mat.append(np.array(newseq))
    return mat

