#This script filters out trivial kmers; 'trivial' is determined by the impact a given kmer has on
#whether or not a sensor / odorant pair will be positive
#Methodology for filtering out kmers is based on 'String-Based Models for Predicting RNA-Protein Interaction',
#https://dl.acm.org/doi/10.1145/3107411.3107508

import Globals

#Import list of accession numbers
acc_list = Globals.initialize_protein_list()

#kmers: set of all possible kmers for the protein
#seqvar: list of Seq class objects (name, sequence = protein sequence, dictionary = kmer freq dictionary)
#pos_counts: key = protein id, value = # pos. interactions with the protein
#neg_counts: key = protein id, value = # neg. interactions with the protein
def richness_protein(kmers, seqvar, pos_counts, neg_counts):
    kmers = list(kmers)

    pos_counts_by_kmer = {}             #key: kmer, value: num. of positive pairs that involve the kmer
    neg_counts_by_kmer = {}             #key: kmer, value: num. of negative pairs that involve the kmer
    pos_prop_by_kmer = {}               #key: kmer, value: proportion of positive pairs that involve the kmer
    neg_prop_by_kmer = {}               #key: kmer, value: proportion of negative pairs that involve the kmer
    for kmer in kmers:
        pos_counts_by_kmer[kmer] = 0
        neg_counts_by_kmer[kmer] = 0

    total_pos = 0                       #total num. of kmers involved in positive pairs
    total_neg = 0                       #total num. of kmers involved in negative pairs

    for seq in seqvar:
        id = seq.name                   #id = accession id of protein
        freq_dict = seq.dictionary      #freq_dict = freq. counts of all known kmers in the protein

        #increase total kmer counts by (num. of pairs that involve the protein) x (num. of kmers that the protein has)
        total_pos += pos_counts[id] * sum(freq_dict.values())
        total_neg += neg_counts[id] * sum(freq_dict.values())

        #increase by_kmer counts by (num. of pairs that involve the protein) x (freq. of a given kmer in the protein)
        for kmer in kmers:
            if (pos_counts[id] > 0) & (kmer in freq_dict):
                pos_counts_by_kmer[kmer] += pos_counts[id] * freq_dict[kmer]
                neg_counts_by_kmer[kmer] += neg_counts[id] * freq_dict[kmer]

    for kmer in kmers:
        pos_prop_by_kmer[kmer] = float(pos_counts_by_kmer[kmer]) / float(total_pos)
        neg_prop_by_kmer[kmer] = float(neg_counts_by_kmer[kmer]) / float(total_neg)

    richness = {}
    #dict richness stores an adapted richness measure for all kmers
    #richness measures whether a kmer was more frequent in positive or negative pairs
    #richness << 1 indicates higher frequency in negative pairs
    #richness >> 1 indicates higher frequency in positive pairs
    #the formula for richness has been adapted to account for imbalances in the dataset
    #(ie, basing it off the prop. of a kmer in pos / neg counts, rather than the pure frequency counts for the kmer)
    for kmer in kmers:
        if neg_counts_by_kmer[kmer] == 0:       #kmer only occurs in positive pairs
            richness[kmer] = 100000
        else:
            richness[kmer] = pos_prop_by_kmer[kmer] / neg_prop_by_kmer[kmer]

    ret = []                #list of kmers that meet filtering conditions; to be used in final matrix
    for kmer in kmers:
        if (richness[kmer] < .125) | (richness[kmer] > 8):
            ret.append(kmer)

    return ret