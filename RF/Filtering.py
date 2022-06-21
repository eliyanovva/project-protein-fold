import SMILE

acc_list = SMILE.create_protein_list()

###WARNING: DONT USE THE RICHNESS_LIGAND FUNCTION YET!!!!!!
def richness_ligand(kmers, freq_dict, pos_counts, neg_counts):
    pos_counts_by_kmer = {}
    neg_counts_by_kmer = {}
    for kmer in kmers:
        pos_counts_by_kmer[kmer] = 0
        neg_counts_by_kmer[kmer] = 0

    total_pos = sum(list(pos_counts.values()))
    total_neg = sum(list(neg_counts.values()))

    for kmer in kmers:
        for lig in freq_dict:
            if (pos_counts[lig] > 0) & (freq_dict[lig][kmer] > 0):
                pos_counts_by_kmer[kmer] += pos_counts[lig]
                neg_counts_by_kmer[kmer] += neg_counts[lig]

    for kmer in kmers:
        orig_pos = pos_counts_by_kmer[kmer]
        orig_neg = neg_counts_by_kmer[kmer]
        pos_counts_by_kmer[kmer] = float(orig_pos) / float(total_pos)
        neg_counts_by_kmer[kmer] = float(orig_neg) / float(total_neg)

    richness = {}
    for kmer in kmers:
        if neg_counts_by_kmer[kmer] == 0:
            richness[kmer] = 10
        else:
            richness[kmer] = pos_counts_by_kmer[kmer] / neg_counts_by_kmer[kmer]

    kmer_to_remove = []
    for kmer in kmers:
        if (richness[kmer] > .8) & (richness[kmer] < 1.2):
            kmer_to_remove.append(kmer)

    for lig in freq_dict:
        for kmer in kmer_to_remove:
            freq_dict[lig].pop(kmer)

    return freq_dict

def richness_protein(kmers, pos_counts, neg_counts):
    pos_counts_by_kmer = {}
    neg_counts_by_kmer = {}
    for kmer in kmers:
        pos_counts_by_kmer[kmer] = 0
        neg_counts_by_kmer[kmer] = 0

    total_pos = sum(list(pos_counts.values()))
    total_neg = sum(list(neg_counts.values()))

    for kmer in kmers:
        for lig in freq_dict:
            if (pos_counts[lig] > 0) & (freq_dict[lig][kmer] > 0):
                pos_counts_by_kmer[kmer] += pos_counts[lig]
                neg_counts_by_kmer[kmer] += neg_counts[lig]

    for kmer in kmers:
        orig_pos = pos_counts_by_kmer[kmer]
        orig_neg = neg_counts_by_kmer[kmer]
        pos_counts_by_kmer[kmer] = float(orig_pos) / float(total_pos)
        neg_counts_by_kmer[kmer] = float(orig_neg) / float(total_neg)

    richness = {}
    for kmer in kmers:
        if neg_counts_by_kmer[kmer] == 0:
            richness[kmer] = 10
        else:
            richness[kmer] = pos_counts_by_kmer[kmer] / neg_counts_by_kmer[kmer]

    kmer_to_remove = []
    for kmer in kmers:
        if (richness[kmer] > .8) & (richness[kmer] < 1.2):
            kmer_to_remove.append(kmer)

    for lig in freq_dict:
        for kmer in kmer_to_remove:
            freq_dict[lig].pop(kmer)

    return freq_dict


"""
protein matrix:
seq1 ~ kmer1 kmer2 kmer3
seq2 ~ kmer1 kmer2 kmer3
seq3 ~ kmer1 kmer2 kmer3
"""
