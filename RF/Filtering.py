import Globals

acc_list = Globals.initialize_protein_list()

def richness_protein(kmers, seqvar, pos_counts, neg_counts):
    pos_counts_by_kmer = {}
    neg_counts_by_kmer = {}
    for kmer in kmers:
        pos_counts_by_kmer[kmer] = 0
        neg_counts_by_kmer[kmer] = 0

    total_pos = sum(list(pos_counts.values()))
    total_neg = sum(list(neg_counts.values()))

    for seq in seqvar:
        id = seq.name
        freq_dict = seq.dictionary
        for kmer in kmers:
            if (pos_counts[id] > 0) & (kmer in freq_dict):
                pos_counts_by_kmer[kmer] += pos_counts[id]
                neg_counts_by_kmer[kmer] += neg_counts[id]

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
    ret = []
    for kmer in kmers:
        if (richness[kmer] < .5) | (richness[kmer] > 1.5):
            ret.append(kmer)
    return ret
