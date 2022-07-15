#This script filters out trivial kmers; 'trivial' is determined by the impact a given kmer has on
#whether or not a sensor / odorant pair will be positive
#Methodology for filtering out kmers is based on 'String-Based Models for Predicting RNA-Protein Interaction',
#https://dl.acm.org/doi/10.1145/3107411.3107508

#richness_protein: returns a list of protein kmers that meet the filtering requirements

#kmers: set of all possible kmers for the protein
#seqvar: key = protein id, value = dict (key: kmer, value: freq. of kmer in protein)
#pos_counts: key = protein id, value = # of positive (pos) pairs with the protein
#neg_counts: key = protein id, value = # of negative (neg) pairs with the protein
def richness_protein(kmers, seqvar, pos_counts, neg_counts, domain):

    pos_counts_by_kmer = {}             #key: kmer, value: freq. of kmer in pos. pairs
    neg_counts_by_kmer = {}             #key: kmer, value: freq. of kmer in neg. pairs
    pos_prop_by_kmer = {}               #key: kmer, value: proportion of times kmer is in pos. pairs
    neg_prop_by_kmer = {}               #key: kmer, value: proportion of times kmer is in neg. pairs

    for kmer in kmers:
        pos_counts_by_kmer[kmer] = 0
        neg_counts_by_kmer[kmer] = 0

    total_pos = 0                       #total num. of kmers found in pos. pairs
    total_neg = 0                       #total num. of kmers found in neg. pairs

    for id in seqvar:                   #id = accession id of protein
        freq_dict = seqvar[id]          #freq_dict = freq. counts of all known kmers in the protein

        #increase counts of total (pos or neg) kmers by:
        #(num. of (pos or neg) pairs that involve the protein) x (num. of kmers that the protein has)
        total_pos += pos_counts[id] * sum(freq_dict.values())
        total_neg += neg_counts[id] * sum(freq_dict.values())
        #increase freq. of kmer in (pos or neg) pairs by:
        #(num. of (pos or neg) pairs that involve the protein) x (freq. of a given kmer in the protein)
        for kmer in freq_dict:
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
    # Will remove kmers with a richness ~1 (occur approx. equally in pos. and neg. pairs)
    #the formula for richness has been adapted to account for imbalances in the dataset
    #(ie, basing it off the proportion of a kmer in pos or neg pairs, instead of the raw frequency counts for the kmer)

    for kmer in kmers:
        if neg_counts_by_kmer[kmer] == 0:       #kmer only occurs in positive pairs
            richness[kmer] = 10000
        else:
            richness[kmer] = pos_prop_by_kmer[kmer] / neg_prop_by_kmer[kmer]

    richness_level = 10      #setting for how strict the richness filter is; increase the number to increase strictness
    ret = []                #list of kmers that meet filtering conditions; to be used in final matrix
    ret2 = []               #for importances list

    for kmer in kmers:
        if (richness[kmer] <= (1/richness_level)) | (richness[kmer] >= richness_level):
        #if (richness[kmer] == 10000) | (richness[kmer] == 0):
            ret.append(kmer)
            ret2.append(kmer + domain)

    return ret, ret2

#richness_ligand: returns an updated version of ligand_counts; the dictionary will now only use kmers
#that meet the filtering requirements

#ligand_counts: key: ligand, value: dict (key: kmer, value: freq. of kmer in the ligand)
#pos_by_lig: key = ligand, value = # of pos. pairs with the ligand
#neg_by_lig: key = ligand, value = # of neg. pairs with the ligand
def richness_ligand(ligand_counts, pos_by_lig, neg_by_lig):
    #kmers = list of all potential kmers for ligand
    kmers = list(ligand_counts['pS6_DE_1p_dimethyltrisulfide.csv'].keys())

    pos_counts_by_kmer = {}             #key: kmer, value: freq. of kmer in pos. pairs
    neg_counts_by_kmer = {}             #key: kmer, value: freq. of kmer in neg. pairs
    pos_prop_by_kmer = {}               #key: kmer, value: proportion of times kmer is in pos. pairs
    neg_prop_by_kmer = {}               #key: kmer, value: proportion of times kmer is in neg. pairs

    for kmer in kmers:
        pos_counts_by_kmer[kmer] = 0
        neg_counts_by_kmer[kmer] = 0

    total_pos = 0                       #total num. of kmers found in positive pairs
    total_neg = 0                       #total num. of kmers found in negative pairs

    for lig in ligand_counts:               #id = accession id of protein
        freq_dict = ligand_counts[lig]      #freq_dict = freq. counts of all known kmers in the protein

        #increase counts of total (pos or neg) kmers by:
        #(num. of (pos or neg) pairs that involve the ligand) x (num. of kmers for the ligand)
        total_pos += pos_by_lig[lig] * sum(freq_dict.values())
        total_neg += neg_by_lig[lig] * sum(freq_dict.values())
        #increase freq of kmer in (pos or neg) pairs by:
        #(num. of (pos or neg) pairs that involve the ligand) x (freq. of a given kmer in the ligand)
        for kmer in freq_dict:
            pos_counts_by_kmer[kmer] += pos_by_lig[lig] * freq_dict[kmer]
            neg_counts_by_kmer[kmer] += neg_by_lig[lig] * freq_dict[kmer]

    for kmer in kmers:
        pos_prop_by_kmer[kmer] = float(pos_counts_by_kmer[kmer]) / float(total_pos)
        neg_prop_by_kmer[kmer] = float(neg_counts_by_kmer[kmer]) / float(total_neg)

    richness = {}
    #key: kmer, value: richness measure of kmer
    #richness measures whether a kmer was more frequent in pos. or neg. pairs
    #richness << 1 indicates higher frequency in neg. pairs
    #richness >> 1 indicates higher frequency in pos. pairs
    #Will remove kmers with a richness ~1 (occur approx. equally in pos. and neg. pairs)
    #the formula for richness has been adapted to account for imbalances in the dataset
    #(ie, basing it off the prop. of a kmer in pos / neg counts, rather than the pure frequency counts for the kmer)

    for kmer in kmers:
        if neg_counts_by_kmer[kmer] == 0:       #kmer only occurs in positive pairs
            richness[kmer] = 10000
        else:
            richness[kmer] = pos_prop_by_kmer[kmer] / neg_prop_by_kmer[kmer]

    richness_level = 10      #setting for how strict the richness filter is; increase the number to increase strictness
    kmers_failed = []       #list of kmers that don't meet the filtering requirements

    for kmer in kmers:
        if (richness[kmer] > (1/richness_level)) & (richness[kmer] < richness_level):
        #if (richness[kmer] != 10000) & (richness[kmer] != 0):
            kmers_failed.append(kmer)

    #Remove all kmers from ligand_counts that didn't meet the filtering requirements
    for lig in ligand_counts:
        for kmer in kmers_failed:
            ligand_counts[lig].pop(kmer)

    return ligand_counts