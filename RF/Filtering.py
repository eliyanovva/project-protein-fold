#This script filters out trivial kmers; 'trivial' is determined by the impact a given kmer has on
#whether or not a sensor / odorant pair will be positive
#Methodology for filtering out kmers is based on 'String-Based Models for Predicting RNA-Protein Interaction',
#https://dl.acm.org/doi/10.1145/3107411.3107508

#Additional coding help from:
#https://www.geeksforgeeks.org/python-program-to-find-sum-of-elements-in-list/
#https://www.geeksforgeeks.org/python-dictionary-values/

def richness_prot_imbalance(kmers, seqvar, pos_counts, neg_counts, domain, richness_level):
    """
    This function returns a list of protein kmers from a given TM domain that meet the filtering requirements.
    It's meant to be used with an imbalanced dataset.
    Args:
        kmers (list): list of kmers found from all proteins in domain
        seqvar (dict): key = (string) id,
            value = (dict) key = (string) kmer, value = (int) freq. of kmer in id
        pos_counts (dict): key = (string) id, value = (int) # of pos protein-ligand pairs with id as the protein
        neg_counts (dict): key = (string) id, value = (int) # of neg protein-ligand pairs with id as the protein
        domain (string): indicates which TM domain the proteins are from
        richness_level (int or string): indicator of how strict the filter should be

    Returns:
        ret (list): list of protein kmers from a given TM domain that meet the filtering requirements
        ret2 (list): list of protein kmers from a given TM domain that meet the filtering requirements, annotated
            with the domain they come from
    """

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

    ret = []                #list of kmers that meet filtering conditions; to be used in final matrix
    ret2 = []               #for importances list

    if richness_level == 'All':
        for kmer in kmers:
            if (richness[kmer] == 10000) | (richness[kmer] == 0):
                ret.append(kmer)
                ret2.append(kmer + domain)
    elif richness_level == 'None':
        for kmer in kmers:
                ret.append(kmer)
                ret2.append(kmer + domain)
    else:
        for kmer in kmers:
            if (richness[kmer] <= (1/richness_level)) | (richness[kmer] >= richness_level):
                ret.append(kmer)
                ret2.append(kmer + domain)

    return ret, ret2

def richness_prot_balance(kmers, seqvar, pos_counts, neg_counts, domain, richness_level):
    """
    This function returns a list of protein kmers from a given TM domain that meet the filtering requirements.
    It's meant to be used with an balanced dataset.
    Args:
        kmers (list): list of kmers found from all proteins in domain
        seqvar (dict): key = (string) id,
            value = (dict) key = (string) kmer, value = (int) freq. of kmer in id
        pos_counts (dict): key = (string) id, value = (int) # of pos protein-ligand pairs with id as the protein
        neg_counts (dict): key = (string) id, value = (int) # of neg protein-ligand pairs with id as the protein
        domain (string): indicates which TM domain the proteins are from
        richness_level (int or string): indicator of how strict the filter should be

    Returns:
        ret (list): list of protein kmers from a given TM domain that meet the filtering requirements
        ret2 (list): list of protein kmers from a given TM domain that meet the filtering requirements, annotated
            with the domain they come from
    """

    pos_counts_by_kmer = {}             #key: kmer, value: freq. of kmer in pos. pairs
    neg_counts_by_kmer = {}             #key: kmer, value: freq. of kmer in neg. pairs

    for kmer in kmers:
        pos_counts_by_kmer[kmer] = 0
        neg_counts_by_kmer[kmer] = 0

    for id in seqvar:                   #id = accession id of protein
        freq_dict = seqvar[id]          #freq_dict = freq. counts of all known kmers in the protein

        #increase freq. of kmer in (pos or neg) pairs by:
        #(num. of (pos or neg) pairs that involve the protein) x (freq. of a given kmer in the protein)
        for kmer in freq_dict:
            pos_counts_by_kmer[kmer] += pos_counts[id] * freq_dict[kmer]
            neg_counts_by_kmer[kmer] += neg_counts[id] * freq_dict[kmer]

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
            richness[kmer] = pos_counts_by_kmer[kmer] / neg_counts_by_kmer[kmer]

    ret = []                #list of kmers that meet filtering conditions; to be used in final matrix
    ret2 = []               #for importances list

    if richness_level == 'All':
        for kmer in kmers:
            if (richness[kmer] == 10000) | (richness[kmer] == 0):
                ret.append(kmer)
                ret2.append(kmer + domain)
    elif richness_level == 'None':
        for kmer in kmers:
                ret.append(kmer)
                ret2.append(kmer + domain)
    else:
        for kmer in kmers:
            if (richness[kmer] <= (1/richness_level)) | (richness[kmer] >= richness_level):
                ret.append(kmer)
                ret2.append(kmer + domain)

    return ret, ret2

def richness_lig_imbalance(ligand_counts, pos_by_lig, neg_by_lig, richness_level, kmers):
    """
    This function returns a list of ligand kmers that meet the filtering requirements.
    It's meant to be used with an imbalanced dataset.
    Args:
        ligand_counts (dict): key = (string) ligand,
            value = (dict) key = (string) kmer, value = (int) freq. of kmer in ligand
        pos_by_lig (dict): key = (string) lig, value = (int) # of pos protein-ligand pairs with lig as the ligand
        neg_by_lig (dict): key = (string) lig, value = (int) # of neg protein-ligand pairs with lig as the ligand
        richness_level: richness_level (int or string): indicator of how strict the filter should be
        kmers: list of kmers found from every ligand in ligand_counts

    Returns:
        ligand_counts (dict): updated version of the parameter ligand_counts; keys are kmers from kmers_success
        kmers_success (list): list of ligand kmers that meet the filtering requirements
    """

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

    for kmer in kmers:
        if neg_counts_by_kmer[kmer] == 0:       #kmer only occurs in positive pairs
            richness[kmer] = 10000
        else:
            richness[kmer] = pos_prop_by_kmer[kmer] / neg_prop_by_kmer[kmer]

    kmers_success = []          #list of kmers that meet the filtering requirements
    kmers_failed = []           #list of kmers that don't meet the filtering requirements

    if richness_level == 'All':
        for kmer in kmers:
            if (richness[kmer] != 10000) & (richness[kmer] != 0):
                kmers_failed.append(kmer)
            else:
                kmers_success.append(kmer)
    elif richness_level == 'None':
        kmers_failed = []
        kmers_success = kmers
    else:
        for kmer in kmers:
            if (richness[kmer] > (1/richness_level)) & (richness[kmer] < richness_level):
                kmers_failed.append(kmer)
            else:
                kmers_success.append(kmer)

    #Remove all kmers from ligand_counts that didn't meet the filtering requirements
    for lig in ligand_counts:
        for kmer in kmers_failed:
            ligand_counts[lig].pop(kmer)

    return ligand_counts, kmers_success

#Ligand filtering for balanced dataset
def richness_lig_balance(ligand_counts, pos_by_lig, neg_by_lig, richness_level, kmers):
    """
    This function returns a list of ligand kmers that meet the filtering requirements.
    It's meant to be used with an balanced dataset.
    Args:
        ligand_counts (dict): key = (string) ligand,
            value = (dict) key = (string) kmer, value = (int) freq. of kmer in ligand
        pos_by_lig (dict): key = (string) lig, value = (int) # of pos protein-ligand pairs with lig as the ligand
        neg_by_lig (dict): key = (string) lig, value = (int) # of neg protein-ligand pairs with lig as the ligand
        richness_level: richness_level (int or string): indicator of how strict the filter should be
        kmers: list of kmers found from every ligand in ligand_counts

    Returns:
        ligand_counts (dict): updated version of the parameter ligand_counts; keys are kmers from kmers_success
        kmers_success (list): list of ligand kmers that meet the filtering requirements

    """

    pos_counts_by_kmer = {}             #key: kmer, value: freq. of kmer in pos. pairs
    neg_counts_by_kmer = {}             #key: kmer, value: freq. of kmer in neg. pairs

    for kmer in kmers:
        pos_counts_by_kmer[kmer] = 0
        neg_counts_by_kmer[kmer] = 0

    for lig in ligand_counts:               #id = accession id of protein
        freq_dict = ligand_counts[lig]      #freq_dict = freq. counts of all known kmers in the protein

        #increase freq of kmer in (pos or neg) pairs by:
        #(num. of (pos or neg) pairs that involve the ligand) x (freq. of a given kmer in the ligand)
        for kmer in freq_dict:
            pos_counts_by_kmer[kmer] += pos_by_lig[lig] * freq_dict[kmer]
            neg_counts_by_kmer[kmer] += neg_by_lig[lig] * freq_dict[kmer]

    richness = {}

    for kmer in kmers:
        if neg_counts_by_kmer[kmer] == 0:       #kmer only occurs in positive pairs
            richness[kmer] = 10000
        else:
            richness[kmer] = pos_counts_by_kmer[kmer] / neg_counts_by_kmer[kmer]

    kmers_success = []          #list of kmers that meet the filtering requirements
    kmers_failed = []           #list of kmers that don't meet the filtering requirements

    if richness_level == 'All':
        for kmer in kmers:
            if (richness[kmer] != 10000) & (richness[kmer] != 0):
                kmers_failed.append(kmer)
            else:
                kmers_success.append(kmer)
    elif richness_level == 'None':
        kmers_failed = []
        kmers_success = kmers
    else:
        for kmer in kmers:
            if (richness[kmer] > (1/richness_level)) & (richness[kmer] < richness_level):
                kmers_failed.append(kmer)
            else:
                kmers_success.append(kmer)

    for lig in ligand_counts:
        for kmer in kmers_failed:
            ligand_counts[lig].pop(kmer)

    return ligand_counts, kmers_success