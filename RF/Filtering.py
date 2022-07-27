#This script filters out trivial kmers; 'trivial' is determined by the impact a given kmer has on
#whether or not a sensor / odorant pair will be positive
#Methodology for filtering out kmers is based on 'String-Based Models for Predicting RNA-Protein Interaction',
#https://dl.acm.org/doi/10.1145/3107411.3107508

#Additional coding help from:
#https://www.geeksforgeeks.org/python-program-to-find-sum-of-elements-in-list/
#https://www.geeksforgeeks.org/python-dictionary-values/

def richness_prot_imbalance(kmers, seqvar, pos_counts, neg_counts, domain, richness_level):
    """
    This function returns a list of protein kmers from a given TM domain with high distinction ability.
    It's meant to be used with an imbalanced dataset.

    Args:
        kmers (list): list of kmers found from all proteins in domain
        seqvar (dict): dictionary mapping a protein id to a frequency dictionary
            ex: seqvar[id]: key = (str) kmer, value = (int) freq. of kmer in id
        pos_counts (dict): dictionary mapping a protein id to the # of positive protein-ligand pairs with id as the protein
            ex: If the protein id only binds with the ligands L1, L2, and L3, then pos_counts[id] = 3
        neg_counts (dict): dictionary mapping a protein id to the # of negative protein-ligand pairs with id as the protein
            ex: If the protein id doesn't bind with the ligands L1 and L4, then neg_counts[id] = 2
        domain (string): indicates which TM domain the proteins are from
        richness_level: indicator of how strict the filter should be

    Returns:
        ret (list): list of protein kmers from a given TM domain with high distinctive ability
        ret2 (list): list of protein kmers from a given TM domain with high distinctive ability,
            annotated with the domain they come from

        Distinctive ability is based on richness, a measure of a kmer's frequency in positive v. negative pairs
        Richness << 1 indicates higher frequency in negative pairs
        Richness >> 1 indicates higher frequency in positive pairs
        Kmers with of richness of ~1 (approx. equal frequency in positive and negative pairs) are removed

        The formula for richness has been adapted to account for imbalances in the dataset.
        Instead of calculating richness(kmer) as (freq. of kmer in positive pairs)
                                                ----------------------------------  ,
                                                 (freq. of kmer in negative pairs)

        richness(kmer) is equivalent to (freq. of kmer in positive pairs) / (# of all kmers in positive pairs)
                                        ----------------------------------------------------------------------
                                        (freq. of kmer in negative pairs) / (# of all kmers in negative pairs)
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

    for kmer in kmers:
        if neg_counts_by_kmer[kmer] == 0:       #kmer only occurs in positive pairs
            richness[kmer] = 10000
        else:
            richness[kmer] = pos_prop_by_kmer[kmer] / neg_prop_by_kmer[kmer]

    ret = []
    ret2 = []

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
    This function returns a list of protein kmers from a given TM domain with high distinction ability.
    It's meant to be used with an balanced dataset.

    Args:
        kmers (list): list of kmers found from all proteins in domain
        seqvar (dict): dictionary mapping a protein id to a frequency dictionary
            ex: seqvar[id]: key = (str) kmer, value = (int) freq. of kmer in id
        pos_counts (dict): dictionary mapping a protein id to the # of positive protein-ligand pairs with id as the protein
            ex: If the protein id only binds with the ligands L1, L2, and L3, then pos_counts[id] = 3
        neg_counts (dict): dictionary mapping a protein id to the # of negative protein-ligand pairs with id as the protein
            ex: If the protein id doesn't bind with the ligands L1 and L4, then neg_counts[id] = 2
        domain (string): indicates which TM domain the proteins are from
        richness_level: indicator of how strict the filter should be

    Returns:
        ret (list): list of protein kmers from a given TM domain with high distinctive ability
        ret2 (list): list of protein kmers from a given TM domain with high distinctive ability, annotated
            with the domain they come from

                Distinctive ability is based on richness, a measure of a kmer's frequency in positive v. negative pairs
        Richness << 1 indicates higher frequency in negative pairs
        Richness >> 1 indicates higher frequency in positive pairs
        Kmers with of richness of ~1 (approx. equal frequency in positive and negative pairs) are removed

        Richness(kmer) is calculated as (freq. of kmer in positive pairs)
                                        ----------------------------------  ,
                                        (freq. of kmer in negative paiirs)
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

    for kmer in kmers:
        if neg_counts_by_kmer[kmer] == 0:       #kmer only occurs in positive pairs
            richness[kmer] = 10000
        else:
            richness[kmer] = pos_counts_by_kmer[kmer] / neg_counts_by_kmer[kmer]

    ret = []
    ret2 = []

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
    This function returns a list of ligand kmers with high distinction ability.
    It's meant to be used with an imbalanced dataset.

    Args:
        ligand_counts (dict): dictionary mapping a (str) ligand to a frequency dictionary
            ex: ligand_counts[lig][kmer] = (int) freq. of kmer in lig
        pos_by_lig (dict): dictionary mapping a ligand 'lig' to the # of positive
            protein-ligand pairs with lig as the ligand
            ex: if a ligand 'lig' binds with the proteins P1 and P2, then pos_by_lig[lig] = 2
        neg_by_lig (dict): dictionary mapping a ligand 'lig' to the # of negative
            protein-ligand pairs with lig as the ligand
            ex: if a ligand 'lig' does not bind with the protein P1, then neg_by_lig[lig] = 1
        richness_level: indicator of how strict the filter should be
        kmers: list of kmers found from every ligand in ligand_counts

    Returns:
        ligand_counts (dict): updated version of the parameter ligand_counts; keys are kmers from kmers_success
        kmers_success (list): list of ligand kmers with high distinctive ability

        Distinctive ability is based on richness, a measure of a kmer's frequency in positive v. negative pairs
        Richness << 1 indicates higher frequency in negative pairs
        Richness >> 1 indicates higher frequency in positive pairs
        Kmers with of richness of ~1 (approx. equal frequency in positive and negative pairs) are removed

        The formula for richness has been adapted to account for imbalances in the dataset.
        Instead of calculating richness(kmer) as (freq. of kmer in positive pairs)
                                                ----------------------------------  ,
                                                 (freq. of kmer in negative pairs)

        richness(kmer) is equivalent to (freq. of kmer in positive pairs) / (# of all kmers in positive pairs)
                                        ----------------------------------------------------------------------
                                        (freq. of kmer in negative pairs) / (# of all kmers in negative pairs)
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
    This function returns a list of ligand kmers with high distinction ability.
    It's meant to be used with an balanced dataset.

    Args:
        ligand_counts (dict): dictionary mapping a (str) ligand to a frequency dictionary
            ex: ligand_counts[lig][kmer] = (int) freq. of kmer in lig
        pos_by_lig (dict): dictionary mapping a ligand 'lig' to the # of positive
            protein-ligand pairs with lig as the ligand
            ex: if a ligand 'lig' binds with the proteins P1 and P2, then pos_by_lig[lig] = 2
        neg_by_lig (dict): dictionary mapping a ligand 'lig' to the # of negative
            protein-ligand pairs with lig as the ligand
            ex: if a ligand 'lig' does not bind with the protein P1, then neg_by_lig[lig] = 1
        richness_level: indicator of how strict the filter should be
        kmers: list of kmers found from every ligand in ligand_counts

    Returns:
        ligand_counts (dict): updated version of the parameter ligand_counts; keys are kmers from kmers_success
        kmers_success (list): list of ligand kmers that meet the filtering requirements

                Distinctive ability is based on richness, a measure of a kmer's frequency in positive v. negative pairs
        Richness << 1 indicates higher frequency in negative pairs
        Richness >> 1 indicates higher frequency in positive pairs
        Kmers with of richness of ~1 (approx. equal frequency in positive and negative pairs) are removed

        Richness(kmer) is calculated as (freq. of kmer in positive pairs)
                                        ----------------------------------  ,
                                        (freq. of kmer in negative pairs)
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