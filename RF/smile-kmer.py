def form_letters(smile):
    letters = []  # list that stores the sectioned off 'letters' of the str smile
    for i in range(0, len(smile)):
        letters.append(0)
    a_index = 0  # index of the latest atom
    s_index = 0  # index of the starting point of the latest side chain
    e_index = 0  # index of the ending point of the latest side chain
    mid_index = 0  # index of the latest side chain that picks up again after a nested side chain
    sides = 0  # current number of nested side chains

    for i in range(0, len(smile)):
        if smile[i] == "(":         #a new side chain has begun
            s_index = i
            e_index = i + 2
            letters[i] = "("
            sides += 1
        elif smile[i] == ")":       #a side chain has terminated
            e_index = i
            sides -= 1
            if s_index > mid_index:     #the current side wasn't nested in other side chains
                orig = letters[s_index]
                letters[s_index] = orig + ")"
            else:
                orig = letters[mid_index]   #the current side was nested in other side chains
                letters[mid_index] = orig + ")"
        elif smile[i].isdigit():            #checks if the current character in string is part of the numerical notation
            if a_index < e_index:           #for some atom
                orig = letters[s_index]
                letters[s_index] = orig + str(smile[i])
                e_index += 1
            elif sides > 0:
                orig = letters[mid_index]
                letters[mid_index] = orig + str(smile[i])
            else:
                orig = letters[a_index]
                letters[a_index] = orig + str(smile[i])
        else:
            a_index = i
            if a_index < e_index:               #current atom is part of a non-nested side chain
                orig = letters[s_index]
                letters[s_index] = orig + smile[i]
                e_index += 1
            elif sides > 0:                     #current atom is following a nested side chain
                if smile[i-1] == ")":
                    letters[i] = (smile[i])
                    mid_index = i
                else:
                    orig = letters[mid_index]
                    letters[mid_index] = orig + smile[i]
            else:                               #current atom is part of the backbone
                letters[i] = (smile[i])

    #remove indices from letters that were not used
    c = letters.count(0)
    for i in range(0, c):
        letters.remove(0)

def smile_list(smile, k):
    kmer_list = []
    letters = form_letters(smile)

    for i in range(0, len(letters) - k + 1):
        k_ster = ""
        for j in range(k):
            k_ster += letters[i+j]
        if kmer_list.count(k_ster) == 0:
            kmer_list.append(k_ster)

    return kmer_list

def smile_dict(smile, k):
    #Option 1 ~  Characters are:
    #atoms in the background
    #the entirety of any side chains
    #double bonds
    kmer_dict = {}
    letters = form_letters(smile)

    for i in range(0, len(letters) - k + 1):
        k_ster = ""
        for j in range(k):
            k_ster += letters[i+j]
        if k_ster not in kmer_dict:
            kmer_dict[k_ster] = 0
        kmer_dict[k_ster] += 1

    return kmer_dict

def find_total_kmers(ligands, k):
    kmers = []
    for lig in ligands:
        k_list = smile_list(ligands[lig], k)
        for kmer in k_list:
            if kmers.count(kmer) == 0:
                kmers.append(kmer)
    return kmers

def ligand_kmer_count(ligands, k):
    ligand_counts = {}
    total_kmers = find_total_kmers(ligands, k)
    for lig in ligands:
        lig_dict = {}
        for kmer in total_kmers:
            lig_dict[kmer] = 0
        freq_dict = smile_dict(ligands[lig], k)
        for kster in freq_dict:
            lig_dict[kster] = freq_dict[kster]
        ligand_counts[lig] = lig_dict
    return ligand_counts

def check_ligand_distinct(ligands, k):
    num_ligands = len(ligands)
    freq_mat = []
    ligand_counts = ligand_kmer_count(ligands, k)
    for lig in ligand_counts:
        row = ''
        freqs = ligand_counts[lig]
        for kmer in freqs:
            row += str(freqs[kmer])
        freq_mat.append(row)
    unique_mat = set(freq_mat)
    if len(unique_mat) == num_ligands:
        print('Ligands are distinct')
    else:
        print('Ligands are not distinct')
