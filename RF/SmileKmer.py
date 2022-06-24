import numpy as np
import Globals

#dict ligand_dict ~ key = name of odorant / ligand, value = SMILE formula
ligand_dict = Globals.initialize_ligand_dict()

#initializes the global variable ligmat to be a matrix of ligand features
#dict ligand_dict ~ key = ligand name, value = SMILE formula
#int k ~ determines the length of the k-mers
#int num_proteins ~ the number of protein types used in the dataset
#Output:ligmat ~ a matrix of ligand features 
def importmatrix(ligand_dict, k, num_proteins):
    global ligmat
    ligmat = ligand_matrix(ligand_dict, k, num_proteins)

#ligand_matrix: initializes a matrix of ligand features
def ligand_matrix(ligand_dict, k, num_proteins):
    ligand_counts = ligand_kmer_count(ligand_dict, k)
    freq_mat = []
    for i in range(num_proteins):
        for lig in ligand_counts:
            freq_mat.append(np.array(list(ligand_counts[lig].values())))
    return np.array(freq_mat)

#create a dictionary of the frequency counts for all kmers
def ligand_kmer_count(ligand_dict, k):
    ligand_counts = {}
    total_kmers = find_total_kmers(ligand_dict, k)
    for lig in ligand_dict:
        lig_dict = {}
        for kmer in total_kmers:
            lig_dict[kmer] = 0
        freq_dict = smile_dict(ligand_dict[lig], k)
        for kster in freq_dict:
            lig_dict[kster] = freq_dict[kster]
        ligand_counts[lig] = lig_dict
    return ligand_counts


#create a list of all kmers that can be found in the ligands
def find_total_kmers(ligand_dict, k):
    kmers = []
    for lig in ligand_dict:
        k_list = smile_list(ligand_dict[lig], k)
        for kmer in k_list:
            if kmers.count(kmer) == 0:
                kmers.append(kmer)
    return kmers

def smile_dict(smile, k):
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
        
#Input: str smile = a SMILE formula for a given ligand
#Output: list letters = a list of substrings of smile; each substring is a partitioned 'letter' of smile that can be used to form k-mers
def form_letters(smile):        
        # Letters are: 
        # ~ atoms in the backbone
        # ~ any bond that isn't a single bond
        # ~ all atoms within a side chain
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

    """
    for j in range(len(letters)):
        curr = letters[j]
        ind_1 = curr.find("(")
        ind_2 = curr.find(")")
        if ind_1 != -1 & ind_2 != -1:
            letters[j] = curr[ind_1 + 1: ind_2]
        elif ind_1 != -1:
            letters[j] = curr[ind_1 + 1:]
        elif ind_2 != -1:
            letters[j] = curr[:ind_2]
    """
    return letters

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
