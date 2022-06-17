import numpy as np
import SMILE

#dict ligand_dict ~ key = name of odorant / ligand, value = SMILE formula
ligand_dict = SMILE.create_ligand_dict()

#importmatrix: initializes the global variable ligmat to be a matrix of ligand features
#Input: dict ligand_dict ~ key = name of odorant / ligand, value = SMILE formula
#       int k ~ determines the length of the k-mers
#       int num_proteins ~ the number of protein types used in the dataset
#Output:ligmat ~ a matrix of ligand features 
def importmatrix(ligand_dict, k, num_proteins):
    global ligmat
    ligmat = ligand_matrix(ligand_dict, k, num_proteins)

#ligand_matrix: initializes a matrix of ligand features
#Input: dict ligand_dict
#       int k
#       int num_proteins
#Output: ligmat ~ a matrix of ligand features 
        #each column refers to a k-mer from the SMILE formulas
        #for a dataset with n ligands, rows 1:n each refer to a different ligand
        #data values represent how many times a given k-mer of length k occurs in a ligand's SMILE formula
        #for a dataset with m proteins, rows 1:n of the matrix will be duplicated m times
def ligand_matrix(ligand_dict, k, num_proteins):
    ligand_counts = find_kmer_counts(ligand_dict, k)
    freq_mat = []
    for i in range(num_proteins):
        for lig in ligand_counts:
            freq_mat.append(np.array(list(ligand_counts[lig].values())))
    return np.array(freq_mat)

#Input: dict ligand_dict
#       int k
#Output: dict ligand_counts ~ key = ligand name, value = dict lig_dict:
#                                                               key = k-mer name, value = # of times the k-mer occurs in the ligand
def find_kmer_counts(ligand_dict, k):
    all_kmers = []
    kmers_byligand = {}
    kmer_counts = {}
    for lig in ligand_dict:
        lig_kmers = []
        letters = form_letters(ligand_dict[lig])
        print(letters)
        for i in range(0, len(letters) - k + 1):
            kmer = ""
            for j in range(k):
                kmer += letters[i+j]
            all_kmers.append(kmer)
            lig_kmers.append(kmer)
        kmers_byligand[lig] = lig_kmers

    total_kmers = np.unique(np.array(all_kmers))

    for lig in ligand_dict:
        kmer_counts[lig] = {}
        for kmer in total_kmers:
            kmer_counts[lig][kmer] = kmers_byligand[lig].count(kmer)

    return kmer_counts

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

    return letters
