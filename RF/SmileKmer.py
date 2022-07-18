#This script generates the ligand features for the final matrix

import numpy as np

#Input Variables:
#ligand_dict ~ key: odorant / ligand name, value = SMILE formula
#k: int, length of the ligand k-mers
#num_proteins: int, number of proteins in the dataset
#smile: string, SMILE formula for a given ligand
#ligands: list of ligands

#initializes the global variable ligmat to be a matrix of ligand features
def importmatrix(ligand_dict, k, num_proteins):
    global ligmat
    ligmat = ligand_matrix(ligand_dict, k, num_proteins)

#initializes a matrix of ligand features
def ligand_matrix(ligand_dict, k, ligands):
    ligand_counts = ligand_kmer_count(ligand_dict, k, ligands)

    ligfeatures = list(ligand_counts['pS6_DE_1p_dimethyltrisulfide.csv'].keys())
    return ligfeatures, ligand_counts

def n_ligand_matrix(ligand_dict, k, ligands, filter_kmers):
    ligand_counts = n_ligand_kmer_count(ligand_dict, k, ligands, filter_kmers)
    return ligand_counts

#create a dictionary of the frequency counts for all kmers
#key: ligand, value: dict (key: kmer, value: freq. of kmer in the ligand)
def ligand_kmer_count(ligand_dict, k, ligands):
    ligand_counts = {}
    total_kmers = find_total_kmers(ligand_dict, k, ligands)      #list of kmers found in ALL the ligands
    for lig in ligands:
        lig_dict = {}                                   #freq. dict of ALL kmers for a given ligand
        #ensures that the ligand has a freq. measure for all kmers, even kmers that weren't
        #found in the ligand; leads to easier formatting for the final matrix
        for kmer in total_kmers:
            lig_dict[kmer] = 0
        freq_dict = smile_dict(ligand_dict[lig], k)     #freq. dict of the kmers that DID occur in the ligand
        for kmer in freq_dict:
            lig_dict[kmer] = freq_dict[kmer]            #if the kmer occured in the ligand, lig_dict is updated accordingly
        ligand_counts[lig] = lig_dict
    return ligand_counts


#key: ligand, value: dict (key: kmer, value: freq. of kmer in the ligand)
def n_ligand_kmer_count(ligand_dict, k, ligands, filter_kmers):
    ligand_counts = {}

    for lig in ligands:
        lig_dict = {}                                   #freq. dict of ALL kmers for a given ligand
        #ensures that the ligand has a freq. measure for all kmers, even kmers that weren't
        #found in the ligand; leads to easier formatting for the final matrix
        for kmer in filter_kmers:
            lig_dict[kmer] = 0
        freq_dict = smile_dict(ligand_dict[lig], k)     #freq. dict of the kmers that DID occur in the ligand
        for kmer in filter_kmers:
            if kmer in freq_dict:
                lig_dict[kmer] = freq_dict[kmer]            #if the kmer occured in the ligand, lig_dict is updated accordingly
        ligand_counts[lig] = lig_dict
    return ligand_counts

#creates a list of all kmers that can be found in the ligands
def find_total_kmers(ligand_dict, k, ligands):
    kmers = []
    # iterates thru all ligands
    for lig in ligands:
        k_list = smile_list(ligand_dict[lig], k)    #list of kmers that can be found in a given ligand
        #creates a unique list of the kmers that can be found in all the ligands
        for kmer in k_list:
            if kmers.count(kmer) == 0:
                kmers.append(kmer)
    return kmers

#creates a frequency dictionary based on a ligand's SMILE formula
#key: kmer, value: frequency of the kmer in the ligand
def smile_dict(smile, k):
    kmer_dict = {}              #stores freq. counts for kmers found in the ligand
    letters = form_letters(smile)
    # iterate thru all possible kmers
    for i in range(0, len(letters) - k + 1):
        kmer = ""
        for j in range(k):
            kmer += letters[i+j]
        # update the frequency counts
        if kmer not in kmer_dict:
            kmer_dict[kmer] = 0
        kmer_dict[kmer] += 1

    return kmer_dict

#creates a list of all kmers found in a given ligand's SMILE formula
def smile_list(smile, k):
    kmer_list = []              #stores all kmers found in ligand
    letters = form_letters(smile)
    #iterate thru all possible kmers
    for i in range(0, len(letters) - k + 1):
        kmer = ""
        for j in range(k):
            kmer += letters[i+j]
        #if the kmer hasn't been seen before, add it to kmer_list
        if kmer_list.count(kmer) == 0:
            kmer_list.append(kmer)
            
    return kmer_list
        
#Paritions a SMILE formula into a list of 'letters'
#Each 'letter' is a substring of the SMILE, and can contain more than 1 character
#Each 'letter' will count as a single character in regards to making the kmers

#An overview of SMILE notation can be found at: https://en.wikipedia.org/wiki/Simplified_molecular-input_line-entry_system
# Alphabetic letters designate an atom
# A set of parenthesis '()' designates a side chain
# '=' designates a double bond
# Lowercase letters such as 'c' designate atoms in an aromatic ring

# By our current defintions, letters are:
#   -atoms in the backbone
#   -atoms within an aromatic ring
#   -any bond that isn't a single bond
#   -all atoms within a side chain

def form_letters(smile):
    letters = []    # list that stores the sectioned off 'letters' of the str smile
    for i in range(0, len(smile)):
        letters.append(0)
    a_index = 0     # index of the latest atom
    s_index = 0     # index of the starting point of the latest side chain
    e_index = 0     # index of the ending point of the latest side chain
    mid_index = 0   # index of the latest side chain that picks up again after a nested side chain
    sides = 0       # current number of nested side chains

    for i in range(0, len(smile)):
        if smile[i] == "(":             #indicates a new side chain has begun
            s_index = i
            e_index = i + 2     #The ending parenthesis must be at least two characters away
            letters[i] = "("
            sides += 1
        elif smile[i] == ")":           #indicates a side chain has terminated
            e_index = i
            sides -= 1
            if s_index > mid_index:     #indicates the just-terminated side chain wasn't nested in other side chains
                orig = letters[s_index]
                letters[s_index] = orig + ")"
            else:
                orig = letters[mid_index]   #indicates the current side WAS nested in other side chains
                letters[mid_index] = orig + ")"
        elif smile[i].isdigit():            #checks if the current char is part of the numerical notation for some atom
            if a_index < e_index:           #add numerical notation to a non-nested side chain atom
                orig = letters[s_index]
                letters[s_index] = orig + str(smile[i])
                e_index += 1
            elif sides > 0:                 #add numerical notation to a nested side chain atom
                orig = letters[mid_index]
                letters[mid_index] = orig + str(smile[i])
            else:                           #add numberical notation to backbone atom
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
        
    return letters