#This script generates the ligand features for the final matrix

import numpy as np

#Additional coding help from:
#https://www.geeksforgeeks.org/remove-all-the-occurrences-of-an-element-from-a-list-in-python/
#https://www.geeksforgeeks.org/python-count-occurrences-element-list/
#https://www.w3schools.com/python/ref_list_sort.asp
#https://www.geeksforgeeks.org/python-convert-set-into-a-list/

def ligand_kmer_count(ligand_dict, k, Ligands):
    """
    This functions creates a dictionary of the frequency counts for all kmers.
    Args:
        ligand_dict (dict): dictionary mapping a (str) ligand to its (str) SMILE formula
        k (int): desired kmer length
        ligands (list): ligands to use in the matrix

    Returns:
        total_kmers (list): sorted list of kmers found in every ligand from Ligands
        ligand_counts (dict): dictionary mapping a (str) ligand to a frequency dictionary
            ex: ligand_counts[lig]: key = (str) kmer, value = (int) freq. of kmer in lig
    """

    ligand_counts = {}
    total_kmers = list(find_total_kmers(ligand_dict, k, Ligands))
    total_kmers.sort()
    for lig in Ligands:
        lig_dict = {}                                   #freq. dict of ALL kmers for a given ligand
        #ensures that the ligand has a freq. measure for all kmers, even kmers that weren't
        #found in the ligand; leads to easier formatting for the final matrix
        for kmer in total_kmers:
            lig_dict[kmer] = 0
        freq_dict = smile_dict(ligand_dict[lig], k)     #freq. dict of the kmers that DID occur in the ligand
        for kmer in freq_dict:
            lig_dict[kmer] = freq_dict[kmer]            #if the kmer occured in the ligand, lig_dict is updated accordingly
        ligand_counts[lig] = lig_dict
    return total_kmers, ligand_counts

def find_total_kmers(ligand_dict, k, Ligands):
    """
    This function creates a set of all kmers found from every ligand in Ligands.
    Args:
        ligand_dict (dict): dictionary mapping a (str) ligand to its (str) SMILE formula
        k (int): desired kmer length
        Ligands (list): ligands to use in the matrix

    Returns:
        total_kmers (set): set of kmers found from every ligand in Ligands
    """
    total_kmers = set()
    # iterates thru all ligands
    for lig in Ligands:
        k_list = smile_list(ligand_dict[lig], k)    #list of kmers that can be found in a given ligand
        #creates a unique list of the kmers that can be found in all the ligands
        for kmer in k_list:
            if kmer not in total_kmers:
                total_kmers.add(kmer)
    return total_kmers

def smile_dict(smile, k):
    """
    This function creates a frequency dictionary of kmers found within smile
    Args:
        smile (string): SMILE formula
        k (int): desired kmer length

    Returns:
        kmer_dict (dictionary): key = (string) kmer, value = (int) freq. of kmer in smile
    """
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
    """
    This function creates a list of all kmers found within smile.
    Args:
        smile (string): SMILE formula
        k (int): desired kmer length

    Returns:
        smile_list (list): list of all kmers found within smile
    """
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
#   -all atoms within a set of brackets
#   -all atoms within a side chain

def form_letters(smile):
    """
    This function partitions smile into a list of 'letters' to be used while forming kmers.
    Args:
        smile (str): SMILE formula

    Returns:
        letters (list): list of partioned substrings ('letters') from smile
    """
    letters = []  # list that stores the sectioned off 'letters' of the str smile
    for i in range(0, len(smile)):
        letters.append(0)
    a_index = 0  # index of the latest atom
    s_index = 0  # index of the starting point of the latest side chain
    e_index = 0  # index of the ending point of the latest side chain
    mid_index = 0  # index of the latest side chain that picks up again after a nested side chain
    sides = 0  # current number of nested side chains
    b_start = 0     # index of the starting point of the latest bracket
    b_end = 0       # index of the ending point of the latest bracket

    for i in range(0, len(smile)):
        if smile[i] == "(":  # indicates a new side chain has begun
            s_index = i
            e_index = i + 2  # The ending parenthesis must be at least two characters away
            letters[i] = "("
            sides += 1
        elif smile[i] == ")":  # indicates a side chain has terminated
            e_index = i
            sides -= 1
            if s_index > mid_index:  # indicates the just-terminated side chain wasn't nested in other side chains
                orig = letters[s_index]
                letters[s_index] = orig + ")"
            else:
                orig = letters[mid_index]  # indicates the current side WAS nested in other side chains
                letters[mid_index] = orig + ")"

        elif smile[i].isdigit():  # checks if the current char is part of the numerical notation for some atom
            if a_index < e_index:  # add numerical notation to a non-nested side chain atom
                orig = letters[s_index]
                letters[s_index] = orig + str(smile[i])
                e_index += 1
            elif sides > 0:  # add numerical notation to a nested side chain atom
                orig = letters[mid_index]
                letters[mid_index] = orig + str(smile[i])
            elif smile[i-1] == ']':
                orig = letters[b_start]
                letters[b_start] = orig + str(smile[i])
            else:  # add numberical notation to backbone atom
                orig = letters[a_index]
                letters[a_index] = orig + str(smile[i])
        else:
            a_index = i
            if sides > 0:
                if a_index < e_index:  # current atom is part of a non-nested side chain
                    orig = letters[s_index]
                    letters[s_index] = orig + smile[i]
                    e_index += 1
            #elif sides > 0:  # current atom is following a nested side chain
                elif smile[i - 1] == ")":
                    letters[i] = (smile[i])
                    mid_index = i
                else:
                    orig = letters[mid_index]
                    letters[mid_index] = orig + smile[i]

            if sides == 0:
                if smile[i] == '[':
                    b_start = i
                    letters[b_start] = '['
                    b_end = i + 2
                elif smile[i] == ']':
                    b_end = i
                    orig = letters[b_start]
                    letters[b_start] = orig + ']'
                elif i < b_end:
                    b_end += 1
                    orig = letters[b_start]
                    letters[b_start] = orig + smile[i]

                else:  # current atom is part of the backbone
                    letters[i] = (smile[i])

    # remove indices from letters that were not used
    c = letters.count(0)
    for i in range(0, c):
        letters.remove(0)

    return letters