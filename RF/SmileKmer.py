import numpy as np

#dict ligand_dict ~ key = name of odorant / ligand, value = SMILE formula
ligand_dict = {"pS6_DE_1p_citronellol.csv":'CC(CCC=C(C)C)CCO',
"pS6_DE_1p_isoamylAcetate.csv":'CC(C)CCOC(=O)C',
"pS6_DE_1p_ethylTiglate.csv":'CCOC(=O)C(=CC)C',
"pS6_DE_1p_bIonone.csv":'CC1=C(C(CCC1)(C)C)C=CC(=O)C',
"pS6_DE_1p_butyricAcid.csv":'CCCC(=O)O',
"pS6_DE_1p_paraCresol.csv":'CC1=CC=C(C=C1)O',
"pS6_DE_1p_bCaryophyllene.csv":'CC1=CCCC(=C)C2CC(C2CC1)(C)C',
"pS6_DE_p1_isovalericAcid.csv":'CC(C)CC(=O)O',
"pS6_DE_1p_Octanal.csv":'CCCCCCCC=O',
"pS6_DE_1p_heptanal.csv":'CCCCCCC=O',
"pS6_DE_1p_tbm.csv":'CC(C)(C)S',
"pS6_DE_1p_bDamascone.csv":'CC=CC(=O)C1=C(CCCC1(C)C)C',
"pS6_DE_1p_pyridine.csv":'C1=CC=NC=C1',
"pS6_DE_1p_propionicAcid.csv":'CCC(=O)O',
"pS6_DE_1p_methylSalicylate.csv":'COC(=O)C1=CC=CC=C1O',
"pS6_DE_p01_e2butene1thiol.csv":'CC=CCS',
"pS6_DE_1p_3methyl1butanethiol.csv":'CC(C)CCS',
"pS6_DE_1p_ethylButyrate.csv":'CCCC(=O)OCC',
"pS6_DE_1p_hexylTiglate.csv":'CCCCCCOC(=O)C(=CC)C',
"pS6_DE_1p_indole.csv":'C1=CC=C2C(=C1)C=CN2',
"pS6_DE_500mM_2propylthietane.csv":'CCCC1CCS1',
"pS6_DE_1p_dimethylSulfide.csv":'CSC',
"pS6_DE_1p_2heptanone.csv":'CCCCCC(=O)C',
"pS6_DE_p01_cyclopentanethiol.csv":'C1CCC(C1)S',
"pS6_DE_1p_dimethyltrisulfide.csv":'CSSSC',
"pS6_DE_1p_guaiacol.csv":'COC1=CC=CC=C1O',
"pS6_DE_1p_Benzaldehyde.csv":'C1=CC=C(C=C1)C=O',
"pS6_DE_p01_citral.csv":'CC(=CCCC(=CC=O)C)C',
"pS6_DE_3mM_androstenone.csv":'CC12CCC3C(C1CC=C2)CCC4C3(CCC(=O)C4)C',
"pS6_DE_100p_ebFarnesene.csv":'CC(=CCCC(=CCCC(=C)C=C)C)C',
"pS6_DE_1p_acetophenone.csv":'CC(=O)C1=CC=CC=C1',
"pS6_DE_1p_transCinnamaldehyde.csv":'C1=CC=C(C=C1)C=CC=O',
"pS6_DE_1p_linalool.csv":'CC(=CCCC(C)(C=C)O)C',
"pS6_DE_1p_2hexanone.csv":'CCCCC(=O)C',
"pS6_DE_1p_isopropylTiglate.csv":'CC=C(C)C(=O)OC(C)C',
"pS6_DE_1p_aPinene.csv":'CC1=CCC2CC1C2(C)C',
"pS6_DE_1p_diacetyl.csv":'CC(=O)C(=O)C',
"pS6_DE_1p_geranoil.csv":'CC(=CCCC(=CCO)C)C',
"pS6_DE_1p_heptanoicAcid.csv":'CCCCCCC(=O)O'}
        
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
    ligand_counts = ligand_kmer_count(ligand_dict, k)
    freq_mat = []
    for i in range(num_proteins):
        for lig in ligand_counts:
            freq_mat.append(np.array(list(ligand_counts[lig].values())))
    return np.array(freq_mat)

#Input: dict ligand_dict
#       int k
#Output: dict ligand_counts ~ key = ligand name, value = dict lig_dict:
#                                                               key = k-mer name, value = # of times the k-mer occurs in the ligand
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


#Input: dict ligand_dict
#       int k
#Output: list kmers ~ list of all k-mers of length k that can be found out of all the ligands in ligand_dict
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
