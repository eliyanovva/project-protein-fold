import numpy as np

ligand_dict = {'2hepatanone': 'CCCCCC(=O)C','2hexanone': 'CCCCC(=O)C','3methyl1butanethiol': 'CC(C)CCS',
'acetophenone': 'CC(=O)C1=CC=CC=C1','aPinene': 'CC1=CCC2CC1C2(C)C','bCaryophyllene': 'CC1=CCCC(=C)C2CC(C2CC1)(C)C',
'bDamascone': 'CC=CC(=O)C1=C(CCCC1(C)C)C','Benzaldehyde': 'C1=CC=C(C=C1)C=O','bIonone': 'CC1=C(C(CCC1)(C)C)C=CC(=O)C',
'butyricAcid': 'CCCC(=O)O','citronellol': 'CC(CCC=C(C)C)CCO','diacetyl': 'CC(=O)C(=O)C','dimethylSulfide': 'CSC',
'dimethyltrisulfide': 'CSSSC','ethylButyrate': 'CCCC(=O)OCC','ethylTiglate': 'CCOC(=O)C(=CC)C','geranoil': 'CC(=CCCC(=CCO)C)C',
'guaiacol': 'COC1=CC=CC=C1O','heptanal': 'CCCCCCC=O','heptanoicAcid': 'CCCCCCC(=O)O','hexylTiglate': 'CCCCCCOC(=O)C(=CC)C',
'indole': 'C1=CC=C2C(=C1)C=CN2','isoamyl acetate': 'CC(C)CCOC(=O)C','isopropyl tiglate': 'CC=C(C)C(=O)OC(C)C',
'linalool': 'CC(=CCCC(C)(C=C)O)C','methylSalicylate': 'COC(=O)C1=CC=CC=C1O','Octanal': 'CCCCCCCC=O','paraCresol': 'CC1=CC=C(C=C1)O',
'propionicAcid': 'CCC(=O)O','pyridine': 'C1=CC=NC=C1','tbm': 'CC(C)(C)S','transCinnamaldehyde': 'C1=CC=C(C=C1)C=CC=O',
'androstenone': 'CC12CCC3C(C1CC=C2)CCC4C3(CCC(=O)C4)C','ebFarnesene': 'CC(=CCCC(=CCCC(=C)C=C)C)C','2proplythietane': 'CCCC1CCS1',
'citral': 'CC(=CCCC(=CC=O)C)C','cyclopentanethiol': 'C1CCC(C1)S','e2butene1thiol': 'CC=CCS','isovalericAcid': 'CC(C)CC(=O)O'}
        
 #importmatrix: initializes the global variable ligmat to be a matrix of ligand features
#Input: dict ligand_dict ~ 
#       int k ~ 
#       int num_proteins ~ 
#Output:ligmat ~ a matrix of ligand features 
def importmatrix(ligand_dict, k, num_proteins):
    global ligmat
    ligmat = ligand_matrix(ligand_dict, k, num_proteins)

#ligand_matrix: initializes a matrix of ligand features
def ligand_matrix(ligands, k, num_proteins):
    ligand_counts = ligand_kmer_count(ligands, k)
    freq_mat = []
    for i in range(num_proteins):
        for lig in ligand_counts:
            freq_mat.append(list(ligand_counts[lig].values()))   
    return np.matrix(freq_mat)

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

def find_total_kmers(ligands, k):
    kmers = []
    for lig in ligands:
        k_list = smile_list(ligands[lig], k)
        for kmer in k_list:
            if kmers.count(kmer) == 0:
                kmers.append(kmer)
    return kmers

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
