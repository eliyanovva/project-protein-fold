#Code adapted from: https://www.datacamp.com/tutorial/random-forests-classifier-python 
import timeit
#imports
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn import metrics
from imblearn.under_sampling import RepeatedEditedNearestNeighbours
import numpy as np
import pandas as pd

def timefunction():
    #Function to create list of protein accessions
    def initialize_protein_list():
        df = pd.read_csv("../data_files/TMdomains/TM.csv")
        protein_list = list(df.iloc[:, 0])
        """
        acc_ids = []
        fr = open("../data_files/AminoAcidSequences/allsequences.fasta", "r")
        lines = fr.readlines()
        for line in lines:
            if line[0] == ">":
                acc_ids.append(line[1:-1])
        fr.close()
        """
        return protein_list
        
        #initializes the global variable ligmat to be a matrix of ligand features
    def importmatrix(ligand_dict, k, num_proteins):
        global ligmat
        ligmat = ligand_matrix(ligand_dict, k, num_proteins)

    #initializes a matrix of ligand features
    def ligand_matrix(ligand_dict, k, num_proteins):
        ligand_counts = ligand_kmer_count(ligand_dict, k)
        freq_mat = []
        for i in range(num_proteins):
            for lig in ligand_counts:
                freq_mat.append(np.array(list(ligand_counts[lig].values())))
        ligfeatures = list(ligand_counts['pS6_DE_1p_citronellol.csv'].keys())
        return np.array(freq_mat) , ligfeatures

    #create a dictionary of the frequency counts for all kmers
    #key: ligand, value: dict (key: kmer, value: freq. of kmer in the ligand)
    def ligand_kmer_count(ligand_dict, k):
        ligand_counts = {}
        total_kmers = find_total_kmers(ligand_dict, k)      #list of kmers found in ALL the ligands
        for lig in ligand_dict:
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

    #creates a list of all kmers that can be found in the ligands
    def find_total_kmers(ligand_dict, k):
        kmers = []
        # iterates thru all ligands
        for lig in ligand_dict:
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
        #Import list of accession numbers
    acc_list = initialize_protein_list()

    #kmers: set of all possible kmers for the protein
    #seqvar: list of Seq class objects (name, sequence = protein sequence, dictionary = kmer freq dictionary)
    #pos_counts: key = protein id, value = # pos. interactions with the protein
    #neg_counts: key = protein id, value = # neg. interactions with the protein
    def richness_protein(kmers, seqvar, pos_counts, neg_counts, domain):
        kmers = list(kmers)

        pos_counts_by_kmer = {}             #key: kmer, value: num. of positive pairs that involve the kmer
        neg_counts_by_kmer = {}             #key: kmer, value: num. of negative pairs that involve the kmer
        pos_prop_by_kmer = {}               #key: kmer, value: proportion of positive pairs that involve the kmer
        neg_prop_by_kmer = {}               #key: kmer, value: proportion of negative pairs that involve the kmer
        counts_by_id = {}
        for kmer in kmers:
            pos_counts_by_kmer[kmer] = 0
            neg_counts_by_kmer[kmer] = 0
            counts_by_id[kmer] = 0

        for kmer in kmers:
            for seq in seqvar:
                dict = seq.dictionary
                if kmer in dict:
                    counts_by_id[kmer] += dict[kmer]

        total_pos = 0                       #total num. of kmers involved in positive pairs
        total_neg = 0                       #total num. of kmers involved in negative pairs

        for seq in seqvar:
            id = seq.name                   #id = accession id of protein
            freq_dict = seq.dictionary      #freq_dict = freq. counts of all known kmers in the protein

            #increase total kmer counts by (num. of pairs that involve the protein) x (num. of kmers that the protein has)
            total_pos += pos_counts[id] * sum(freq_dict.values())
            total_neg += neg_counts[id] * sum(freq_dict.values())

            #increase by_kmer counts by (num. of pairs that involve the protein) x (freq. of a given kmer in the protein)
            for kmer in kmers:
                if (pos_counts[id] > 0) & (kmer in freq_dict):
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
        #the formula for richness has been adapted to account for imbalances in the dataset
        #(ie, basing it off the prop. of a kmer in pos / neg counts, rather than the pure frequency counts for the kmer)
        for kmer in kmers:
            if neg_counts_by_kmer[kmer] == 0:       #kmer only occurs in positive pairs
                richness[kmer] = 100000
            else:
                richness[kmer] = pos_prop_by_kmer[kmer] / neg_prop_by_kmer[kmer]

        ret = []                #list of kmers that meet filtering conditions; to be used in final matrix
        ret2 = []               #for importances list
        lowest_in_top10 = 1000000000
        max = 0
        max_kmer = ""

        for kmer in kmers:

            total_kmer_freq = pos_counts_by_kmer[kmer] + neg_counts_by_kmer[kmer]
            protein_half = len(acc_list) / 2

            if (richness[kmer] <= .125) | (richness[kmer] >= 8):
                ret.append(kmer)
                ret2.append(kmer + domain)
            """

            if counts_by_id[kmer] > max:
                max = counts_by_id[kmer]
                max_kmer = kmer
            """
        return ret, ret2 #, max, max_kmer
        #Creating protein sequence class
    class Seq:
        def __init__(self, name, sequence, dictionary):
            self.name = name
            self.sequence = sequence
            self.dictionary = dictionary #key: k-mer; value: frequency
        def __repr__(self):
            return self.name
        def __eq__(self, other):
            return self.name == other.name

    #function to extract kmers and kmer-frequency from protein sequences
    def featurize(seq,k,feat):
        dict = {}
        for i in range(0, len(seq) - k + 1):
            kmer = ""
            for j in range(k):
                kmer += seq[i + j]

            if kmer not in dict:
                dict[kmer] = 0
                feat.add(kmer)
            dict[kmer] += 1
        return dict

    def make_seqvar_TMS(TM_dict, TM_num, k, seqvar, feat):
        for id in TM_dict:
            name = id
            seq = TM_dict[name][TM_num]
            seqvar.append(Seq(name, seq, featurize(seq, k, feat)))
        return seqvar, feat

    #Return: List of sequence objects representing each protein; List of k-mers found in the protein
    #fasta = file to read
    #seqvar = list to be populated with sequence objects
    #feat = list to be populated with k-mers
    def make_seqvar(fasta, seqvar, feat):
        i = 0
        j = 0
        for line in fasta:
            #Identify and add new sequence objects
            if line[0] == '>':
                name = line.replace('\n','')
                name = name.replace('>', '')
                seqvar.append(Seq(name, '', featurize('',7,feat)))
                i += 1
            #Update sequence objects to store sequences and frequency dictionaries
            else:
                sequence = line.replace('\n','')
                seqvar[i-1].sequence = seqvar[i-1].sequence + sequence
                seqvar[i-1].dictionary = featurize(seqvar[i-1].sequence,7,feat)
            j += 1

        return seqvar, feat

    #Return: 2 dimensional matrix (accession number by k-mer) of frequency values
    #seqvar = list of sequence objects
    #feat = list of k-mers
    #mat = empty matrix to be populated with frequency values
    def makematrix(seqvar, feat, mat):
        for seq in seqvar:
            newseq = []

            for kmer in feat:
                #For kmers not found in the protein, populate the matrix with zeros
                if kmer not in seq.dictionary:
                    seq.dictionary[kmer] = 0
                #Add the frequency value of the kmer
                newseq.append(seq.dictionary.get(kmer))
            #Add a frequency array for each protein
            mat.append(np.array(newseq))

        return mat

    #Function to create a dictionary of ligands matched to SMILES strings
    def initialize_ligand_dict():
        ligand_dict = {}
        df = pd.read_csv('../Ligands_withSMILE/ligand_SMILEs.csv')
        files = df['ligand file'].tolist()
        smiles = df['SMILE'].tolist()
        for i in range(len(files)):
            ligand_dict[files[i]] = smiles[i]

        return ligand_dict

    #List of filenames for ligands that we have matched SMILES strings to
    #Cannot distinguish chirality
    def initialize_ligand_list():
        ligands = ['pS6_DE_1p_citronellol.csv', 'pS6_DE_1p_isoamylAcetate.csv', 'pS6_DE_1p_ethylTiglate.csv',
                'pS6_DE_1p_bIonone.csv', 'pS6_DE_1p_butyricAcid.csv',
                'pS6_DE_1p_paraCresol.csv', 'pS6_DE_1p_bCaryophyllene.csv', 'pS6_DE_p1_isovalericAcid.csv',
                'pS6_DE_1p_Octanal.csv', 'pS6_DE_1p_heptanal.csv',
                'pS6_DE_1p_tbm.csv', 'pS6_DE_1p_bDamascone.csv', 'pS6_DE_1p_pyridine.csv', 'pS6_DE_1p_propionicAcid.csv',
                'pS6_DE_1p_methylSalicylate.csv',
                'pS6_DE_p01_e2butene1thiol.csv', 'pS6_DE_1p_3methyl1butanethiol.csv', 'pS6_DE_1p_ethylButyrate.csv',
                'pS6_DE_1p_hexylTiglate.csv', 'pS6_DE_1p_indole.csv',
                'pS6_DE_500mM_2propylthietane.csv', 'pS6_DE_1p_2heptanone.csv', 'pS6_DE_p01_cyclopentanethiol.csv',
                'pS6_DE_1p_dimethyltrisulfide.csv',
                'pS6_DE_1p_guaiacol.csv', 'pS6_DE_1p_Benzaldehyde.csv', 'pS6_DE_p01_citral.csv',
                'pS6_DE_3mM_androstenone.csv', 'pS6_DE_100p_ebFarnesene.csv',
                'pS6_DE_1p_acetophenone.csv', 'pS6_DE_1p_transCinnamaldehyde.csv', 'pS6_DE_1p_linalool.csv',
                'pS6_DE_1p_2hexanone.csv', 'pS6_DE_1p_isopropylTiglate.csv',
                'pS6_DE_1p_aPinene.csv', 'pS6_DE_1p_diacetyl.csv', 'pS6_DE_1p_geranoil.csv',
                'pS6_DE_1p_heptanoicAcid.csv', 'pS6_DE_1p_2e3mp.csv', 'pS6_DE_1p_2hac.csv', 'pS6_DE_1p_2m2p.csv', 
                'pS6_DE_1p_2m2t.csv', 'pS6_DE_1p_2phenylAlcohol.csv', 'pS6_DE_1p_4methylAC.csv', 'pS6_DE_1p_25dmp.csv',
                'pS6_DE_1p_nCarvone.csv', 'pS6_DE_1p_nDihydrocarveol.csv', 'pS6_DE_1p_nMenthol.csv', 'pS6_DE_1p_ntmt.csv', 
                'pS6_DE_1p_pCarvone.csv', 'pS6_DE_1p_pDihydrocarveol.csv', 'pS6_DE_1p_pLimonene.csv', 'pS6_DE_1p_pMenthol.csv', 'pS6_DE_1p_sbt.csv', 'pS6_DE_1p_tmt.csv']
        return ligands


    def initialize_AA_dict():
        df = pd.read_csv("../data_files/TMdomains/TM.csv")
        protein_list = initialize_protein_list()

        TMs_by_id = {}

        for i in range(len(protein_list)):
            TMs = [str(df.iloc[i, 1]), str(df.iloc[i, 4]), str(df.iloc[i, 7]), str(df.iloc[i, 10])]
            TMs_by_id[protein_list[i]] = TMs

        return categorize(TMs_by_id)

    def initialize_indices():
        df = pd.read_csv("../data_files/TMdomains/TM.csv")
        protein_list = initialize_protein_list()

        TM_indices = {}
        for i in range(len(protein_list)):
            indices = [int(df.iloc[i, 2]), int(df.iloc[i, 3]), int(df.iloc[i, 5]), int(df.iloc[i, 6]), int(df.iloc[i, 8]),
                    int(df.iloc[i, 9]), int(df.iloc[i, 11]), int(df.iloc[i, 12]), ]
            TM_indices[protein_list[i]] = indices

        return TM_indices

    def initialize_3Di_dict():
        TM_indices = initialize_indices()
        Di_dict = {}
        Di = open("../data_files/3DiSequences/fullset_ss.fasta", "r")
        lines = Di.readlines()
        for i in range(len(lines)):
            if i % 2 == 0:
                id = lines[i][1:-1]
                seq = lines[i + 1][:-1]
                Di_dict[id] = seq
        Di.close()

        Di_TMs = {}
        for id in TM_indices:
            TMs = []
            seq = Di_dict[id]

            for i in range(4):
                start = TM_indices[id][2 * i]
                end = TM_indices[id][(2 * i) + 1]
                TMs.append(seq[start - 1:end])

            Di_TMs[id] = TMs

        return Di_TMs

    def categorize(TM_dict):
        categorize_dict = {}
        for id in TM_dict:
            categorize_TMs = []
            for TM in TM_dict[id]:
                TM = TM.replace('A', 'a').replace('G', 'a').replace('V', 'a')
                TM = TM.replace('I', 'b').replace('L', 'b').replace('F', 'b').replace('P', 'b')
                TM = TM.replace('Y', 'c').replace('M', 'c').replace('T', 'c').replace('S', 'c')
                TM = TM.replace('H', 'd').replace('N', 'd').replace('Q', 'd').replace('W', 'd')
                TM = TM.replace('R', 'e').replace('K', 'e')
                TM = TM.replace('D', 'f').replace('E', 'f')
                TM = TM.replace('C', 'g')
                categorize_TMs.append(TM)
            categorize_dict[id] = categorize_TMs
        return categorize_dict

    """
    fr = open('TMs.txt', "r")
    lines = fr.readlines()
    fr.close()

    fw = open('../data_files/TMdomains/TM.csv', 'w')
    fw.write('protein,TM3,s3,e3,TM5,s5,e5,TM6,s6,e6,TM7,s7,e7\n')

    for i in range(len(lines)):
        line = lines[i][:-1]
        #print(line)
        if (i+1) % 6 != 0:
            if line[0] == ">":
                full_line = line[1:] + "," + lines[i+1][:-1] + "," + lines[i+2][:-1]\
                            + "," + lines[i+3][:-1] + "," + lines[i+4][:-1] + "\n"
                fw.write(full_line)
    fw.close()
    """
        #Generate lists of proteins and ligands
    acc_ids = initialize_protein_list()
    csvs = initialize_ligand_list()

    #Initialize variables
    num_proteins = len(acc_ids)
    num_ligands = len(csvs)
    logFC_byID = {}
    FDR_byID = {}
    cit_logFC = {}
    cit_FDR = {}

    #id = accession number of a given protein
    for id in acc_ids:
            logFC_byID[id] = {}
            FDR_byID[id] = {}

    #returns dictionaries of the logFC and p-values
    #key: protein id, value: dict (key: ligand file name, value: data label)
    def labels():
        fas_df = pd.read_csv('uniprot_ensemble.csv', index_col='accession number')
        
        #Read each csv file for the corresponding ligand
        for csv in csvs:
            file_name = '../olfr_de/'+csv
            curr_df = pd.read_csv(file_name, index_col='ensembl_gene_id')

            for id in acc_ids:
                ensem_id = fas_df.loc[id]['ensembl_gene_id'] #The ENSEMBLE id corresponding to the accession number
                logFC_byID[id][csv] = (curr_df.loc[ensem_id]['logFC']) #Find logFC for the ligand-protein pair
                FDR_byID[id][csv] = (curr_df.loc[ensem_id]['FDR']) #Find p-value for the ligand-protein pair

        #Return dictionaries with protein-ligand pair keys and logFC and p-value values
        return logFC_byID, FDR_byID

    #Create a classification dictionary with protein-ligand pair keys and bind (1) or not bind (0) as values
    def classified_logFC_FDR(logFC_byID, FDR_byID):
        classified = {}
        pos_counts = {} #key: protein id, value: number of positive protein interactions
        neg_counts = {} #key: protein id, value: number of negative protein interactions

        class_by_CSV = {}

        for csv in csvs:
            class_by_CSV[csv] = 0

        for id in FDR_byID:
            pos = 0
            neg = 0
            classified[id] = {}
            for csv in csvs:
                if (logFC_byID[id][csv] >= 1) & (FDR_byID[id][csv] <= .05): #The protein and ligand bind
                    classified[id][csv] = 1
                    pos += 1
                    class_by_CSV[csv] += 1
                else: #The protein and ligand do not bind
                    classified[id][csv] = 0
                    neg += 1
            pos_counts[id] = pos
            neg_counts[id] = neg

        return classified, pos_counts, neg_counts

    #Create classification dictionary
    logFC, FDR = labels()
    classified, pos_counts, neg_counts = classified_logFC_FDR(logFC, FDR)

    #Initialize Variables
    #categorized variables
    categorized_features_TM3 = set()
    categorized_seqs_TM3 = []
    categorized_matrix_TM3 = []
    categorized_features_TM5 = set()
    categorized_seqs_TM5 = []
    categorized_matrix_TM5 = []
    categorized_features_TM6 = set()
    categorized_seqs_TM6 = []
    categorized_matrix_TM6 = []
    categorized_features_TM7 = set()
    categorized_seqs_TM7 = []
    categorized_matrix_TM7 = []
    #3Di variables
    di_features_TM3 = set()
    di_seqs_TM3 = []
    di_matrix_TM3 = []
    di_features_TM5 = set()
    di_seqs_TM5 = []
    di_matrix_TM5 = []
    di_features_TM6 = set()
    di_seqs_TM6 = []
    di_matrix_TM6 = []
    di_features_TM7 = set()
    di_seqs_TM7 = []
    di_matrix_TM7 = []

    """
    #Creating output for categorized amino acids
    #Read fasta file
    fasta1 = open("../data_files/AminoAcidSequences/fully_categorized.fasta")
    #Create kmer frequency dictionary
    seqvar1, features1 = ReadingFasta.make_seqvar(fasta1, categorized_seqs, categorized_features)
    #Remove insignificant kmers
    filter_feat = Filtering.richness_protein(features1, seqvar1, pos_counts, neg_counts)
    # Make the matrix
    AA_mat = ReadingFasta.makematrix(seqvar1, filter_feat, categorized_matrix)
    """

    #Create AA output for TMs 3,5,6,7
    AA_dict = initialize_AA_dict()
    AA_seqvar_TM3, AA_features_TM3 = make_seqvar_TMS(AA_dict, 0, 5, categorized_seqs_TM3, categorized_features_TM3)
    AA_seqvar_TM5, AA_features_TM5 = make_seqvar_TMS(AA_dict, 1, 5, categorized_seqs_TM5, categorized_features_TM5)
    AA_seqvar_TM6, AA_features_TM6 = make_seqvar_TMS(AA_dict, 2, 5, categorized_seqs_TM6, categorized_features_TM6)
    AA_seqvar_TM7, AA_features_TM7 = make_seqvar_TMS(AA_dict, 3, 5, categorized_seqs_TM7, categorized_features_TM7)

    AA_filter_TM3, feat1 = richness_protein(AA_features_TM3, AA_seqvar_TM3, pos_counts, neg_counts, "TM3")
    AA_filter_TM5, feat2 = richness_protein(AA_features_TM5, AA_seqvar_TM5, pos_counts, neg_counts, "TM5")
    AA_filter_TM6, feat3 = richness_protein(AA_features_TM6, AA_seqvar_TM6, pos_counts, neg_counts, "TM6")
    AA_filter_TM7, feat4 = richness_protein(AA_features_TM7, AA_seqvar_TM7, pos_counts, neg_counts, "TM7")

    print('AA_TM3 kmers: ' + str(len(AA_filter_TM3)))
    print('AA_TM5 kmers: ' + str(len(AA_filter_TM5)))
    print('AA_TM6 kmers: ' + str(len(AA_filter_TM6)))
    print('AA_TM7 kmers: ' + str(len(AA_filter_TM7)))

    AA_mat_TM3 = makematrix(AA_seqvar_TM3, AA_filter_TM3, categorized_matrix_TM3)
    AA_mat_TM5 = makematrix(AA_seqvar_TM5, AA_filter_TM5, categorized_matrix_TM5)
    AA_mat_TM6 = makematrix(AA_seqvar_TM6, AA_filter_TM6, categorized_matrix_TM6)
    AA_mat_TM7 = makematrix(AA_seqvar_TM7, AA_filter_TM7, categorized_matrix_TM7)

    AA_matrix = np.concatenate((np.array(AA_mat_TM3, dtype = np.uint8), np.array(AA_mat_TM5, dtype = np.uint8),
                                np.array(AA_mat_TM6, dtype = np.uint8), np.array(AA_mat_TM7, dtype = np.uint8)) , axis = 1)

    #318 + 623 + 544 + 375 = 1860
    #183 + 312 + 393 + 280 = 1168

    """
    #Creating output for 3Di sequences
    # Read fasta file
    fasta2 = open("../data_files/3DiSequences/fullset_ss.fasta")
    #Create kmer frequency dictionary
    seqvar2, features2 = ReadingFasta.make_seqvar(fasta2, di_seqs, di_features)
    #Remove insignificant kmers
    filter_feat2 = Filtering.richness_protein(features2, seqvar2, pos_counts, neg_counts)
    # Make the matrix
    Di_mat = ReadingFasta.makematrix(seqvar2, filter_feat2, di_matrix)
    """

    #Create 3Di output for Tms 3,5,6,7
    Di_dict = initialize_3Di_dict()
    Di_seqvar_TM3, Di_features_TM3 = make_seqvar_TMS(Di_dict, 0, 5, di_seqs_TM3, di_features_TM3)
    Di_seqvar_TM5, Di_features_TM5 = make_seqvar_TMS(Di_dict, 1, 5, di_seqs_TM5, di_features_TM5)
    Di_seqvar_TM6, Di_features_TM6 = make_seqvar_TMS(Di_dict, 2, 5, di_seqs_TM6, di_features_TM6)
    Di_seqvar_TM7, Di_features_TM7 = make_seqvar_TMS(Di_dict, 3, 5, di_seqs_TM7, di_features_TM7)

    Di_filter_TM3, feat5 = richness_protein(Di_features_TM3, Di_seqvar_TM3, pos_counts, neg_counts, "TM3")
    Di_filter_TM5, feat6 = richness_protein(Di_features_TM5, Di_seqvar_TM5, pos_counts, neg_counts, "TM5")
    Di_filter_TM6, feat7 = richness_protein(Di_features_TM6, Di_seqvar_TM6, pos_counts, neg_counts, "TM6")
    Di_filter_TM7, feat8 = richness_protein(Di_features_TM7, Di_seqvar_TM7, pos_counts, neg_counts, "TM7")

    print('Di_TM3 kmers: ' + str(len(Di_filter_TM3)))
    print('Di_TM5 kmers: ' + str(len(Di_filter_TM5)))
    print('Di_TM6 kmers: ' + str(len(Di_filter_TM6)))
    print('Di_TM7 kmers: ' + str(len(Di_filter_TM7)))

    Di_mat_TM3 = makematrix(Di_seqvar_TM3, Di_filter_TM3, di_matrix_TM3)
    Di_mat_TM5 = makematrix(Di_seqvar_TM5, Di_filter_TM5, di_matrix_TM5)
    Di_mat_TM6 = makematrix(Di_seqvar_TM6, Di_filter_TM6, di_matrix_TM6)
    Di_mat_TM7 = makematrix(Di_seqvar_TM7, Di_filter_TM7, di_matrix_TM7)

    Di_matrix = np.concatenate((np.array(Di_mat_TM3, dtype = np.uint8), np.array(Di_mat_TM5, dtype = np.uint8),
                                np.array(Di_mat_TM6, dtype = np.uint8), np.array(Di_mat_TM7, dtype = np.uint8)) , axis = 1)

    #Concatenate AA and 3Di matrices
    #intermed_matrix = np.concatenate((np.array(AA_mat, dtype = np.uint8), np.array(Di_mat, dtype = np.uint8)) , axis = 1)
    intermed_matrix = np.concatenate((np.array(AA_matrix, dtype = np.uint8), np.array(Di_matrix, dtype = np.uint8)) , axis = 1)
    #Expand the protein matrix to account for the ligands
    ligand_count = 55
    proteins_matrix = np.repeat(intermed_matrix, repeats = ligand_count, axis = 0)

    #Import dictionary matching ligands to SMILES String
    ligand_dict = initialize_ligand_dict()
    #Create ligands matrix
    #ligand_matrix, ligand_features = SmileKmer.ligand_matrix(ligand_dict, 5, 1084)
    ligand_matrix, ligand_features = ligand_matrix(ligand_dict, 5, 796)

    #Concatenate protein and ligand matrices
    final_matrix = np.concatenate((proteins_matrix, np.array(ligand_matrix, dtype = np.uint8)), axis = 1)

    #Create Classification Vector
    #proteins = seqvar1
    proteins = initialize_protein_list()
    logFCmat = []
    for protein in proteins:
        for ligand in list(ligand_dict.keys()):
            logFCmat.append(float(classified[protein][ligand]))

#Fixed Classification Model 
    def train(features, labels):
        #define features and labels
        X = features #Kmers
        y = labels #Binds or not

        #split into training and test set
        X_train, X_test, y_train, y_test = train_test_split(X, y, stratify=y,test_size=0.1) # 90% training and 10% test

        #Oversampling was necessary, because most ligand/receptor pairs do not bind in our dataset
        enn = RepeatedEditedNearestNeighbours()

        X_res, y_res = enn.fit_resample(X_train, y_train)

        #Create a Gaussian Regression
        clf=RandomForestClassifier(n_estimators=100)

        #Train the model
        clf.fit(X_res,y_res)

        #Form predictions
        y_pred=clf.predict_proba(X_test)[:,1]

        precision, recall, thresholds = metrics.precision_recall_curve(y_test, y_pred)

        acc = metrics.roc_auc_score(y_test, y_pred)
        rec = metrics.auc(recall,precision)

        #Print accuracy of the model
        print("Accuracy:",acc)
        print("Recall:",rec)

        return acc,rec

    #Train and Test

    acc,rec = train(final_matrix, logFCmat)

t = timeit.timeit(lambda: timefunction(), number = 10, setup = """
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn import metrics
from imblearn.under_sampling import RepeatedEditedNearestNeighbours
import numpy as np
import pandas as pd""")
print(t)