import tensorflow as tf
import numpy as np
from sklearn.model_selection import train_test_split
import os

import data_prep.constants as constants
# train-test split upon uniprot-ligand-logfc csv

class GraphCNN:
    def __init__(self):
        pass
    
    
    def trainTestSplit(self):
        ### TODO: Remove entries in the comb_file_lines for which we don't have the protein or the ligand!!!
        with open('uniprot_ligand_logfc_pvalue.csv', 'r') as comb_file:
            comb_file_lines = comb_file.readlines()
            x, y = [], []
        for row in comb_file_lines:
            row = row[:-1].split(',')
            if len(row) == 4:
                protein_ligand_list = []
                protein_ligand_list.append(row[0]) # protein-ligand(file name from Hiro's lab) tuple
                comp_name_index = row[1].rfind('_')
                protein_ligand_list.append(row[1][comp_name_index + 1:])
                x.append(protein_ligand_list)
                y.append(row[2]) # logfc score only; doesn't make sense to predict pvalues

        X_train, X_test, y_train, y_test = train_test_split(x, y, test_size=0.3)
        return X_train, X_test, y_train, y_test


    def __fetchProteinData(self, protein_name):
        # the protein name should be from the train/test X data.
        protein_adjacency_matrix = np.load(os.path.join(constants.PROTEIN_ADJACENCY_PATH, protein_name + '_adj_mat.npy'))
        protein_feature_matrix = np.load(os.path.join(constants.PROTEIN_FEATURE_PATH, protein_name + '_feat_mat.npy'))
        return protein_adjacency_matrix, protein_feature_matrix

    def __fetchLigandData(self, ligand_file_name):
        ligand_file_names = os.listdir(constants.MOL_ADJACENCY_PATH)
        data_file_name = self.__getLigandName(ligand_file_names, ligand_file_name)
        ligand_adjacency_matrix = np.load(os.path.join(constants.MOL_ADJACENCY_PATH, data_file_name))
        return ligand_adjacency_matrix


    def __getLigandName(self, files_list, ligand_name):
        res = None
        print(ligand_name)
        for name in files_list:
            #print(name)
            if ligand_name in name:
                print('TUK SUM!!!!!!!')
                res = name
                
        return res
        
g = GraphCNN()
a, b, c, d = g.trainTestSplit()

print(len(a))
print(g.fetchProteinData(a[4000][0]))

#print(b)
#print(d)
#with open('testf.txt', 'a') as test_file:
##    test_file.write(','.join(','.join(b)))
  #  test_file.write(','.join(d))

