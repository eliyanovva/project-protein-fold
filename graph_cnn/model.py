import tensorflow as tf
import numpy as np
import pandas as pd
import logging as log
from sklearn.model_selection import train_test_split
import os

import data_prep.log_config
from protein_handlers import ProteinAdjacencyData, ProteinFeatureData
from ligand_handlers import LigandAdjacencyData
import data_prep.constants as constants


class GraphCNN:
    def __init__(self):
        self.prot_adj_data_handler = ProteinAdjacencyData()
        self.prot_feat_data_handler = ProteinFeatureData()
        self.lig_adj_data_handler = LigandAdjacencyData()
        
        log.info('initialized GraphCNN class')
        self.dataset = pd.read_pickle(os.path.join(constants.MATRIX_DATA_FILES_PATH, 'dataset_part_0.pkl'))
        log.info('Loaded dataset')
        self.y = self.dataset.iloc[:,3]
        log.info('Loaded target columns')
        self.dataset = self.dataset.drop(self.dataset.columns[[3]],axis = 1)
        log.info('Removed y column from dataset')
        
        
#        self.dataset = tf.data.experimental.make_csv_dataset(
#            './dataset.csv',
#            column_names=['ProteinAdjacency', 'ProteinFeatures', 'LigandAdjacency', 'logFC'],
#            batch_size=4, # Artificially small to make examples easier to show.
#            label_name='logFC',
#            num_epochs=1,
#            header=True
            #ignore_errors=True
#            )
        #chunksize = 10 ** 6
        #textFileReader = pd.read_csv('dataset.csv', chunksize=chunksize)
        #full_data = pd.concat(textFileReader, ignore_index=True)
        #self.dataset = pd.read_csv('dataset.csv')
        #print('finished reading csv')
        #self.y = self.dataset.iloc[:, 2]
        #self.dataset.drop([2], axis=1)
    
    
    def trainTestSplit(self):
        #with open('uniprot_ligand_logfc_pvalue.csv', 'r') as comb_file:
        #    comb_file_lines = self.__getValidEntries(comb_file.readlines())
        #    x, y = [], []
        #for row in comb_file_lines:
        #    row = row[:-1].split(',')
        #    if len(row) == 4:
        #        protein_ligand_list = []
        #        protein_ligand_list.append(row[0]) # protein-ligand(file name from Hiro's lab) tuple
        #        comp_name_index = row[1].rfind('_')
        #        protein_ligand_list.append(row[1][comp_name_index + 1:])
        #        x.append(protein_ligand_list)
        #        y.append(row[2]) # logfc score only; doesn't make sense to predict pvalues

        X_train, X_test, y_train, y_test = train_test_split(
            self.dataset, self.y, test_size=0.3
            )
        log.info('Split data in training and testing sets.')
        return X_train, X_test, y_train, y_test


    def createModel(self):
        #undid '#' infront of these three lines
        prot_adj_model = self.prot_adj_data_handler.createModel()
        prot_feat_model = self.prot_feat_data_handler.createModel()
        lig_adj_model = self.lig_adj_data_handler.createModel()
        
        prot_adj_in = tf.keras.layers.Input(
            shape=(constants.PROTEIN_ADJACENCY_MAT_SIZE, ),
            name='Protein-Adjacency-Matrix'
        )
        
        prot_feat_in = tf.keras.layers.Input(
            shape=(constants.PROTEIN_ADJACENCY_MAT_SIZE, ),
            name='Protein-Feature-Matrix'
        )

        ligand_adj_in = tf.keras.layers.Input(
            shape=(constants.LIGAND_ADJACENCY_MAT_SIZE, ),
            name='Ligand-Adjacency-Matrix'
        )


        concat_layer = tf.keras.layers.Concatenate([
            prot_adj_model.output,
            prot_feat_model.output,
            lig_adj_model.output
        ])

        output = tf.keras.layers.Dense(1024, activation='relu')(concat_layer)
        output = tf.keras.layers.Dense(512, activation='relu')(output)
        output = tf.keras.layers.Dense(128, activation='relu')(output)
        output = tf.keras.layers.Dense(1, activation='relu')(output)    
        
        model = tf.keras.models.Model([
            prot_adj_model.input,
            prot_feat_model.input,
            lig_adj_model.input
            ], 
            outputs=output
        )

        model.compile(optimizer='adam', loss='mean_absolute_error')
        return model


    def __getValidEntries(self, comb_file_lines):
        # ensures that we have all necessary mol and pdb files for the training and testing data.
        all_protein_files = os.listdir(constants.PROTEIN_ADJACENCY_PATH)
        all_ligand_files = os.listdir(constants.MOL_ADJACENCY_PATH)

        for line in list(comb_file_lines):
            line_list = line.split(',')
            
            protein_match = False
            ligand_match = False
            for protein in all_protein_files:
                if line_list[0] in protein:
                    protein_match = True
                    break
            
            for ligand in all_ligand_files:
                comp_name_index = line_list[1].rfind('_')
                
                if line_list[1][comp_name_index + 1 :] in ligand:
                    ligand_match = True
                    break
            if not (protein_match and ligand_match):
                comb_file_lines.remove(line)
        return comb_file_lines


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
        for name in files_list:
            if ligand_name in name:
                res = name
                
        return res
        
g = GraphCNN()
X_train, X_test, y_train, y_test = g.trainTestSplit()
print(y_test)
#model = g.createModel()
#print(len(a))
#print(g.fetchProteinData(a[4000][0]))

#print(b)
#print(d)
#with open('testf.txt', 'a') as test_file:
##    test_file.write(','.join(','.join(b)))
  #  test_file.write(','.join(d))

