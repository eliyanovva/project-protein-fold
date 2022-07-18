import logging as log
from operator import mod
import os
from contextlib import redirect_stdout

import config
import matplotlib.pyplot as plt
import numpy as np
import tensorflow as tf
from matplotlib import pyplot as plt
from sklearn.model_selection import train_test_split

from graph_cnn.ligand_handlers import LigandAdjacencyData, LigandFeatureData
from graph_cnn.protein_handlers import ProteinAdjacencyData, ProteinFeatureData

log.info(os.getpid())

class GraphCNN:
    def __init__(self):
        pass


    def initialize(self):
        self.prot_adj_data_handler = ProteinAdjacencyData()
        self.prot_feat_data_handler = ProteinFeatureData()
        self.lig_adj_data_handler = LigandAdjacencyData()
        self.lig_feat_data_handler = LigandFeatureData()
        
        log.info('initialized GraphCNN class')


    def trainTestSplit(self, model_test_size=0.3, batch_size=-1):
        with open(config.MATRIX_DATA_FILES_PATH + '/uniprot_ligand_logfc_pvalue.csv', 'r') as comb_file:
            comb_file_lines = comb_file.readlines()

            if batch_size == -1:
                batch_size = len(comb_file_lines)
 
            x, y = [], []
        for row in comb_file_lines[:batch_size]:
            x_new, y_new = self.__extractRowDataProteinLigandPair(row)
            x.append(x_new)
            y.append(y_new) # logfc score only; doesn't make sense to predict pvalues

        X_train, X_test, y_train, y_test = train_test_split(
            x, y, test_size=model_test_size
            )

        log.info('Split data in training and testing sets.')

        return X_train, X_test, y_train, y_test


    def createModel(self, hp_optimizer='adagrad'):
        
        prot_adj_in = tf.keras.layers.Input(
            shape=(config.PROTEIN_ADJACENCY_MAT_SIZE, config.PROTEIN_ADJACENCY_MAT_SIZE),
            name='Protein-Adjacency-Matrix'
        )
        
        prot_feat_in = tf.keras.layers.Input(
            shape=(config.PROTEIN_ADJACENCY_MAT_SIZE, config.PROTEIN_FEATURES_COUNT),
            name='Protein-Feature-Matrix'
        )

        ligand_adj_in = tf.keras.layers.Input(
            shape=(config.LIGAND_ADJACENCY_MAT_SIZE, config.LIGAND_ADJACENCY_MAT_SIZE),
            name='Ligand-Adjacency-Matrix'
        )

        ligand_feat_in = tf.keras.layers.Input(
            shape=(config.LIGAND_ADJACENCY_MAT_SIZE, config.LIGAND_FEATURES_COUNT),
            name='Ligand-Feature-Matrix'
        )

        
        x = tf.keras.layers.Conv1D(filters=1024, kernel_size=3, activation='relu')(prot_adj_in)
        x = tf.keras.layers.MaxPooling1D(pool_size=(2))(x)
        x = tf.keras.layers.Conv1D(filters=512, kernel_size=3, activation='relu')(x)
        x = tf.keras.layers.MaxPooling1D(pool_size=(2))(x)
        x = tf.keras.layers.Conv1D(filters=256, kernel_size=3, activation='relu')(x)
        x = tf.keras.layers.MaxPooling1D(pool_size=(2))(x)
        x = tf.keras.layers.Flatten()(x)
        x = tf.keras.layers.Dense(1024, activation="relu")(x)
        x = tf.keras.layers.Dense(512, activation="relu")(x)
        x = tf.keras.Model(inputs=prot_adj_in, outputs=x)
        
        y = tf.keras.layers.Flatten()(prot_feat_in)
        y = tf.keras.layers.Dense(512, activation="relu")(y)
        y = tf.keras.layers.Dense(64, activation="relu")(y)
        y = tf.keras.Model(inputs=prot_feat_in, outputs=y)

        z = tf.keras.layers.Flatten()(ligand_adj_in)
        z = tf.keras.layers.Dense(64, activation="relu")(z)
        z = tf.keras.layers.Dense(16, activation="relu")(z)
        z = tf.keras.Model(inputs=ligand_adj_in, outputs=z)
        
        z1 = tf.keras.layers.Flatten()(ligand_feat_in)
        z1 = tf.keras.layers.Dense(256, activation="relu")(z1)
        z1 = tf.keras.layers.Dense(64, activation="relu")(z1)
        z1 = tf.keras.Model(inputs=ligand_feat_in, outputs=z1)

        combined = tf.keras.layers.concatenate([x.output, y.output, z.output, z1.output])
        
        out = tf.keras.layers.Dense(1024, activation="relu")(combined)
        out = tf.keras.layers.Dense(512, activation="relu")(out)
        out = tf.keras.layers.Dense(64, activation="relu")(out)
        out_regression = tf.keras.layers.Dense(1, activation="linear")(out)
        #FIXME: Add the classification layers
        #out_classification = tf.keras.layers.Dense(1, activation='softmax')(out)
        
        model = tf.keras.Model(
            inputs=[x.input, y.input, z.input, z1.input],
            outputs= out_regression#, out_classification]
        )
        
        with open('results.txt', 'a') as res_log:
            with redirect_stdout(res_log):
                model.summary()
            res_log.write('\n')

        model.compile(
            optimizer=hp_optimizer,
            loss=tf.keras.losses.MeanSquaredLogarithmicError(),
            metrics=[tf.keras.metrics.LogCoshError(),
                tf.keras.metrics.RootMeanSquaredError(),
                tf.keras.metrics.MeanSquaredError()
                ]
        )
        return model


    def getTensors(self, X, y):
        X_proteins = [row[0] for row in X]
        X_ligands = [row[1] for row in X]

        tf_prot_adjacency = self.__getTensor(self.prot_adj_data_handler,
            X_proteins, 
            (config.PROTEIN_ADJACENCY_MAT_SIZE, config.PROTEIN_ADJACENCY_MAT_SIZE),
            np.int8
        )

        tf_prot_features = self.__getTensor(self.prot_feat_data_handler,
            X_proteins, 
            (config.PROTEIN_ADJACENCY_MAT_SIZE, config.PROTEIN_FEATURES_COUNT),
            float
        )

        tf_lig_adjacency = self.__getTensor(self.lig_adj_data_handler,
            X_ligands, 
            (config.LIGAND_ADJACENCY_MAT_SIZE, config.LIGAND_ADJACENCY_MAT_SIZE),
            np.int8
        )
        
        tf_lig_features = self.__getTensor(self.lig_feat_data_handler,
            X_ligands, 
            (config.LIGAND_ADJACENCY_MAT_SIZE, config.LIGAND_FEATURES_COUNT),
            float
        )

        tf_logfc = tf.strings.to_number(
            tf.convert_to_tensor(np.array(y)),
            out_type=tf.dtypes.float32
        )
        
        return [tf_prot_adjacency, tf_prot_features, tf_lig_adjacency, tf_lig_features], tf_logfc


    def __getTensor(self, obj, data, size, data_type):
        obj.initialize(
            data, 
            size,
            data_type
        )
        tf_obj_features = obj.getTensor()
        return tf_obj_features
    
    
    def __extractRowDataProteinLigandPair(self, row):
        protein_ligand_list = []
        logfc = None

        row = row[:-1].split(',')
        if len(row) == 4:
            protein_ligand_list.append(row[0]) # protein-ligand(file name from Hiro's lab) tuple
            comp_name_index = row[1].rfind('_')
            protein_ligand_list.append(row[1][comp_name_index + 1:])
            logfc = row[2]

        return protein_ligand_list, logfc
