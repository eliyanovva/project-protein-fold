import logging as log
import os
from stringprep import in_table_c8

import numpy as np
import pandas as pd
import tensorflow as tf
from sklearn.model_selection import train_test_split

import data_prep.constants as constants
import data_prep.log_config
from data_prep.data_handlers import DataHandlers
from ligand_handlers import LigandAdjacencyData, LigandFeatureData
from protein_handlers import ProteinAdjacencyData, ProteinFeatureData

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


    def trainTestSplit(self):
        with open(constants.MATRIX_DATA_FILES_PATH + '/uniprot_ligand_logfc_pvalue.csv', 'r') as comb_file:
            comb_file_lines = self.__getValidEntries(comb_file.readlines())
            x, y = [], []
            # FIX SIZE LATER, SET TO 100 FOR A SMALLER BATCH TRY
        for row in comb_file_lines[:100]:
            row = row[:-1].split(',')
            if len(row) == 4:
                protein_ligand_list = []
                protein_ligand_list.append(row[0]) # protein-ligand(file name from Hiro's lab) tuple
                comp_name_index = row[1].rfind('_')
                protein_ligand_list.append(row[1][comp_name_index + 1:])
                x.append(protein_ligand_list)
                y.append(row[2]) # logfc score only; doesn't make sense to predict pvalues

        X_train, X_test, y_train, y_test = train_test_split(
            x, y, test_size=0.3
            )

        log.info('Split data in training and testing sets.')

        return X_train, X_test, y_train, y_test


    def createModel(self):
        
        prot_adj_in = tf.keras.layers.Input(
            shape=(constants.PROTEIN_ADJACENCY_MAT_SIZE,
            constants.PROTEIN_ADJACENCY_MAT_SIZE ),
            name='Protein-Adjacency-Matrix'
        )
        
        prot_feat_in = tf.keras.layers.Input(
            shape=(constants.PROTEIN_ADJACENCY_MAT_SIZE,
            constants.PROTEIN_FEATURES_COUNT),
            name='Protein-Feature-Matrix'
        )

        ligand_adj_in = tf.keras.layers.Input(
            shape=(constants.LIGAND_ADJACENCY_MAT_SIZE,
            constants.LIGAND_ADJACENCY_MAT_SIZE),
            name='Ligand-Adjacency-Matrix'
        )

        ligand_feat_in = tf.keras.layers.Input(
            shape=(constants.LIGAND_ADJACENCY_MAT_SIZE,
            constants.LIGAND_FEATURES_COUNT),
            name='Ligand-Feature-Matrix'
        )
 
        x = tf.keras.layers.Conv1D(
            2, 3, activation='relu'
        )(prot_adj_in)
        x = tf.keras.layers.MaxPooling1D(pool_size=(2))(x)
        x = tf.keras.layers.Conv1D(
            2, 3, activation='relu'
        )(x)
        x = tf.keras.layers.MaxPooling1D(pool_size=(3))(x)
        
        x = tf.keras.layers.Flatten()(x)
        x = tf.keras.layers.Dense(64, activation="relu")(x)
        x = tf.keras.Model(inputs=prot_adj_in, outputs=x)
        
        y = tf.keras.layers.Flatten()(prot_feat_in)
        y = tf.keras.layers.Dense(64, activation="relu")(y)
        y = tf.keras.layers.Dense(32, activation="relu")(y)
        y = tf.keras.Model(inputs=prot_feat_in, outputs=y)

        z = tf.keras.layers.Flatten()(ligand_adj_in)
        z = tf.keras.layers.Dense(32, activation="relu")(z)
        z = tf.keras.layers.Dense(8, activation="relu")(z)
        z = tf.keras.Model(inputs=ligand_adj_in, outputs=z)
        
        z1 = tf.keras.layers.Flatten()(ligand_feat_in)
        z1 = tf.keras.layers.Dense(16, activation="relu")(z1)
        z1 = tf.keras.layers.Dense(8, activation="relu")(z1)
        z1 = tf.keras.Model(inputs=ligand_feat_in, outputs=z1)

        combined = tf.keras.layers.concatenate([x.output, y.output, z.output, z1.output])
        

        out = tf.keras.layers.Dense(64, activation="relu")(combined)
        out = tf.keras.layers.Dense(1, activation="linear")(out)
        
        model = tf.keras.Model(inputs=[x.input, y.input, z.input, z1.input], outputs=out)

        print(model.summary())

        model.compile(
            optimizer='adam',
            loss=tf.keras.losses.MeanSquaredLogarithmicError(),
            metrics=[tf.keras.metrics.RootMeanSquaredError(), tf.keras.metrics.Accuracy()]
        )
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


    def getTensors(self, X, y):
        X_proteins = [row[0] for row in X]
        X_ligands = [row[1] for row in X]
        self.prot_adj_data_handler.initialize(
            X_proteins, 
            (constants.PROTEIN_ADJACENCY_MAT_SIZE, constants.PROTEIN_ADJACENCY_MAT_SIZE),
            np.int8
        )
        tf_prot_adjacency = self.prot_adj_data_handler.getTensor()

        self.prot_feat_data_handler.initialize(
            X_proteins, 
            (constants.PROTEIN_ADJACENCY_MAT_SIZE, constants.PROTEIN_FEATURES_COUNT),
            float
        )
        tf_prot_features = self.prot_feat_data_handler.getTensor()

        self.lig_adj_data_handler.initialize(
            X_ligands, 
            (constants.LIGAND_ADJACENCY_MAT_SIZE, constants.LIGAND_ADJACENCY_MAT_SIZE),
            np.int8
        )
        tf_lig_adjacency = self.lig_adj_data_handler.getTensor()
        
        self.lig_feat_data_handler.initialize(
            X_ligands, 
            (constants.LIGAND_ADJACENCY_MAT_SIZE, constants.LIGAND_FEATURES_COUNT),
            float
        )
        tf_lig_features = self.lig_feat_data_handler.getTensor()


        tf_logfc = tf.strings.to_number(
            tf.convert_to_tensor(np.array(y)),
            out_type=tf.dtypes.float32
        )
        
        return [tf_prot_adjacency, tf_prot_features, tf_lig_adjacency, tf_lig_features], tf_logfc




g = GraphCNN()
g.initialize()
X_train_labels, X_test_labels, y_train_labels, y_test_labels = g.trainTestSplit()

X_train, y_train = g.getTensors(X_train_labels, y_train_labels)
X_test, y_test = g.getTensors(X_test_labels, y_test_labels)

model = g.createModel()

log.info('model fitting started')

mod_history = model.fit(X_train, y_train, epochs=10, verbose=True, batch_size=1)
log.info ('model fitting finished successfully')

log.info('model evaluation started')
print(model.evaluate(X_test, y_test, verbose=True))
log.info('model evaluation completed')


#import matplotlib.pyplot as plt

# Plot the results
#acc = mod_history.history['accuracy']
#loss = mod_history.history['loss']

#epochs = range(len(acc))

#plt.plot(epochs, acc, 'r', label='Training accuracy')
#plt.plot(epochs, val_acc, 'b', label='Validation accuracy')
#plt.title('Training accuracy')
#plt.legend(loc=0)
#plt.figure()

#plt.savefig('training_accuracy.png')

