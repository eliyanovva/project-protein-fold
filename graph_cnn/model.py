import logging as log
import os
import sys

import time
import numpy as np
import pandas as pd
import tensorflow as tf
from sklearn.model_selection import train_test_split
from matplotlib import pyplot as plt
from contextlib import redirect_stdout
import matplotlib.pyplot as plt

# FIXME: VERY UGLY FIX FOR THE config.py
# fix upon packaging and creating the CLI!!!!!
MAIN_PACKAGE_DIR = os.path.abspath(os.curdir)
last_dir = os.path.split(os.curdir)[-1]
while last_dir != "project-protein-fold":
    os.chdir("..")
    MAIN_PACKAGE_DIR = os.path.abspath(os.curdir)
    last_dir = os.path.split(MAIN_PACKAGE_DIR)[-1]
    
sys.path.append(MAIN_PACKAGE_DIR)

import config
from data_prep.data_handlers import DataHandlers
from ligand_handlers import LigandAdjacencyData, LigandFeatureData
from protein_handlers import ProteinAdjacencyData, ProteinFeatureData

log.info(os.getpid())

from tensorflow.keras import backend as K

def coeff_determination(y_true, y_pred):
    SS_res =  K.sum(K.square( y_true-y_pred )) 
    SS_tot = K.sum(K.square( y_true - K.mean(y_true) ) ) 
    return ( 1 - SS_res/(SS_tot + K.epsilon()) )


class GraphCNN:
    def __init__(self):
        pass


    def initialize(self):
        self.prot_adj_data_handler = ProteinAdjacencyData()
        self.prot_feat_data_handler = ProteinFeatureData()
        self.lig_adj_data_handler = LigandAdjacencyData()
        self.lig_feat_data_handler = LigandFeatureData()
        
        log.info('initialized GraphCNN class')


    def trainTestSplit(self, model_test_size=0.3):
        with open(config.MATRIX_DATA_FILES_PATH + '/uniprot_ligand_logfc_pvalue.csv', 'r') as comb_file:
            comb_file_lines = comb_file.readlines()
            x, y = [], []
            # FIX SIZE LATER, SET TO 100 FOR A SMALLER BATCH TRY
        for row in comb_file_lines:
            x_new, y_new = self.__extractRowDataProteinLigandPair(row)
            x.append(x_new)
            y.append(y_new) # logfc score only; doesn't make sense to predict pvalues

        X_train, X_test, y_train, y_test = train_test_split(
            x, y, test_size=model_test_size
            )

        log.info('Split data in training and testing sets.')

        return X_train, X_test, y_train, y_test


    def createModel(self):
        
        prot_adj_in = tf.keras.layers.Input(
            shape=(config.PROTEIN_ADJACENCY_MAT_SIZE,
            config.PROTEIN_ADJACENCY_MAT_SIZE ),
            name='Protein-Adjacency-Matrix'
        )
        
        prot_feat_in = tf.keras.layers.Input(
            shape=(config.PROTEIN_ADJACENCY_MAT_SIZE,
            config.PROTEIN_FEATURES_COUNT),
            name='Protein-Feature-Matrix'
        )

        ligand_adj_in = tf.keras.layers.Input(
            shape=(config.LIGAND_ADJACENCY_MAT_SIZE,
            config.LIGAND_ADJACENCY_MAT_SIZE),
            name='Ligand-Adjacency-Matrix'
        )

        ligand_feat_in = tf.keras.layers.Input(
            shape=(config.LIGAND_ADJACENCY_MAT_SIZE,
            config.LIGAND_FEATURES_COUNT),
            name='Ligand-Feature-Matrix'
        )
 
        x = tf.keras.layers.Conv1D(
            filters=1024, kernel_size=3, activation='relu'
        )(prot_adj_in)
        x = tf.keras.layers.MaxPooling1D(pool_size=(2))(x)
        x = tf.keras.layers.Conv1D(
            filters=512, kernel_size=3, activation='relu'
        )(x)
        x = tf.keras.layers.MaxPooling1D(pool_size=(2))(x)
        x = tf.keras.layers.Conv1D(
            filters=256, kernel_size=3, activation='relu'
        )(x)
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
        #out_classification = tf.keras.layers.Dense(1, activation='softmax')(out)
        
        model = tf.keras.Model(
            inputs=[x.input, y.input, z.input, z1.input],
            outputs= out_regression#, out_classification]
        )
        
        print(model.summary())

        string = model.summary()
        with open('results.txt', 'a') as res_log:
            with redirect_stdout(res_log):
                model.summary()
            res_log.write('\n')


        model.compile(
            optimizer=tf.keras.optimizers.Adagrad(
                learning_rate=0.0001,
                initial_accumulator_value=0.1,
                epsilon=1e-07),
            loss=tf.keras.losses.MeanSquaredLogarithmicError(),
            metrics=[tf.keras.metrics.LogCoshError(),
                coeff_determination,
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


    def __visualizeOutput(self, y):
        log.info('The range for y is [', min(y), ', ', max(y), ']')
        x = np.linspace(0, len(y) - 1, len(y))
        plt.scatter(x, y) 
        plt.savefig('visuals/y_distribution.png')
        plt.close()
        with open(os.path.join(config.MATRIX_DATA_FILES_PATH, 'outputs.csv'), 'w') as res_file:
            res_file.write(','.join([y_mem for y_mem in y]))
        

g = GraphCNN()
g.initialize()

start_train_test_split = time.time()
X_train_labels, X_test_labels, y_train_labels, y_test_labels = g.trainTestSplit()
end_train_test_split = time.time()

start_train_data_load = time.time()
X_train, y_train = g.getTensors(X_train_labels, y_train_labels)
end_train_data_load = time.time()

start_test_data_load = time.time()
X_test, y_test = g.getTensors(X_test_labels, y_test_labels)
end_test_data_load = time.time()

callbacks = [
    tf.keras.callbacks.EarlyStopping(
        monitor = 'val_loss',
        patience = 2,
        restore_best_weights = True,
    )
]

start_model_fitting = time.time()
model = g.createModel()
log.info('model fitting started')
#look into ideal batch size
mod_history = model.fit(X_train, y_train, epochs=10, verbose=True, batch_size=100, callbacks = callbacks, validation_split=0.2)
end_model_fitting = time.time()

log.info ('model fitting finished successfully')

timing_measures = [
    end_train_test_split - start_train_test_split,
    end_train_data_load - start_train_data_load,
    end_test_data_load - start_test_data_load,
    end_model_fitting - start_model_fitting,
    end_model_fitting - start_train_test_split
]

log.info('model evaluation started')
with open('results.txt', 'a') as res_log:
    results = model.evaluate(X_test, y_test, verbose=True)
    res_log.write(' '.join([str(r) for r in results]) + ' \n')
    res_log.write('Timing Benchmarks:\n')
    res_log.write(' '.join([str(r) for r in timing_measures]) + '\n')

print(results)
log.info('model evaluation completed')
#returns loss value (MeanSquaredLogarithmicError) and metric values, currently LogCoshError and coeff_determination
#and RootMeanSquaredError and MeanSquaredError

#plot the loss curve: test vs training

fig, ax1 = plt.subplots(1, figsize=(15, 5))

ax1.plot(mod_history.history["loss"])
ax1.plot(mod_history.history["val_loss"])
ax1.legend(["train", "test"], loc="upper right")
ax1.set_xlabel("Epochs")
ax1.set_ylabel("Loss")

plt.savefig('visuals/loss_graph.png')
plt.close()

# Plot the results
#print(mod_history.history.keys())
#acc = mod_history.history['accuracy']
loss = mod_history.history['loss']
val_loss = mod_history.history['val_loss']
epochs = range(len(loss))

plt.scatter(epochs[1:], loss[1:], c='red', label='Training MSE')
plt.scatter(epochs[1:], val_loss[1:], c='blue', label='Validation MSE')
plt.title('Training accuracy')
plt.legend(loc=0)

plt.savefig('visuals/training_accuracy.png')
plt.close()

