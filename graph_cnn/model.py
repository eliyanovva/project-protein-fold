import tensorflow as tf
import numpy as np
import pandas as pd
import logging as log
from sklearn.model_selection import train_test_split
from data_prep.data_handlers import DataHandlers
import os

import data_prep.log_config
from protein_handlers import ProteinAdjacencyData, ProteinFeatureData
from ligand_handlers import LigandAdjacencyData
import data_prep.constants as constants

from tensorflow.compat.v1 import ConfigProto
from tensorflow.compat.v1 import InteractiveSession

config = ConfigProto()
config.gpu_options.per_process_gpu_memory_fraction = 0.4
session = InteractiveSession(config=config)


log.info(os.getpid())

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
        with open(constants.MATRIX_DATA_FILES_PATH + '/uniprot_ligand_logfc_pvalue.csv', 'r') as comb_file:
            comb_file_lines = self.__getValidEntries(comb_file.readlines())
            x, y = [], []
            # FIX SIZE LATER, SET TO 250 FOR A SMALLER BATCH TRY
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
        #log.info('Split data in training and testing sets.')
        return X_train, X_test, y_train, y_test


    def createModel(self):
        #prot_adj_model = self.prot_adj_data_handler.createModel()
        #prot_feat_model = self.prot_feat_data_handler.createModel()
        #lig_adj_model = self.lig_adj_data_handler.createModel()
        
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
 
        # the first branch operates on the first input
        x = tf.keras.layers.Conv1D(
            2, 3, activation='relu'
        )(prot_adj_in)
        x = tf.keras.layers.AveragePooling1D(pool_size=(2))(x)
        x = tf.keras.layers.Flatten()(x)
        x = tf.keras.layers.Dense(128, activation="relu")(x)
        x = tf.keras.Model(inputs=prot_adj_in, outputs=x)
        # the second branch opreates on the second input
        
        y = tf.keras.layers.Flatten()(prot_feat_in)
        y = tf.keras.layers.Dense(32, activation="relu")(y)
        #y = tf.keras.layers.Dense(4, activation="relu")(y)
        y = tf.keras.Model(inputs=prot_feat_in, outputs=y)

        z = tf.keras.layers.Flatten()(ligand_adj_in)
        z = tf.keras.layers.Dense(32, activation="relu")(z)
        #z = tf.keras.layers.Dense(4, activation="relu")(z)
        z = tf.keras.Model(inputs=ligand_adj_in, outputs=z)
        # combine the output of the two branches
        combined = tf.keras.layers.concatenate([x.output, y.output, z.output])
        # apply a FC layer and then a regression prediction on the
        # combined outputs
        out = tf.keras.layers.Dense(2, activation="relu")(combined)
        out = tf.keras.layers.Dense(1, activation="linear")(out)
        # our model will accept the inputs of the two branches and
        # then output a single value
        model = tf.keras.Model(inputs=[x.input, y.input, z.input], outputs=out)

        print(model.summary())
        #concat_layer = tf.keras.layers.Concatenate([
        #    prot_adj_model.output,
        #    prot_feat_model.output,
        #    lig_adj_model.output
       # ])

        #output = tf.keras.layers.Dense(1024, activation='relu')(concat_layer)
        #output = tf.keras.layers.Dense(512, activation='relu')(output)
        #output = tf.keras.layers.Dense(128, activation='relu')(output)
        #output = tf.keras.layers.Dense(1, activation='relu')(output)    
        
        #model = tf.keras.models.Model([
        #    prot_adj_model.input,
        #    prot_feat_model.input,
        #    lig_adj_model.input
        #    ], 
        #    outputs=output
        #)

        model.compile(optimizer='adam', loss='mean_squared_error', metrics=['accuracy'])
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

np_prot_adjacency = np.zeros((
    len(X_train),
    constants.PROTEIN_ADJACENCY_MAT_SIZE,
    constants.PROTEIN_ADJACENCY_MAT_SIZE
    ), dtype = 'int8'
)

np_prot_features = np.zeros((
    len(X_train),
    constants.PROTEIN_ADJACENCY_MAT_SIZE,
    constants.PROTEIN_FEATURES_COUNT)
)

np_ligand_adjacency = np.zeros((
    len(X_train),
    constants.LIGAND_ADJACENCY_MAT_SIZE,
    constants.LIGAND_ADJACENCY_MAT_SIZE), dtype='int8'
)

np_logfc = np.array(y_train)

log.info('Began initializing the arrays')
for i in range(len(X_train)):
    dh = DataHandlers()

    log.info('NP arrays start loading')
    protein_data = dh.fetchProteinData(X_train[i][0])
    ligand_data = dh.fetchLigandData(X_train[i][1])
    log.info('np arrays loading finished')

    protein_adjacency_matrix = np.pad(
        np.array(protein_data[0], ndmin=2),
        ((0, constants.PROTEIN_ADJACENCY_MAT_SIZE - len(protein_data[0])),
        (0, constants.PROTEIN_ADJACENCY_MAT_SIZE - len(protein_data[0][0]))),
        'constant',
        constant_values=(0)
    )

    protein_feature_matrix = np.pad(
        np.array(protein_data[1], ndmin=2),
        ((0, constants.PROTEIN_ADJACENCY_MAT_SIZE - len(protein_data[1])),
        (0, constants.PROTEIN_FEATURES_COUNT - len(protein_data[1][0]))),
        'constant',
        constant_values=(0)
    )

    ligand_adjacency_matrix = np.pad(
        np.array(ligand_data, ndmin=2),
        ((0, constants.LIGAND_ADJACENCY_MAT_SIZE - len(ligand_data)),
        (0, constants.LIGAND_ADJACENCY_MAT_SIZE - len(ligand_data[0]))),
        'constant',
        constant_values=(0)
    )
        
    np_prot_adjacency[i] = protein_adjacency_matrix
    np_prot_features[i] = protein_feature_matrix
    np_ligand_adjacency[i] = ligand_adjacency_matrix

    log.info('Done with initializing the arrays ')

log.info('start converting numpy data to tensors')
tf_prot_adjacency = tf.convert_to_tensor(
    np_prot_adjacency
    #shape=[len(X_train), constants.PROTEIN_ADJACENCY_MAT_SIZE, constants.PROTEIN_ADJACENCY_MAT_SIZE],
    #dtype=int
)
tf_prot_features = tf.convert_to_tensor(
    np_prot_features
    #shape=(len(X_train), constants.PROTEIN_ADJACENCY_MAT_SIZE, constants.PROTEIN_FEATURES_COUNT),
    #dtype=float
)
tf_ligand_adjacency = tf.convert_to_tensor(
    np_ligand_adjacency
#    shape=(len(X_train), constants.LIGAND_ADJACENCY_MAT_SIZE, constants.LIGAND_ADJACENCY_MAT_SIZE),
#    dtype=int
)
print(y_train[0])

tf_logfc = tf.strings.to_number(tf.convert_to_tensor(np_logfc), out_type=tf.dtypes.float32)

print(tf_ligand_adjacency.dtype, tf_logfc.dtype, tf_prot_adjacency.dtype, tf_prot_features.dtype)

log.info('completed tensor generation')
model = g.createModel()
#print(X_train)
#model = g.createModel()

log.info('model fitting started')

model.fit([tf_prot_adjacency, tf_prot_features, tf_ligand_adjacency],
          tf_logfc,
          epochs=10, verbose=True, batch_size=1)

log.info ('model fitting finished successfully')

#log.info('model evaluation started')
#print(model.evaluate([X_test.iloc[:,0], X_test.iloc[:,1], X_test.iloc[:,2]],
#          y_test, verbose=False))
#log.info('model evaluation completed')"""
#print(len(a))
#print(g.fetchProteinData(a[4000][0]))

#print(b)
#print(d)
#with open('testf.txt', 'a') as test_file:
##    test_file.write(','.join(','.join(b)))
  #  test_file.write(','.join(d))

