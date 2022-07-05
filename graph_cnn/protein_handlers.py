import tensorflow as tf
import numpy as np
import os

import data_prep.constants as constants


class ProteinAdjacencyData:
    def __init__(self):
        #self.protein_name = protein_name
        pass


    def __fetchData(self):
        # the protein name should be from the train/test X data.
        protein_data = os.listdir(constants.PROTEIN_ADJACENCY_PATH)
        for file in protein_data:
            protein_adjacency_matrix = np.load(os.path.join(constants.PROTEIN_ADJACENCY_PATH, file))
        return protein_adjacency_matrix
    

    def createModel(self):
        input_layer = tf.keras.layers.Input()
        conv_layer1 = tf.keras.layers.Conv2D(
            filters=2, kernel_size=3, activation='relu', padding='same'
        )(input_layer)
        pool_layer1 = tf.keras.layers.AveragePooling2D(
            pool_size=(2,2), strides=(1,1), padding='valid'
        )(conv_layer1)
        conv_layer2 = tf.keras.layers.Conv2D(
            filters=2, kernel_size=3, activation='relu', padding='same'
        )(pool_layer1)
        pool_layer2 = tf.keras.layers.AveragePooling2D(
            pool_size=(2,2), strides=(1,1), padding='valid'
        )(conv_layer2)
        flat_layer = tf.keras.layers.Flatten()(pool_layer2)

        adjacency_model = tf.keras.models.Model(
            input_layer,
            conv_layer1,
            pool_layer1,
            conv_layer2,
            pool_layer2,
            flat_layer,
            name='Protein-Adjacency-Model'
        )

        return adjacency_model



class ProteinFeatureData:
    def __init__(self):
        #self.protein_name = protein_name
        pass


    def __fetchData(self):
        # the protein name should be from the train/test X data.
        protein_data = os.listdir(constants.PROTEIN_FEATURE_PATH)
        for file in protein_data:
            protein_feature_matrix = np.load(os.path.join(constants.PROTEIN_FEATURE_PATH, file))
        return protein_feature_matrix
    

    def createModel(self):
        input_layer = tf.keras.layers.Input()
        flat_layer = tf.keras.layers.Flatten()(input_layer)

        feature_model = tf.keras.models.Model(
            input_layer,
            flat_layer,
            name='Protein-Feature-Model'
        )

        return feature_model

class ProteinFeatureDataPDB:
    def __init__(self):
        pass

    
    def __fetchData(self):
        protein_data = os.listdir(constants.PROTEIN_FEATURE_PATH_PDB)
        for file in protein_data:
            protein_feature_matrix = np.load(os.path.join(constants.PROTEIN_FEATURE_PATH_PDB, file))
        return protein_feature_matrix