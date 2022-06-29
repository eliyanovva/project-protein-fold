import tensorflow as tf
import numpy as np
import os

import data_prep.constants as constants


class LigandAdjacencyData:
    def __init__(self):
        #self.protein_name = protein_name
        pass


    def __fetchData(self):
        # the protein name should be from the train/test X data.
        ligand_data = os.listdir(constants.MOL_ADJACENCY_PATH)
        for file in ligand_data:
            ligand_adjacency_matrix = np.load(os.path.join(constants.MOL_ADJACENCY_PATH, file))
        return ligand_adjacency_matrix
    

    def createModel(self):
        input_layer = tf.keras.layers.Input()
        conv_layer1 = tf.keras.layers.Conv2D(
            filters=2, kernel_size=3, activation='relu', padding='same'
        )(input_layer)
        pool_layer1 = tf.keras.layers.AveragePooling2D(
            pool_size=(2,2), strides=(1,1), padding='valid'
        )(conv_layer1)
        flat_layer = tf.keras.layers.Flatten()(pool_layer1)

        adjacency_model = tf.keras.models.Model(
            input_layer,
            conv_layer1,
            pool_layer1,
            flat_layer,
            name='Ligand-Adjacency-Model'
        )

        return adjacency_model

