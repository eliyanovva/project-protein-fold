import tensorflow as tf
import numpy as np
import logging as log
import os
import abc

import config

from graph_cnn.data_prep.data_handlers import DataHandlers


class ProteinAdjacencyData(DataHandlers):
    
    def loadDataSingleMatrix(self, label_name):
        # the protein name should be from the train/test X data.
        protein_adjacency_matrix = np.load(
            os.path.join(self.folder, label_name + '_adj_mat.npy')
        )
        return protein_adjacency_matrix


class ProteinFeatureData(DataHandlers):

    def loadDataSingleMatrix(self, label_name):
        # the protein name should be from the train/test X data.
        protein_feature_matrix = np.load(
            os.path.join(
                self.folder,
                label_name + '_feat_mat.npy')
        )
        return protein_feature_matrix

