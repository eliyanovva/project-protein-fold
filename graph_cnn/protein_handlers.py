import tensorflow as tf
import numpy as np
import logging as log
import os
import abc

import data_prep.constants as constants
from data_prep.data_handlers import DataHandlers


class ProteinAdjacencyData(DataHandlers):
    
    #@abc.abstractmethod
    def loadDataSingleMatrix(self, label_name):
        # the protein name should be from the train/test X data.
        protein_adjacency_matrix = np.load(
            os.path.join(constants.PROTEIN_ADJACENCY_PATH, label_name + '_adj_mat.npy')
        )
        return protein_adjacency_matrix


class ProteinFeatureData(DataHandlers):

    def loadDataSingleMatrix(self, label_name):
        # the protein name should be from the train/test X data.
        protein_feature_matrix = np.load(
            os.path.join(constants.PROTEIN_FEATURE_PATH, label_name + '_feat_mat.npy')
        )
        return protein_feature_matrix


class ProteinFeatureDataPDB:
    def __init__(self):
        pass

    
    def __fetchData(self):
        protein_data = os.listdir(constants.PROTEIN_FEATURE_PATH_PDB)
        for file in protein_data:
            protein_feature_matrix = np.load(os.path.join(constants.PROTEIN_FEATURE_PATH_PDB, file))
        return protein_feature_matrix