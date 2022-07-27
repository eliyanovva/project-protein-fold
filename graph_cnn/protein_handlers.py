import numpy as np
import logging as log
import os
import re

import config

from graph_cnn.data_prep.data_handlers import DataHandlers


class ProteinAdjacencyData(DataHandlers):
    
    def loadDataSingleMatrix(self, label_name):
        # the protein name should be from the train/test X data.
        mylist = os.listdir(self.folder)
        r = re.compile(".*"+ label_name + ".*npy")
        adjacency_matrix_filename = list(filter(r.match, mylist))[0]
        
        protein_adjacency_matrix = np.load(
            os.path.join(self.folder, adjacency_matrix_filename)
        )

        return protein_adjacency_matrix


class ProteinFeatureData(DataHandlers):

    def loadDataSingleMatrix(self, label_name):
        # the protein name should be from the train/test X data.
        mylist = os.listdir(self.folder)
        r = re.compile(".*"+ label_name + ".*npy")
        feature_matrix_filename = list(filter(r.match, mylist))[0]
        
        protein_feature_matrix = np.load(
            os.path.join(
                self.folder,
                feature_matrix_filename)
        )

        return protein_feature_matrix

