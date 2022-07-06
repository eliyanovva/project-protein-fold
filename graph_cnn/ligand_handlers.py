import tensorflow as tf
import numpy as np
import os
import re
import abc

import data_prep.constants as constants
from data_prep.data_handlers import DataHandlers


class LigandAdjacencyData(DataHandlers):

    #@abc.abstractmethod
    def loadDataSingleMatrix(self, label_name):
        # the protein name should be from the train/test X data.
        adjacency_file_name = self.__getLigandFileNames(label_name)
        ligand_adjacency_matrix = np.load(
            os.path.join(constants.MOL_ADJACENCY_PATH, adjacency_file_name)
        )
        return ligand_adjacency_matrix 
    
    
    def __getLigandFileNames(self, ligand_name):
        """The function returns the filename of the ligand adjacency matrix corresponding to the
        ligand name from the ligand column in the uniprot-ligand-logFC-pValue csv.

        Args:
            ligand_name (_type_): _description_
        """

        mylist = os.listdir(constants.MOL_ADJACENCY_PATH)
        left_index = ligand_name.rfind('_')
        r = re.compile(".*"+ ligand_name[left_index:] + "_adj_mat.npy")
        adjacency_matrix_filename = list(filter(r.match, mylist))
        return adjacency_matrix_filename[0]


class LigandFeatureData(DataHandlers):

    def loadDataSingleMatrix(self, label_name):
        # the protein name should be from the train/test X data.
        features_file_name = self.__getLigandFileNames(label_name)
        ligand_features_matrix = np.load(
            os.path.join(constants.LIGAND_FEATURE_PATH, features_file_name)
        )
        return ligand_features_matrix 
    
    
    def __getLigandFileNames(self, ligand_name):
        """The function returns the filename of the ligand adjacency matrix corresponding to the
        ligand name from the ligand column in the uniprot-ligand-logFC-pValue csv.

        Args:
            ligand_name (_type_): _description_
        """

        mylist = os.listdir(constants.LIGAND_FEATURE_PATH)
        print(mylist)
        left_index = ligand_name.rfind('_')
        print(ligand_name)
        r = re.compile(".*"+ ligand_name[left_index:] + ".*_feat_mat.npy")
        adjacency_matrix_filename = list(filter(r.match, mylist))
        print(adjacency_matrix_filename)
        return adjacency_matrix_filename[0]
