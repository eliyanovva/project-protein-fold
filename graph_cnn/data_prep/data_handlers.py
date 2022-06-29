import numpy as np
import pandas as pd
import logging as log
import os
import re

import constants

class DataHandlers:
    def __init__(self):
        pass

    def createDatasetCSV(self):
        with open('uniprot_ligand_logfc_pvalue.csv') as couples_csv:
            couples_list = couples_csv.readlines()
            cur_dataset_list = couples_list
            with open ('dataset.csv', 'a') as dataset:
                dataset.write('ProteinAdjacency, ProteinFeatures, LigandAdjacency, logFC')
                for i in range(len(cur_dataset_list)):
                    line_list = cur_dataset_list[i].split(',')
                    new_list = [0, 1, 2, 3]
                    new_list[0], new_list[1] = self.__fetchProteinData(line_list[0])
                    log.info('first add successfull for row ', i, ' out of ', len(cur_dataset_list))
                    new_list[2] = self.__fetchLigandData(line_list[1])
                    log.info('second add successfull for row ', i, ' out of ', len(cur_dataset_list))
                    new_list[3] = (line_list[2])
                    log.info('third add successfull for row ', i, ' out of ', len(cur_dataset_list))

                    dataset.write(','.join(new_list))

    def concatenateDatasets(self):
        data_prep_files = os.listdir('.')
        for file in data_prep_files:
            if file.endswith('.pkl'):
                pass

# TODO: fix the fact that this function exists in 2 other files :)
    def __fetchProteinData(self, protein_name):
        # the protein name should be from the train/test X data.
        protein_adjacency_matrix = np.load(
            os.path.join(constants.PROTEIN_ADJACENCY_PATH, protein_name + '_adj_mat.npy')
        )
        protein_feature_matrix = np.load(
            os.path.join(constants.PROTEIN_FEATURE_PATH, protein_name + '_feat_mat.npy')
        )
        return np.array2string(protein_adjacency_matrix, precision=4, separator=',',
                      suppress_small=True), np.array2string(protein_feature_matrix, precision=4, separator=',',
                      suppress_small=True)


    def __fetchLigandData(self, ligand_name):
        # the protein name should be from the train/test X data.
        adjacency_file_name = self.__getLigandFileNames(ligand_name)
        ligand_adjacency_matrix = np.load(
            os.path.join(constants.MOL_ADJACENCY_PATH, adjacency_file_name)
        )
        return np.array2string(ligand_adjacency_matrix, precision=4, separator=',',
                      suppress_small=True)
    
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


dh = DataHandlers()
dh.createDatasetCSV()