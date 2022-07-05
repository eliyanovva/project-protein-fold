import numpy as np
import pandas as pd
import logging as log
import os
import sys
import re
from . import log_config
from . import constants

class DataHandlers:
    def __init__(self):
        pass

    def createDatasetCSV(self):
        with open('/home/users/tep18/new_ppp/project-protein-fold/graph_cnn/data_prep/uniprot_ligand_logfc_pvalue.csv', 'r') as couples_csv:
            couples_list = couples_csv.readlines()
            cur_dataset_list = couples_list
            with open ('/home/users/tep18/new_ppp/project-protein-fold/graph_cnn/data_prep/dataset.csv', 'a') as dataset:
                dataset.truncate(0)
                #dataset.write('ProteinAdjacency, ProteinFeatures, LigandAdjacency, logFC\n')
                for i in range(3):#len(cur_dataset_list)):
                    line_list = cur_dataset_list[i].split(',')
                    
                    new_list = [0, 1, 2, 3]

                    protein_data = self.fetchProteinData(line_list[0])
                    ligand_data = self.fetchLigandData(line_list[1])
                    new_list[0], new_list[1] = self.__convertNumpyToString(protein_data[0]), self.__convertNumpyToString(protein_data[1])
                    log.info('first add successfull for row ' + str(i) + ' out of ' + str(len(cur_dataset_list)))
                    new_list[2] = self.__convertNumpyToString(self.__fetchLigandData(line_list[1]))
                    log.info('second add successfull for row ' + str(i) + ' out of ' + str(len(cur_dataset_list)))
                    new_list[3] = (line_list[2])
                    log.info('third add successfull for row ' + str(i) + ' out of ' + str(len(cur_dataset_list)))

                    
                    #print(','.join(new_list))
                    #dataset.write('\t'.join(new_list) + '\n')

    
    def createDatasetPickles(self):
        with open('/home/users/tep18/new_ppp/project-protein-fold/graph_cnn/data_prep/uniprot_ligand_logfc_pvalue.csv', 'r') as couples_csv:
            couples_list = couples_csv.readlines()
            cur_dataset_list = couples_list
            row_count = len(couples_list)
            dataset_count = row_count // 250 + 1
            print(dataset_count)

            for i in range(dataset_count):
                dataset_size = min(250, row_count - i * 250)
                df = pd.DataFrame(index=range(dataset_size),columns=range(4))
                
                for j in range(i*250, i*250 + dataset_size):
                    line_list = cur_dataset_list[j].split(',')
                    protein_data = self.fetchProteinData(line_list[0])
                    ligand_data = self.fetchLigandData(line_list[1])
                    logFc = line_list[2]

                    log.info('data extraction completed' + str(j))

                    df.iat[j % 250, 0] = protein_data[0].resize((constants.PROTEIN_ADJACENCY_MAT_SIZE, constants.PROTEIN_ADJACENCY_MAT_SIZE))
                    
                    log.info('protein adjacency matrix resized and added to df')
                    df.iat[j % 250, 1] = protein_data[1].resize((constants.PROTEIN_ADJACENCY_MAT_SIZE, constants.PROTEIN_FEATURES_COUNT))
                    log.info('protein feature matrix resized and added to df')
                    df.iat[j % 250, 2] = ligand_data.resize((constants.LIGAND_ADJACENCY_MAT_SIZE, constants.LIGAND_ADJACENCY_MAT_SIZE))
                    log.info('ligand adjacency matrix resized and added to df')
                    df.iat[j % 250, 3] = logFc
                    log.info('dataframe load complete' +  str(j))

                df.to_pickle(os.path.join(constants.MATRIX_DATA_FILES_PATH, 'dataset_part_' + str(i) + '.pkl'))


    def concatenateDatasets(self):
        data_prep_files = os.listdir('.')
        for file in data_prep_files:
            if file.endswith('.pkl'):
                pass

# TODO: fix the fact that this function exists in 2 other files :)
    def fetchProteinData(self, protein_name):
        # the protein name should be from the train/test X data.
        protein_adjacency_matrix = np.load(
            os.path.join(constants.PROTEIN_ADJACENCY_PATH, protein_name + '_adj_mat.npy')
        )
        protein_feature_matrix = np.load(
            os.path.join(constants.PROTEIN_FEATURE_PATH, protein_name + '_feat_mat.npy')
        )
        return protein_adjacency_matrix, protein_feature_matrix


    def fetchLigandData(self, ligand_name):
        # the protein name should be from the train/test X data.
        adjacency_file_name = self.__getLigandFileNames(ligand_name)
        ligand_adjacency_matrix = np.load(
            os.path.join(constants.MOL_ADJACENCY_PATH, adjacency_file_name)
        )
        return ligand_adjacency_matrix
    

    def __convertNumpyToString(self, np_matrix):
        log.info('started converting numpy file to string')
        arr_str = "["
        for row in np_matrix:
            row_str = '['
            for item in row:
                row_str = row_str + str(item) + ","
            row_str = row_str[:-1]
            row_str = row_str + '],'
            arr_str = arr_str + row_str
        arr_str = arr_str[:-1]
        arr_str = arr_str + "]"
        log.info('finished numpy file conversion')
        #print(arr_str)
        return arr_str


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


#dh = DataHandlers()
#dh.createDatasetPickles()