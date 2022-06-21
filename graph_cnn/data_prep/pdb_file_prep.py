import numpy as np
import math
import constants
import logging as log

# so far creates an adjacency matrix from a single pdb file
class PDB_Datafile:
    def __init__(self, pdb_filename):
        log.info('Class initialized!')
        self.pdb_filename = pdb_filename
        self.atom_count = 0
        self.protein_name = ''


    def getAdjacencyMatrix(self):
    # the adjacency matrix is created with the size of the protein, and will later go through a
    # convolution until it reaches a certain size
        log.info('Initiated creation of PDB Adjacency matrix for protein ' + self.protein_name)
        pdb_data = self.__getPDBFileAtomData()
        adjacency_matrix = np.zeros((self.atom_count, self.atom_count))
        for i in range(self.atom_count):
            for j in range(self.atom_count):
                adjacency_matrix[i][j] = self.__distanceBetweenTwoAtoms(
                    pdb_data[i][1:4], pdb_data[j][1:4]
                )
        log.info('Initiated saving of adjacency matrix')
        np.save('pdb_adjacency_data/' + self.protein_name + '_adj_mat', adjacency_matrix)
        log.info('The adjacency matrix has been saved!')

    def getFeatureMatrix(self):
        # current number of features - 5: atom type, coordinates, confidence score
        log.info('Initiated creation of PDB Feature matrix for protein ' + self.protein_name)
        feature_matrix = np.zeros((self.atom_count, 5))
        pdb_features = self.__getFeaturesFromPDB()
        for i in range(self.atom_count):
            for j in range(5):
                feature_matrix[i][j] = pdb_features[i][j]
        log.info('Initiated saving of feature matrix')
        np.save('pdb_features_data/' + self.protein_name + '_feat_mat', feature_matrix)
        log.info('The features matrix has been saved!')

    def __getFeaturesFromPDB(self):
        pdb_data = self.__getPDBFileAtomData()
        data_stripped_rows = np.delete(pdb_data, np.s_[self.atom_count:], 0)
        data_stripped_cols = np.delete(data_stripped_rows, np.s_[self.atom_count:], 1)
        return data_stripped_cols        

    def __setProteinName(self):
        left_index = self.pdb_filename.rfind('/') + 1
        right_index = self.pdb_filename.find('-model')
        self.protein_name = self.pdb_filename[left_index : right_index]


    def __distanceBetweenTwoAtoms(self, coordinates_A, coordinates_B):
        distance = math.sqrt(
            (coordinates_A[0] - coordinates_B[0])**2 + \
            (coordinates_A[1] - coordinates_B[1])**2 + \
            (coordinates_A[2] - coordinates_B[2])**2
        )
        return distance


    def __getDataFromFileLine(self, data_line):
        data_line = data_line.strip()
        numerical_data = []
        if data_line[0:4] == 'ATOM' and data_line[-1] != 'H':
            numerical_data = data_line[27:67].split(' ')
            numerical_data = [float(x) for x in numerical_data if x]
        return numerical_data


    def __getPDBFileAtomData(self):
        log.info('Started data extraction from ', self.pdb_filename)
        data = np.zeros((constants.ATOMS_COUNT, 5))
        self.__setProteinName()
        with open(self.pdb_filename) as file:
            log.info('PDB file has been succesfully opened')
            pdb_lines = file.readlines()
            data_index = 0
            for line in pdb_lines:
                numerical_data = self.__getDataFromFileLine(line)
                # TODO: ugly fix for the strip here, fix it later
                atom_type = line.strip()[-1]
                if len(numerical_data) > 0:
                    self.atom_count += 1
                    entry = np.array(
                        [
                            constants.ATOM_DICT[atom_type], # atom type (C, O, N, S)
                            numerical_data[0], # x coordinate
                            numerical_data[1], # y_coordinate
                            numerical_data[2], # z_coordinate
                            numerical_data[4]  # confidence score
                        ]
                    )
                    data[data_index] = entry
                    data_index += 1
        log.info('Data has been succesfully generated')
        return data





    # create feature matrix