import numpy as np
import math
import constants
import logging as log

# so far creates an adjacency matrix from a single pdb file
class PDB_Datafile:
    def __init__(self, pdb_filename):
        log.info('Class initialized!')
        self.pdb_filename = pdb_filename
        self.protein_name = ''


    def getAdjacencyMatrix(self):
    # the adjacency matrix is created with the size of the constant set up in constants.py
        log.info('Initiated creation of PDB Adjacency matrix for protein ' + self.protein_name)
        adjacency_matrix = np.zeros((constants.ATOMS_COUNT, constants.ATOMS_COUNT))
        pdb_data = self.__getPDBFileAtomData()
        for i in range(constants.ATOMS_COUNT):
            for j in range(constants.ATOMS_COUNT):
                adjacency_matrix[i][j] = self.__distanceBetweenTwoAtoms(
                    pdb_data[i][1:4], pdb_data[j][1:4]
                )
        log.info('Initiated saving of adjacency matrix')
        np.save('pdb_adjacency_data/' + self.protein_name + '_adj_mat', adjacency_matrix)
        log.info('The adjacency matrix has been saved!')


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