import numpy as np
import math
import constants
import logging as log
import log_config

# so far creates an adjacency matrix from a single pdb file
# TODO: Figure out how/whether you need to save adjacencies to empty
# places as 0s or as current distances or as -1

class MolDataFile:
    def __init__(self, mol_filename):
        log.info('Class initialized!')
        self.mol_filename = mol_filename
        self.compound_name = ''


    def getAdjacencyMatrix(self):
    # the adjacency matrix is created with the size of the constant set up in constants.py
        log.info('Initiated creation of Mol Adjacency matrix for the compound ' + self.compound_name)
        adjacency_matrix = np.zeros((constants.MOL_ATOMS_COUNT, constants.MOL_ATOMS_COUNT))
        mol_data = self.__getMolFileAtomData()
        for i in range(constants.MOL_ATOMS_COUNT):
            for j in range(constants.MOL_ATOMS_COUNT):
                adjacency_matrix[i][j] = self.__distanceBetweenTwoAtoms(
                    mol_data[i][0:3], mol_data[j][0:3]
                )
        log.info('Initiated saving of adjacency matrix')
        file_name = 'mol_adjacency_data/' + self.compound_name + '_adj_mat'
        np.save(file_name, adjacency_matrix)
        log.info('The adjacency matrix has been saved!')


    def __setCompoundName(self):
        with open(self.mol_filename, 'r') as file:
            self.compound_name = file.readline()[:-1]


    def __distanceBetweenTwoAtoms(self, coordinates_A, coordinates_B):
        distance = math.sqrt(
            (coordinates_A[0] - coordinates_B[0])**2 + \
            (coordinates_A[1] - coordinates_B[1])**2 + \
            (coordinates_A[2] - coordinates_B[2])**2
        )
        return distance


    def __getDataFromFileLine(self, data_line):
        # TODO might have to change this function to extract feature data from the mol file.
        data_line = data_line.strip()
        data_list = [i for i in data_line.split() if i]
        numerical_data = []
        if len(data_list) > 4 and data_list[3] not in "0123456789":
            numerical_data = data_list[0:3]
            numerical_data = [float(x) for x in numerical_data if x]
        return numerical_data


    def __getMolFileAtomData(self):
        log.info('Started data extraction from ', self.mol_filename)
        data = np.zeros((constants.MOL_ATOMS_COUNT, 3))
        self.__setCompoundName()
        with open(self.mol_filename) as file:
            log.info('Mol file has been succesfully opened')
            mol_lines = file.readlines()[4:]
            data_index = 0
            for line in mol_lines:
                numerical_data = self.__getDataFromFileLine(line)
                # TODO: ugly fix for the strip here, fix it later
                if len(numerical_data) > 0:
                    entry = np.array(
                        [
                            #constants.ATOM_DICT[atom_type], # atom type (C, O, N, S)
                            numerical_data[0], # x coordinate
                            numerical_data[1], # y_coordinate
                            numerical_data[2], # z_coordinate
                        #    numerical_data[4]  # confidence score
                        ]
                    )
                    data[data_index] = entry
                    data_index += 1
        log.info('Data has been succesfully generated')
        return data




m_class = MolDataFile('/home/users/tep18/new_ppp/project-protein-fold/moleculekit/mol_files/Undecanal.mol')
m_class.getAdjacencyMatrix()

a = np.load('/home/users/tep18/new_ppp/project-protein-fold/graph_cnn/data_prep/mol_adjacency_data/C11H22O_adj_mat.npy')
print(a)
    # create feature matrix