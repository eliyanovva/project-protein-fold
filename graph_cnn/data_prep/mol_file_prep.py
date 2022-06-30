import numpy as np
import math
import os
import constants
import logging as log
import log_config
from feature_matrices import getMatrix

# so far creates an adjacency matrix from a single pdb file
# TODO: Figure out how/whether you need to save adjacencies to empty
# places as 0s or as current distances or as -1

class MolDataFile:
    def __init__(self, mol_filename):
        log.info('Class initialized!')
        self.mol_filename = mol_filename
        self.__setCompoundName()
        self.__getSizes()
        

    def getAdjacencyMatrix(self):
    # the adjacency matrix is created with the size of the constant set up in constants.py
        log.info('Initiated creation of Mol Adjacency matrix for the compound ' + self.compound_name)
        adjacency_matrix = np.zeros((self.atom_count, self.atom_count))
        bond_data = self.__getMolFileBondData()
        for i in range(self.bond_count):
            adjacency_matrix[bond_data[i][0] - 1][bond_data[i][1] - 1] = bond_data[i][2]
        
        log.info('Initiated saving of adjacency matrix')
        file_name = os.path.join('mol_adjacency_data/', self.compound_name + '_adj_mat')
        np.save(file_name, adjacency_matrix)
        log.info('The adjacency matrix has been saved!')

    def getFeatureMatrix(self):
    #building feature matrix for mol file
        log.info('Initiated creation of Mol Feature matrix for the compound ' + self.compound_name)
        getMatrix('mol', self.mol_filename)


    def __setCompoundName(self):
        right_index = self.mol_filename.rfind('/')
        left_index = self.mol_filename.find('.')
        self.compound_name = self.mol_filename[right_index + 1 : left_index]


    def __distanceBetweenTwoAtoms(self, coordinates_A, coordinates_B):
        distance = math.sqrt(
            (coordinates_A[0] - coordinates_B[0])**2 + \
            (coordinates_A[1] - coordinates_B[1])**2 + \
            (coordinates_A[2] - coordinates_B[2])**2
        )
        return distance


    def __get3DDataFromFileLine(self, data_line):
        # TODO might have to change this function to extract feature data from the mol file.
        data_line = data_line.strip()
        data_list = [i for i in data_line.split() if i]
        numerical_data = []
        if len(data_list) > 4 and data_list[3] not in "0123456789":
            numerical_data = data_list[0:3]
            numerical_data = [float(x) for x in numerical_data if x]
        return numerical_data


    def __getBondDataFromFileLine(self, data_line):
        data_line = data_line.strip()
        data_list = [i for i in data_line.split() if i]
        numerical_data = data_list[0:3]
        numerical_data = [int(x) for x in numerical_data if x]
        return numerical_data


    def __getMolFile3DData(self):
        log.info('Started 3d coordinates data extraction from ', self.mol_filename)
        data = np.zeros((self.atom_count, 3))
        with open(self.mol_filename) as file:
            log.info('Mol file has been succesfully opened')
            mol_lines = file.readlines()[4:]
            data_index = 0
            for line in mol_lines:
                numerical_data = self.__get3DDataFromFileLine(line)
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


    def __getSizes(self):
        """This function generates the count of atoms and the count of bonds in the mol file.
        It refers to the mol file format formulated in 
        https://chem.libretexts.org/Courses/University_of_Arkansas_Little_Rock/ChemInformatics_(2017)%3A_Chem_4399_5399/2.2%3A_Chemical_Representations_on_Computer%3A_Part_II/2.2.2%3A_Anatomy_of_a_MOL_file
        """
        with open(self.mol_filename, 'r') as mol_file:
            lines = mol_file.readlines()
            if len(lines) == 0:
                raise Exception('The mol file ' + self.compound_name + ' is empty')
            size_line = lines[3].split()
            self.atom_count = int(size_line[0])
            self.bond_count = int(size_line[1])


    def __getMolFileBondData(self):
        log.info('Started bond data extraction from ' + self.mol_filename)
        data = np.zeros((self.bond_count, 3), dtype='int32')
        with open(self.mol_filename) as file:
            log.info('Mol file has been succesfully opened')
            mol_lines = file.readlines()[4 + self.atom_count:4 + self.atom_count + self.bond_count]
            data_index = 0
            for line in mol_lines:
                numerical_data = self.__getBondDataFromFileLine(line)
                if len(numerical_data) > 0:
                    entry = np.array(
                        numerical_data,
                        dtype='int32'
                    )
                    data[data_index] = entry
                    data_index += 1
        log.info('Bond data has been succesfully generated')
        return data




#m_class = MolDataFile('/home/users/tep18/new_ppp/project-protein-fold/moleculekit/mol_files/Undecanal.mol')
#m_class.getAdjacencyMatrix()

#a = np.load('/home/users/tep18/new_ppp/project-protein-fold/graph_cnn/data_prep/mol_adjacency_data/C11H22O_adj_mat.npy')
#print(a)
    # create feature matrix