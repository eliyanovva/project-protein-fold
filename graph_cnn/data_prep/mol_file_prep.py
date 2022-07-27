import sys
import logging as log
import math
import os
CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(os.path.dirname(CURRENT_DIR)))

import numpy as np
from openbabel import pybel

import config

class MolDataFile:
    def __init__(self, mol_filename):
        log.info('Class initialized!')
        self.mol_filename = mol_filename
        self.__setCompoundName()
        self.__getSizes()
        

    def getAdjacencyMatrix(self, target_folder=config.MOL_ADJACENCY_PATH):
    # the adjacency matrix is created with the size of the constant set up in config.py
        log.info('Initiated creation of Mol Adjacency matrix for the compound ' + self.compound_name)
        adjacency_matrix = np.zeros((self.atom_count, self.atom_count))
        bond_data = self.__getMolFileBondData()
        for i in range(self.bond_count):
            adjacency_matrix[bond_data[i][0] - 1][bond_data[i][1] - 1] = bond_data[i][2]
        
        log.info('Initiated saving of adjacency matrix')
        file_name = os.path.join(target_folder, self.compound_name + '_adj_mat')
        np.save(file_name, adjacency_matrix)
        log.info('The adjacency matrix has been saved!')


    def getFeatureMatrix(self, target_folder=config.LIGAND_FEATURE_PATH):
    #building feature matrix for mol file
        log.info('Initiated creation of Mol Feature matrix for the compound ' + self.compound_name)
        feature_matrix = np.zeros((self.atom_count, config.LIGAND_FEATURES_COUNT), dtype='float')
        atom_types = self.__getAtomTypes()
       
        prefile = next(pybel.readfile('mol', self.mol_filename))
        molecule = pybel.Molecule(prefile)

        atom_index = 0
        for atom in molecule:
            feature_matrix[atom_index][:-1] = self.__getMoleculeDataRow(atom)
            feature_matrix[atom_index][-1] = normalize(atom_types[atom_index], 4.0)
            atom_index += 1

        log.info('Initiated saving of feature matrix')
        np.save(os.path.join(target_folder, self.compound_name + '_feat_mat'), feature_matrix)
        log.info('The features matrix has been saved!')


    def __setCompoundName(self):
        right_index = max(self.mol_filename.rfind('/'), self.mol_filename.rfind('\\'))
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


    def __getMoleculeDataRow(self, atom):
        features_arr = np.zeros((config.LIGAND_FEATURES_COUNT - 1))
        features_arr[0] = normalize(atom.atomicmass, 32.065)
        features_arr[1] = normalize(atom.exactmass, 31.972071)
        #features_arr[2] = atom.formalcharge
        features_arr[2] = normalize(atom.heavydegree, 4.0)
        features_arr[3] = normalize(atom.heterodegree, 3.0)
        features_arr[4] = normalize(atom.hyb, 3.0)
        # FIXME: create a classification index which is based on where in the protein can we find the atom
        features_arr[5] = normalize(atom.idx, 50.0)
        #features_arr[7] = atom.isotope
        features_arr[6] = normalize(atom.partialcharge, 0.6, -0.6)
        #features_arr[9] = atom.spin
        #FIXME: EXTRACT TYPE FROM ELSEWHERE
        #features_arr[10] = x.type
        features_arr[7] = normalize(atom.degree, 4.0)
        return features_arr


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
                            #config.ATOM_DICT[atom_type], # atom type (C, O, N, S)
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


    def __getAtomTypes(self):
        log.info('Initiated extracting of atom types from MOL file for ' + self.compound_name)
        atom_types = np.zeros(self.atom_count, dtype='float')
        with open(self.mol_filename) as mol_file:
            lines = mol_file.readlines()
            index = 0
            for line in lines[4:4 + self.atom_count]:
                line_list = [x for x in line.split() if x != '']
                #FIXME: add check if it is a valid key; assign non-valid to a key for others
                atom_types[index] = config.ATOM_DICT[line_list[3]]
                index += 1
        log.info('Extraction of atom types from MOL completed!')
        return atom_types


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
    

def normalize(x, max, min=0.0):
    y = (float(x) - float(min))/(float(max) - float(min))
    return y




#m_class = MolDataFile('/home/users/tep18/new_ppp/project-protein-fold/moleculekit/mol_files/Undecanal.mol')
#m_class.getFeatureMatrix()

#a = np.load('/home/users/tep18/new_ppp/project-protein-fold/graph_cnn/data_prep/mol_adjacency_data/C11H22O_adj_mat.npy')
#print(a)
    # create feature matrix
