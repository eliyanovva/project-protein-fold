import logging as log
import os

import numpy as np
from openbabel import pybel

import config


#adjacency matrix from bgf file, feature matrix from pdb and bgf file
class BGFDataFile:
    def __init__(self, bgf_filename, pdb_filename):
        self.bgf_filename = bgf_filename
        self.pdb_filename = pdb_filename
        self.__setProteinName()
        self.__getSizes()
        log.info('Class initialized!')
        

    def getAdjacencyMatrix(self, target_folder=config.PROTEIN_ADJACENCY_PATH):
        log.info('Initiated creation of BGF Adjacency matrix for protein ' + self.protein_name)
        adjacency_matrix = np.zeros((self.atom_count, self.atom_count))
        with open(self.bgf_filename, 'r') as bgf_file:
            bgf_lines = bgf_file.readlines()[5:]
            for i in range(self.atom_count + 2, self.atom_count + self.bond_count + 2, 2):
                data_lines = bgf_lines[i:i+2]
                data = self.__getDataFromConectOrderLines(data_lines)
                for bond in data:
                    adjacency_matrix[int(bond[0])][int(bond[1])] = bond[2]    
        log.info('Initiated saving of adjacency matrix')
        np.save(os.path.join(target_folder, self.protein_name + '_adj_mat'), adjacency_matrix)
        log.info('The adjacency matrix has been saved!')


    def getFeatureMatrix(self, target_folder=config.PROTEIN_FEATURE_PATH):
        log.info('Initiated creation of BGF Feature matrix for protein ' + self.protein_name)
        feature_matrix = np.zeros((self.atom_count, config.PROTEIN_FEATURES_COUNT), dtype='float')
        confidence_scores = self.__getConfidenceScores()
        atom_types = self.__getAtomTypes()
       
        prefile = next(pybel.readfile('pdb', self.pdb_filename))
        molecule = pybel.Molecule(prefile)

        atom_index = 0
        for atom in molecule:
            feature_matrix[atom_index][:-2] = self.__getMoleculeDataRow(atom)
            feature_matrix[atom_index][-2] = normalize(atom_types[atom_index], 4.0)
            feature_matrix[atom_index][-1] = normalize(confidence_scores[atom_index], 100.0)
            atom_index += 1

        log.info('Initiated saving of feature matrix')
        np.save(os.path.join(target_folder, self.protein_name + '_feat_mat'), feature_matrix)
        log.info('The features matrix has been saved!')


    def __getMoleculeDataRow(self, atom):
        features_arr = np.zeros((config.PROTEIN_FEATURES_COUNT - 2))
        features_arr[0] = normalize(atom.atomicmass, 32.065)
        features_arr[1] = normalize(atom.exactmass, 31.972071)
        features_arr[2] = normalize(atom.heavydegree, 4.0)
        features_arr[3] = normalize(atom.heterodegree, 3.0)
        features_arr[4] = normalize(atom.hyb, 3.0)
        # FIXME: create a classification index which is based on where in the protein can we find the atom
        features_arr[5] = normalize(atom.idx, 3500.0)
        features_arr[6] = normalize(atom.partialcharge, 0.6, -0.6)
        features_arr[7] = normalize(atom.degree, 4.0)
        return features_arr


    def __setProteinName(self):
        left_index = max(self.bgf_filename.rfind('/'), self.bgf_filename.rfind('\\'))
        right_index = self.bgf_filename.rfind('.')
        self.protein_name = self.bgf_filename[max(0, left_index + 1) : right_index]  
        

    def __getSizes(self):
        """This function generates the count of atoms and the count of bonds in the bgf file.
        It refers to the bgf file format formulated in 
        http://www.chem.cmu.edu/courses/09-560/docs/msi/modenv/D_Files.html#944609
        """
        with open(self.bgf_filename, 'r') as bgf_file:
            lines = bgf_file.readlines()
            self.atom_count = 0
            self.bond_count = 0
            if len(lines) == 0:
                raise Exception('The bgf file ' + self.bgf_filename + ' is empty')
            for line in lines:
                if line.startswith('HETATM'):
                    self.atom_count += 1
                elif line.startswith('CONECT'):
                    self.bond_count += 1


    def __getDataFromConectOrderLines(self, data_lines):
        """
        Collects bond data from CONECT-ORDER paor of lines in a bgf file.

        Args:
            data_lines (List[str]): a CONECT line and a ORDER line.

        Returns:
            List: A list of lists in this format: [bonding atom 1, bonding atom 2, bond val]
        """
        conect_line = data_lines[0].split()[1:]
        order_line = data_lines[1].split()[1:]
        base_atom_index = conect_line[0]
        data_line = []
        for i in range(1, len(conect_line)):
            data_line.append([base_atom_index, conect_line[i], order_line[i]])
        return data_line


    def __getConfidenceScores(self):
        log.info('Initiated extracting of confidence scores from PDB file for ' + self.protein_name)
        confidence_scores = np.zeros(self.atom_count, dtype='float')
        with open(self.pdb_filename) as pdb_file:
            lines = pdb_file.readlines()
            index = 0
            for line in lines:
                if line.startswith('ATOM'):
                    line_list = [x for x in line.split() if x != '']
                    confidence_scores[index] = float(line_list[-2])
                    index += 1
        log.info('Extraction of confidence scores from PDB completed!')
        return confidence_scores


    def __getAtomTypes(self):
        log.info('Initiated extracting of atom types from PDB file for ' + self.protein_name)
        atom_types = np.zeros(self.atom_count, dtype='float')
        with open(self.pdb_filename) as pdb_file:
            lines = pdb_file.readlines()
            index = 0
            for line in lines:
                if line.startswith('ATOM'):
                    line_list = [x for x in line.split() if x != '']
                    #FIXME: add check if it is a valid key; assign non-valid to a key for others
                    atom_types[index] = config.ATOM_DICT[line_list[-1]]
                    index += 1
        log.info('Extraction of atom types from PDB completed!')
        return atom_types


def normalize(x, max, min=0.0):
    y = (float(x) - float(min))/(float(max) - float(min))
    return y