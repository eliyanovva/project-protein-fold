import logging as log
import os

import numpy as np
from openbabel import pybel

import constants
import log_config


# so far creates an adjacency matrix from a single pdb file
class BGFDataFile:
    def __init__(self, bgf_filename):
        self.bgf_filename = bgf_filename
        self.__setProteinName()
        self.__getSizes()
        log.info('Class initialized!')
        

    def getAdjacencyMatrix(self):
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
        np.save(os.path.join('pdb_adjacency_data/', self.protein_name + '_adj_mat'), adjacency_matrix)
        log.info('The adjacency matrix has been saved!')


    def getFeatureMatrix(self):
        # current number of features - 5: 
        # atom type, max covalent bonds, number of lone pairs, atomic charge, alphafold score
        log.info('Initiated creation of BGF Feature matrix for protein ' + self.protein_name)
        feature_matrix = np.zeros((self.atom_count, constants.PROTEIN_FEATURES_COUNT), dtype='float')
        confidence_scores = self.__getConfidenceScores()
        atom_types = self.__getAtomTypes()
       
        prefile = next(pybel.readfile('pdb', self.pdb_file_name))
        molecule = pybel.Molecule(prefile)

        atom_index = 0
        for atom in molecule:
            feature_matrix[atom_index][:-2] = self.__getMoleculeDataRow(atom)
            feature_matrix[atom_index][-2] = atom_types[atom_index]
            feature_matrix[atom_index][-1] = confidence_scores[atom_index]
            atom_index += 1

        log.info('Initiated saving of feature matrix')
        np.save(os.path.join(constants.PROTEIN_FEATURE_PATH, self.protein_name + '_feat_mat'), feature_matrix)
        log.info('The features matrix has been saved!')


    def __getMoleculeDataRow(self, atom):
        features_arr = np.zeros((constants.PROTEIN_FEATURES_COUNT - 2))
        features_arr[0] = atom.atomicmass
        features_arr[1] = atom.exactmass
        features_arr[2] = atom.formalcharge
        features_arr[3] = atom.heavydegree
        features_arr[4] = atom.heterodegree
        features_arr[5] = atom.hyb
        # FIXME: create a classification index which is based on where in the protein can we find the atom
        features_arr[6] = atom.idx
        features_arr[7] = atom.isotope
        features_arr[8] = atom.partialcharge
        features_arr[9] = atom.spin
        #FIXME: EXTRACT TYPE FROM ELSEWHERE
        #features_arr[10] = x.type
        features_arr[10] = atom.degree
        return features_arr


    def __setProteinName(self):
        left_index = self.bgf_filename.rfind('/AF-') + 4
        right_index = self.bgf_filename.find('-F1')
        self.protein_name = self.bgf_filename[left_index : right_index]
        self.pdb_file_name = os.path.join(
                constants.PDB_FILES_PATH,
                'AF-' + self.protein_name + '-F1-model_v2.pdb'
            )


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
                raise Exception('The bgf file ' + self.compound_name + ' is empty')
            for line in lines:
                if line.startswith('HETATM'):
                    self.atom_count += 1
                elif line.startswith('CONECT'):
                    self.bond_count += 1


    def __getDataFromHetAtmLine(self, data_line):
        data_line = data_line.strip().split()
        data_line = [
            int(data_line[1]), # atom count
            (constants.ATOM_DICT[data_line[2][0]]), # atom type
            int(data_line[-3]), # max number of covalent bonds
            int(data_line[-2]), # number of lone pairs
            float(data_line[-1])  # atomic charge
        ]
        return data_line


    def __getDataFromConectOrderLines(self, data_lines):
        """Collects bond data from CONECT-ORDER paor of lines in a bgf file.

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
        with open(self.pdb_file_name) as pdb_file:
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
        with open(self.pdb_file_name) as pdb_file:
            lines = pdb_file.readlines()
            index = 0
            for line in lines:
                if line.startswith('ATOM'):
                    line_list = [x for x in line.split() if x != '']
                    #FIXME: add check if it is a valid key; assign non-valid to a key for others
                    atom_types[index] = constants.ATOM_DICT[line_list[-1]]
                    index += 1
        log.info('Extraction of atom types from PDB completed!')
        return atom_types


m_class = BGFDataFile(os.path.join(constants.BGF_FILES_PATH, 'AF-Q0VAX9-F1-model_v2.bgf'))

m_class.getFeatureMatrix()

#a = np.load(os.join(constants.PROTEIN_FEATURE_PATH, )'/home/users/tep18/new_ppp/project-protein-fold/graph_cnn/data_prep/mol_adjacency_data/C11H22O_adj_mat.npy')
#print(a)
    # create feature matrix


    # create feature matrix
