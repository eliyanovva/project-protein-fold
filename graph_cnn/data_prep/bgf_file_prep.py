import numpy as np
import math
import constants
import logging as log
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
        left_index = self.bgf_filename.rfind('/') + 1
        right_index = self.bgf_filename.find('-model')
        self.protein_name = self.bgf_filename[left_index : right_index]


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


    def __distanceBetweenTwoAtoms(self, coordinates_A, coordinates_B):
        distance = math.sqrt(
            (coordinates_A[0] - coordinates_B[0])**2 + \
            (coordinates_A[1] - coordinates_B[1])**2 + \
            (coordinates_A[2] - coordinates_B[2])**2
        )
        return distance


    def __getDataFromHetAtmLine(self, data_line):
        data_line = data_line.strip().split()
        data_line = [
            float(data_line[1]), # atom count
            constants.ATOM_DICT[data_line[2][0]], # atom type
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


log.info('log test')
m_class = BGFDataFile('/home/users/tep18/new_ppp/project-protein-fold/bgf_files/AF-Q0VAX9-F1-model_v2.bgf')
m_class.getAdjacencyMatrix()

#a = np.load('/home/users/tep18/new_ppp/project-protein-fold/graph_cnn/data_prep/mol_adjacency_data/C11H22O_adj_mat.npy')
#print(a)
    # create feature matrix


    # create feature matrix