import os
import constants
from bgf_file_prep import BGFDataFile
from mol_file_prep import MolDataFile


def generateBGFAdjacencyMatrices(folder_path=constants.BGF_FILES_PATH):
    for filename in os.listdir(folder_path):
        if filename.endswith('.bgf'):
            pdb_file = BGFDataFile(os.path.join(folder_path, filename))
            pdb_file.getAdjacencyMatrix()


def generateBGFFeatureMatrices(folder_path=constants.BGF_FILES_PATH):
    for filename in os.listdir(folder_path):
        if filename.endswith('.bgf'):
            bgf_file = BGFDataFile(os.path.join(folder_path, filename))
            bgf_file.getFeatureMatrix()


def generateMolAdjacencyMatrices(folder_path=constants.MOL_FILES_PATH):
    for filename in os.listdir(folder_path):
        if filename.endswith('.mol'):
            mol_file = MolDataFile(os.path.join(folder_path, filename))
            mol_file.getAdjacencyMatrix()


#generateMolAdjacencyMatrices('/home/users/tep18/new_ppp/project-protein-fold/mol_files/')
generateBGFAdjacencyMatrices()
generateBGFFeatureMatrices()