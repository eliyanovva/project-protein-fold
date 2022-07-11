import os
import config
from bgf_file_prep import BGFDataFile
from mol_file_prep import MolDataFile


def generateProteinAdjacencyMatrices(folder_path=config.BGF_FILES_PATH):
    for filename in os.listdir(folder_path):
        if filename.endswith('.bgf'):
            pdb_file = BGFDataFile(os.path.join(folder_path, filename))
            pdb_file.getAdjacencyMatrix()


def generateProteinFeatureMatrices(folder_path=config.BGF_FILES_PATH):
    for filename in os.listdir(folder_path):
        if filename.endswith('.bgf'):
            bgf_file = BGFDataFile(os.path.join(folder_path, filename))
            bgf_file.getFeatureMatrix()


def generateLigandAdjacencyMatrices(folder_path=config.MOL_FILES_PATH):
    for filename in os.listdir(folder_path):
        if filename.endswith('.mol'):
            mol_file = MolDataFile(os.path.join(folder_path, filename))
            mol_file.getAdjacencyMatrix()


def generateLigandFeatureMatrices(folder_path=config.MOL_FILES_PATH):
    for filename in os.listdir(folder_path):
        if filename.endswith('.mol'):
            mol_file = MolDataFile(os.path.join(folder_path, filename))
            mol_file.getFeatureMatrix()


#generateLigandAdjacencyMatrices()
#generateLigandFeatureMatrices()
#generateProteinAdjacencyMatrices()
#generateProteinFeatureMatrices()