import sys
sys.path.append('../')
import os
import config
import subprocess
from graph_cnn.data_prep.bgf_file_prep import BGFDataFile
from graph_cnn.data_prep.mol_file_prep import MolDataFile

def createBGFfile(pdb_protein_path, bgf_directory):
    protein_name = pdb_protein_path[max(pdb_protein_path.rfind('/') + 1, 0) : pdb_protein_path.rfind('.')]
    subprocess.run(
        ["obabel " + pdb_protein_path + " -O " + os.path.join(bgf_directory, protein_name  + ".bgf")],
        shell=True
    )

    return os.path.join(bgf_directory, protein_name + ".bgf")

def generateProteinMatrices(
    pdb_path=config.PDB_FILES_PATH,
    bgf_path=config.BGF_FILES_PATH,
    target_adj_path=config.PROTEIN_ADJACENCY_PATH,
    target_feat_path=config.PROTEIN_FEATURE_PATH
):
    for filename in os.listdir(pdb_path):
        if filename.endswith('.pdb'):
            bgf_file_path = createBGFfile(os.path.join(pdb_path, filename), bgf_path)
            protein_handler = BGFDataFile(bgf_file_path, os.path.join(pdb_path, filename))
            protein_handler.getAdjacencyMatrix(target_folder=target_adj_path)
            protein_handler.getFeatureMatrix(target_folder=target_feat_path)


def generateLigandMatrices(
    mol_path=config.MOL_FILES_PATH,
    target_adj_path=config.MOL_ADJACENCY_PATH,
    target_feat_path=config.LIGAND_FEATURE_PATH
):
    for filename in os.listdir(mol_path):
        if filename.endswith('.mol'):
            mol_file = MolDataFile(os.path.join(mol_path, filename))
            mol_file.getAdjacencyMatrix(target_folder=target_adj_path)
            mol_file.getFeatureMatrix(target_folder=target_feat_path)

