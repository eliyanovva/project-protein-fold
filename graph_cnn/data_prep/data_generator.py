import os
from pdb_file_prep import PDB_Datafile
from mol_file_prep import MolDataFile
def generatePDBAdjacencyMatrices(folder_path):
    for filename in os.listdir(folder_path):
        # TODO: add regex check to ensure it is a pdb file
        if filename.endswith('.pdb'):
            pdb_file = PDB_Datafile(folder_path + filename)
            pdb_file.getAdjacencyMatrix()


def generateMolAdjacencyMatrices(folder_path):
    for filename in os.listdir(folder_path):
        if filename.endswith('.mol'):
            mol_file = MolDataFile(folder_path + filename)
            mol_file.getAdjacencyMatrix()


generateMolAdjacencyMatrices('/home/users/tep18/new_ppp/project-protein-fold/mol_files/')
#generatePDBAdjacencyMatrices('/home/users/tep18/new_ppp/project-protein-fold/pdb_data_files/')