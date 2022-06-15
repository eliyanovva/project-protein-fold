import os
from pdb_file_prep import PDB_Datafile

def generatePDBAdjacencyMatrices(folder_path):
    for filename in os.listdir(folder_path):
        # TODO: add regex check to ensure it is a pdb file
        if filename.endswith('.pdb'):
            pdb_file = PDB_Datafile(folder_path + filename)
            pdb_file.getAdjacencyMatrix()

generatePDBAdjacencyMatrices('/home/users/tep18/new_ppp/project-protein-fold/pdb_data_files/')