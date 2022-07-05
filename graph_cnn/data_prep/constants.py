import os

PROTEIN_ADJACENCY_MAT_SIZE = 3500
LIGAND_ADJACENCY_MAT_SIZE = 70
PROTEIN_FEATURES_COUNT = 5
ATOM_DICT = {'C':0, 'O':1, 'N':2, 'S':3}

#absolute path for a file
script_dir = os.path.abspath(os.curdir)
#final directory in filepath
last_dir = os.path.split(os.curdir)[-1]
#ensures the final directory is project-protein-fold
while last_dir != "project-protein-fold":
    os.chdir("..")
    script_dir = os.path.abspath(os.curdir)
    last_dir = os.path.split(script_dir)[-1]

#script_dir allows DATA_FILES_PATH to be called from any computer
DATA_FILES_PATH = os.path.join(script_dir, 'data_files')
PDB_FILES_PATH = os.path.join(DATA_FILES_PATH, 'pdb_data_files')
BGF_FILES_PATH = os.path.join(DATA_FILES_PATH, 'bgf_files')
MOL_FILES_PATH = os.path.join(DATA_FILES_PATH, 'mol_files')
SMILES_FILES_PATH = os.path.join(DATA_FILES_PATH, 'smiles_files')

MATRIX_DATA_FILES_PATH = os.path.join(script_dir, 'graph_cnn', 'data_prep')
MOL_ADJACENCY_PATH = os.path.join(MATRIX_DATA_FILES_PATH, 'mol_adjacency_data')
MOL_FEATURE_PATH = os.path.join(MATRIX_DATA_FILES_PATH, 'mol_features_data')
PROTEIN_ADJACENCY_PATH = os.path.join(MATRIX_DATA_FILES_PATH, 'pdb_adjacency_data')
PROTEIN_FEATURE_PATH = os.path.join(MATRIX_DATA_FILES_PATH, 'pdb_features_data')
PROTEIN_FEATURE_PATH_PDB = os.path.join(MATRIX_DATA_FILES_PATH, 'pdb_features_data_local')


#TODO: fix constants and logging setup throughout the entire package
PVALUE_THRESHOLD = 0.05
# TODO: explain data file naming conventions somewhere


