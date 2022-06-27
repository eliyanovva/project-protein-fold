import os

ATOMS_COUNT = 3000
MOL_ATOMS_COUNT = 70
PROTEIN_FEATURES_COUNT = 5
ATOM_DICT = {'C':0, 'O':1, 'N':2, 'S':3}

DATA_FILES_PATH = os.path.join('/home', 'users', 'tep18', 'new_ppp', 'project-protein-fold', 'data_files')
PDB_FILES_PATH = os.path.join(DATA_FILES_PATH, 'pdb_data_files')
BGF_FILES_PATH = os.path.join(DATA_FILES_PATH, 'bgf_files')
MOL_FILES_PATH = os.path.join(DATA_FILES_PATH, 'mol_files')
SMILES_FILES_PATH = os.path.join(DATA_FILES_PATH, 'smiles_files')

# TODO: explain data file naming conventions somewhere
