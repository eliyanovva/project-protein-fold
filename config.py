import logging
import sys
import os
from tensorboard.plugins.hparams import api as hp

# Create logging directory
if not os.path.exists('./logs'):
    os.makedirs('./logs')

    
logging.basicConfig(
    encoding='utf-8', 
    level=logging.INFO, 
    format='%(asctime)s %(message)s',
    handlers=[
        logging.FileHandler("logs/script.log"),
        logging.StreamHandler(sys.stdout)
    ]
)


PROTEIN_ADJACENCY_MAT_SIZE = 3000
LIGAND_ADJACENCY_MAT_SIZE = 70
PROTEIN_FEATURES_COUNT = 10
LIGAND_FEATURES_COUNT = 9
ATOM_DICT = {'C':0, 'O':1, 'N':2, 'S':3, 'H':4}

#absolute path for a file
MAIN_PACKAGE_DIR = os.path.abspath(os.curdir)
#final directory in filepath
last_dir = os.path.split(os.path.abspath(os.curdir))[-1]
#ensures the final directory is project-protein-fold
while last_dir != "project-protein-fold":
    os.chdir("..")
    MAIN_PACKAGE_DIR = os.path.abspath(os.curdir)
    last_dir = os.path.split(MAIN_PACKAGE_DIR)[-1]

#script_dir allows DATA_FILES_PATH to be called from any computer
DATA_FILES_PATH = os.path.join(MAIN_PACKAGE_DIR, 'data_files')
PDB_FILES_PATH = os.path.join(DATA_FILES_PATH, 'pdb_data_files')
BGF_FILES_PATH = os.path.join(DATA_FILES_PATH, 'bgf_files')
MOL_FILES_PATH = os.path.join(DATA_FILES_PATH, 'mol_files')
SMILES_FILES_PATH = os.path.join(DATA_FILES_PATH, 'smiles_files')

MATRIX_DATA_FILES_PATH = os.path.join(MAIN_PACKAGE_DIR, 'graph_cnn', 'data_prep')
# MOL is for LIGAND
MOL_ADJACENCY_PATH = os.path.join(MATRIX_DATA_FILES_PATH, 'mol_adjacency_data')
LIGAND_FEATURE_PATH = os.path.join(MATRIX_DATA_FILES_PATH, 'mol_feature_data')
PROTEIN_ADJACENCY_PATH = os.path.join(MATRIX_DATA_FILES_PATH, 'pdb_adjacency_data')
PROTEIN_FEATURE_PATH = os.path.join(MATRIX_DATA_FILES_PATH, 'pdb_features_data')
MATRIX_FOLDERS = [
    PROTEIN_ADJACENCY_PATH, PROTEIN_FEATURE_PATH, MOL_ADJACENCY_PATH, LIGAND_FEATURE_PATH
]

#TODO: fix constants and logging setup throughout the entire package
PVALUE_THRESHOLD = 0.05
# TODO: explain data file naming conventions somewhere

HP_BATCH_SIZE = hp.HParam('batch_size', hp.Discrete([64, 128, 256, 512]))
HP_DROPOUT = hp.HParam('dropout', hp.Discrete([0.15, 0.2, 0.25]))
HP_OPTIMIZER = hp.HParam('optimizer', hp.Discrete(['adam', 'adamax']))
HP_LOSS = hp.HParam('loss', hp.Discrete(['MeanSquaredLogarithmicError', 'MeanAbsoluteError', 'MeanSquaredError']))
HP_LEARNINGRATE = hp.HParam('learning_rate', hp.Discrete([0.001, 0.0001, 0.00001]))
HP_VALIDATION_SPLIT = hp.HParam('validation_split', hp.Discrete([0.1, 0.15, 0.2]))
HP_TEST_TRAIN_SPLIT = hp.HParam('test_train_split', hp.Discrete([0.1, 0.15, 0.2]))

METRIC_ACCURACY = 'accuracy'
