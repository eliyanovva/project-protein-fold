import os
import shutil
import sys
import subprocess
import string

# Add repo to path!
MAIN_PACKAGE_DIR = os.path.abspath(os.curdir)    
sys.path.append(MAIN_PACKAGE_DIR)

from graph_cnn.data_prep import data_generator
from cli_arguments import ModelingParser
from graph_cnn.run_model import runModel
try:
    from graph_cnn.hp_model import optimizeHyperparameters
except:
    pass
import config

def createTemporaryDirectories():
    os.mkdir('temp_ligand_mol')
    os.mkdir('temp_protein_bgf')
    os.mkdir('temp_ligand_adj_npy')
    os.mkdir('temp_ligand_feat_npy')
    os.mkdir('temp_protein_adj_npy')
    os.mkdir('temp_protein_feat_npy')


def removeTemporaryDirectories():
    # FIXME make it catch exception when the directory doesn't exist
    shutil.rmtree('temp_ligand_mol')
    shutil.rmtree('temp_protein_bgf')
    shutil.rmtree('temp_ligand_adj_npy')
    shutil.rmtree('temp_ligand_feat_npy')
    shutil.rmtree('temp_protein_adj_npy')
    shutil.rmtree('temp_protein_feat_npy')


def createTempBGFfile(pdb_protein_path):
    subprocess.run(
        ["obabel " + pdb_protein_path + " -O " + os.path.join("temp_protein_bgf", "input_protein.bgf")],
        shell=True
    )


def createTempMOLfile(mol_ligand_path):
    subprocess.run(
        ["cp " + mol_ligand_path + " " + os.path.join("temp_ligand_mol")],
        shell=True
    )


def ppp():    
    parser = ModelingParser()
    parser.setup_arguments()
    args = parser.parse_args()

    if args.batch_size:
        batch_size = args.batch_size
    else:
        batch_size = -1
    
    if args.model == 'gnn':
        if args.gnn_mode == 'hptuning':
            optimizeHyperparameters()
        elif args.gnn_mode == 'eval_tuple':
            protein_file = args.pdb_file
            ligand_file = args.mol_file
            createTemporaryDirectories()
            createTempBGFfile(protein_file)
            createTempMOLfile(ligand_file)
            data_generator.generateProteinAdjacencyMatrices(
                folder_path='temp_protein_bgf',
                folder_path='temp_protein_adj_npy'
            )
            data_generator.generateProteinFeatureMatrices(
                folder_path='temp_protein_bgf',
                folder_path='temp_protein_feat_npy'
            )

            data_generator.generateLigandAdjacencyMatrices(
                folder_path='temp_ligand_mol',
                folder_path='temp_ligand_adj_npy'
            )
            data_generator.generateLigandFeatureMatrices(
                folder_path='temp_ligand_mol',
                folder_path='temp_ligand_feat_npy'
            )
            
        elif args.gnn_mode == 'eval_protein':
            pass
        elif args.gnn_mode == 'eval_ligand':
            pass
        else:
            runModel(batch_size)
    
    elif args.model == 'cnn':
        print('CNN CLI is not implemented yet!')
    
    elif args.model == 'rf':
        print('RF CLI is not implemented yet!')

#createTemporaryDirectories()
#createTempMOLfile("/home/users/tep18/new_ppp/project-protein-fold/data_files/mol_files/pyridine_pyridine.mol")

removeTemporaryDirectories()
#ppp()
    