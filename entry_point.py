import os
import sys

# Add repo to path!
MAIN_PACKAGE_DIR = os.path.abspath(os.curdir)    
sys.path.append(MAIN_PACKAGE_DIR)

from cli_arguments import ModelingParser
from graph_cnn.run_model import runModel
try:
    from graph_cnn.hp_model import optimizeHyperparameters
except:
    pass
import config


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
            pass
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

ppp()
    