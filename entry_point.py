import os
import sys

# Add repo to path!
MAIN_PACKAGE_DIR = os.path.abspath(os.curdir)    
sys.path.append(MAIN_PACKAGE_DIR)

from cli_arguments import ModelingParser
from graph_cnn.run_model import runModel
from graph_cnn.hp_model import optimizeHyperparameters


#def ppp():
    
parser = ModelingParser()
parser.setup_arguments()
args = parser.parse_args()

batch_size = args.batch_size

if args.model == 'gnn':
    if args.hptuning:
        print('HERE')
        optimizeHyperparameters()
    else:
        runModel()

#ppp()
    