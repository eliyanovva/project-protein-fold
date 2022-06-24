#This script tests the model on only one protein with all ligands.
#Protein is A0A1L1SQF6
#This script must be run from outside the sample_runs folder. ie: While pwd is RF, call python3 sample_runs/OneProteinRF.py

import sys
sys.path.append('../../project-protein-fold/RF/')
import CombineLigandsProteins
import FixedClassificationModel
import SmileKmer
import Globals

#Create ligand matrix
ligand_dict = Globals.initialize_ligand_dict()
SmileKmer.importmatrix(ligand_dict, 5, 1)
ligand_matrix = SmileKmer.ligmat

#Create classification vector
CombineLigandsProteins.import_final()
classified = CombineLigandsProteins.dictionary
logFCmat = []
for ligand in list(ligand_dict.keys()):
    logFCmat.append(float(classified[str('A0A1L1SQF6')][ligand]))

FixedClassificationModel.train(ligand_matrix, logFCmat)
