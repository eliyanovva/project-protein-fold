#This script tests the algorithm accuracy without any protein structural information
#This script must be run from outside the sample_runs folder. ie: While pwd is RF, call python3 sample_runs/No3Di.py

import sys
sys.path.append('../../project-protein-fold/RF/')  
import CombineLigandsProteins
import FixedClassificationModel
import numpy as np
import SmileKmer
import Globals

#Create protein matrix
CombineLigandsProteins.import_final()
protmat = CombineLigandsProteins.AA
ligand_count = 38
proteins_matrix = np.repeat(protmat, repeats = ligand_count, axis = 0)

#Import ligands matrix
ligand_dict = Globals.initialize_ligand_dict()
SmileKmer.importmatrix(ligand_dict, 5, 1084)
ligand_matrix = SmileKmer.ligmat

#Combine ligands and proteins 
final_matrix = np.concatenate((np.array(proteins_matrix, dtype = np.uint8), np.array(ligand_matrix, dtype = np.uint8)), axis = 1)

classification = CombineLigandsProteins.Y

FixedClassificationModel.train(final_matrix, classification)