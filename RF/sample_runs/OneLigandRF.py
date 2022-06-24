#This script tests the model on only one ligand with all proteins.
#Ligand is citronellol
#This script must be run from outside the sample_runs folder. ie: While pwd is RF, call python3 sample_runs/OneLigandRF.py

import sys
sys.path.append('../../project-protein-fold/RF/')
import CombineLigandsProteins
import FixedClassificationModel

#Import protein matrix
CombineLigandsProteins.import_final()
protmat = CombineLigandsProteins.protmat

#Create logFC vector
proteins = CombineLigandsProteins.proteins
classified = CombineLigandsProteins.dictionary
logFCmat = []
for protein in proteins:
    logFCmat.append(float(classified[str(protein.name)]['citronellol']))

FixedClassificationModel.train(protmat, logFCmat)