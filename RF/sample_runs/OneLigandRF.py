#This script tests the model on only one ligand with all proteins.

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