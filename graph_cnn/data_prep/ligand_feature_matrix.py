import openbabel
from openbabel import pybel
import constants

# importing ligand
filepath = "/Users/austin/Documents/ProteinFoldWork/project-protein-fold/data_files/mol_files/(E)-2-butene-1-thiol_e2butene1thiol.mol"
prefile = next(pybel.readfile("mol", filepath))

#ligand = openbabel.OBMol(prefile)
ligand = pybel.Molecule(prefile)

#will eventually follow same/similar framework as protein_feature_matrix