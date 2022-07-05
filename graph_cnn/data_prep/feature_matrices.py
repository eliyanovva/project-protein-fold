import openbabel
from openbabel import pybel
from constants import PDB_FILES_PATH
import logging as log
import log_config
import os

def getMatrix(filetype, filepath):
    #reading pdb or mol file
    filepath = os.path.join(PDB_FILES_PATH, filepath)
    log.info('filetype in format ' + filetype)
    log.info('filepath is ' + filepath)
    prefile = next(pybel.readfile(filetype, filepath))

    #creating a Molecule object
    molecule = pybel.Molecule(prefile)

    #building the feature matrix
    # rows are features, columns are the atoms in the molecule
    #using all the features - will figure out which are impactful later
    feature_matrix = []
    for x in molecule:
        row = []
        atomic_mass = x.atomicmass
        row.append(atomic_mass)
        #coordinateidx - possible to do coordinate index,
        #molecules are not lined up, so I am not adding it
        #plus location is represented in the adjacency matrix
        #notably other location based features are present,
        #similarly I am not adding them
        exact_mass = x.exactmass
        row.append(exact_mass)
        #unsure how this is different than atomic mass,
        #will look into it
        formal_charge = x.formalcharge
        row.append(formal_charge)
        heavy_degree = x.heavydegree
        row.append(heavy_degree)
        #number of non-hydrogen atoms attached
        hetero_degree = x.heterodegree
        row.append(hetero_degree)
        #number of heteroatoms(non Carbons or Hydrogens) attached
        hybridization = x.hyb
        row.append(hybridization)
        sequence = x.idx
        row.append(sequence)
        #interestingly/horrifyingly the index starts at one
        isotope = x.isotope
        row.append(isotope)
        #isotope for the atom, 0 otherwise
        partial_charge = x.partialcharge
        row.append(partial_charge)
        spin = x.spin
        row.append(spin)
        #spin multiplicity - complex orientation thing
        #link if curious to learn more
        #https://chem.libretexts.org/Bookshelves/Physical_and_Theoretical_Chemistry_Textbook_Maps/Supplemental_Modules_(Physical_and_Theoretical_Chemistry)/Electronic_Structure_of_Atoms_and_Molecules/Evaluating_Spin_Multiplicity
        type = x.type
        row.append(type)
        #atom type
        degree = x.degree
        row.append(degree)
        #number of explicit connections
        feature_matrix.append(row)
    return feature_matrix


#Notes:
# feature ideas - in binding pocket(a boolean)
#use help(protein) to learn about the molecule object
#always creates error "Open Babel Warning  in PerceiveBondOrders
#  Failed to kekulize aromatic bonds in OBMol::PerceiveBondOrders"
#error already reported on github https://github.com/volkamerlab/teachopencadd/issues/180
#does not appear to cause problems
