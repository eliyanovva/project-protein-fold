import pandas as pd
import numpy as np
from moleculekit.molecule import Molecule
from moleculekit.tools.voxeldescriptors import getVoxelDescriptors, viewVoxelFeatures
from moleculekit.tools.atomtyper import prepareProteinForAtomtyping
from moleculekit.smallmol.smallmol import SmallMol
from moleculekit.home import home
import os
import re
from openbabel import pybel

def run():
    os.chdir('/home/users/jvs15/project-protein-fold/moleculekit')
    or_dat = pd.read_csv('odorants_ORs_paper.csv')
    ligand_dat = pd.read_csv('odorants_paper.csv', encoding='latin-1')
    dir = 'mol_cid_files'

    for index, row in ligand_dat.iterrows():
        name = row['PubChem CID']
        new_mol(f'lig_{name}', row['SMILES'], 'mol_cid_files')


def new_mol(name, smile, dir):
    if isinstance(smile, str):
        mol = pybel.readstring("smi", smile)
        mol.make3D()
        with open(f'{dir}/{name}.mol', 'w') as f:
            f.write(mol.write('sdf'))

if __name__ == "__main__":
    run()