from moleculekit.molecule import Molecule
from moleculekit.tools.voxeldescriptors import getVoxelDescriptors, viewVoxelFeatures
from moleculekit.tools.atomtyper import prepareProteinForAtomtyping
from moleculekit.smallmol.smallmol import SmallMol
from moleculekit.home import home
from moleculekit.tools.atomtyper import metal_atypes
import os
import torch
import torch.nn as nn
import torch.nn.functional as F
import numpy as np
from torch.utils.data import Dataset, DataLoader, random_split
import re
import pandas as pd
import sys

# Does all preprocessing for inputs--just need to be converted 
# into tensors with desired settings when imported


def voxelize_ligands():
    """
    Makes a dictionary of ligand names and voxelized ligand objects
    """
    x1 = {}
    os.chdir('/home/users/jvs15/project-protein-fold')

    # Specifies box size and center for the ligand
    box = [20, 20, 20]
    cent = [10, 10, 10]

    for filename in os.listdir('mol_data_files'):
        # Following preprocessing step on Acellera page: 
        # https://software.acellera.com/docs/latest/moleculekit/tutorials/voxelization_tutorial.html?highlight=voxelization

        # Voxelize ligands and make tensors
        chem_name = filename[0:len(filename)-4]
        ligand = SmallMol(f'mol_data_files/{filename}', force_reading=True)
        lig_vox, lig_centers, lig_N = getVoxelDescriptors(ligand, voxelsize = 0.5, buffer=1, boxsize=box, center=cent)
        lig_vox_t = lig_vox.transpose().reshape([1, 8, lig_N[0], lig_N[1], lig_N[2]])
        lig_vox_t = torch.tensor(lig_vox_t.astype(np.float32))
        x1[chem_name] = lig_vox_t

    os.chdir('/home/users/jvs15/project-protein-fold/moleculekit')
    x1np = np.asarray(x1)
    np.save('ligand_data_tensors.npy', x1np)

    print('Wrote ligand data')
    return

def voxelize_proteins():
    """
    Makes a dictionary of protein accessions and voxelized protein objects
    """
    x2 = {}
    num = 1
    os.chdir('/home/users/jvs15/project-protein-fold')

    box = [40, 40, 40]
    cent = [20, 20, 20]

    # Chunk .pdb files for processing
    filenames = os.listdir('pdb_data_files')
    groups = list(chunks(filenames, 100))

    for chunk in groups:
        for filename in chunk:
            # Check that file is in right format
            if (filename[0:3] == 'AF-' and filename[len(filename)-4:len(filename)] == '.pdb'):
                a, b = filename.find('AF-')+3, filename.find('-F1')
                uniprot = filename[a:b]

                # Voxelize protein
                prot = Molecule(f'pdb_data_files/{filename}')
                prot.filter("protein or water or element {}".format(" ".join(metal_atypes)))
                prot = prepareProteinForAtomtyping(prot)
                prot_vox, prot_centers, prot_N = getVoxelDescriptors(prot, buffer = 1, boxsize=box, center=cent)

                # Make tensor from voxelized format
                prot_vox_t = prot_vox.transpose().reshape([1, 8, prot_N[0], prot_N[1], prot_N[2]])
                prot_vox_t = torch.tensor(prot_vox_t.astype(np.float32))

                # Add tensor to dict
                x2[uniprot] = prot_vox_t

        os.chdir('/home/users/jvs15/project-protein-fold/moleculekit/prot_files')
        x2np = np.asarray(x2)
        np.save(f'protein_data_tensors{num}.npy', x2np)
        x2 = {}
        os.chdir('/home/users/jvs15/project-protein-fold')
        num += 1
    return

def gen_outputs(fdr, logfc):
    # There are about 49k total interactions from listed ligands and proteins so far;
    # 2% of those have FDR < 0.05 and 4% have FDR < 0.25
    """
    @fdr: fdr cutoff threshold
    @logfc: logfc cutoff threshold
    Makes a dataframe containing protein and ligand names, their tensors, and output
    """

    os.chdir('/home/users/jvs15/project-protein-fold/moleculekit')
    lig = np.load('ligand_data_tensors.npy', allow_pickle=True)
    lig = lig.item()
    print('Ligand data loaded')

    prot = {}
    for filename in os.listdir('prot_files'):
        d = np.load(f'prot_files/{filename}', allow_pickle=True)
        d = d.item()
        prot = {**d, **prot}
    print('Protein data loaded')

    de = np.load('de_data.npy', allow_pickle=True)
    de = de.item()
    print('Expression data loaded')

    ret = []
    num = 1
    os.chdir('/home/users/bmp40/project-protein-fold/cnn_scripts/model_files')
    for i in lig.keys():
        for j in prot.keys():
            a = de.get((i, j))
            if a.iloc[0]['FDR'] < fdr:
                if abs(a.iloc[0]['logFC']) > logfc:
                    y = 1
                else:
                    y = 0
                
                entry = {'Ligand Name' : i,
                'Ligand Data' : lig.get(i),
                'Protein Name' : j,
                'Protein Data' : prot.get(j),
                'Binding' : y}

                ret.append(entry)
                if len(ret) == 100:
                    ret_np = np.asarray(ret)
                    np.save(f'all_data{num}.npy', ret_np)
                    ret = []
                    num+=1

    ret_np = np.asarray(ret)
    np.save(f'all_data{num}.npy', ret_np)

    return        

def chunks(lst, n):
    for i in range(0, len(lst), n):
        yield lst[i:i+n]


if __name__ == "__main__":
    gen_outputs(0.05, 1)
    #voxelize_proteins()
    #print(list(chunks([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10], 4)))