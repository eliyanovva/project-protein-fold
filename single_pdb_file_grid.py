""" This module extracts the coordinates from a single PDB file and populates a 4D
grid containing 3D grids for the locations of all Carbon, Oxygen, Nitrogen, and Sulphur atoms.
"""

from pkg_resources import require
import numpy as np
import matplotlib.pyplot as plt
import logging as log
import argparse

import log_config


"""
    TODO: Remove this part later.
    NOTES:
    -> Decided to use only numpy, without pandas because iteration over pandas
    is not a good practice.
    -> numpy can work only with integer indices
    -> pass data folder as a parameter to this script
    -> move argparse setting to a separate file, just like logging
    -> current issue - the PDB grid is too big and too slowly processed
"""
ATOM_DICT = {'C':0, 'O':1, 'N':2, 'S':3}
# MAKE SURE TO USE THE COMPLETE PATH TO THE PDB FILE
PDB_FILE_NAME = 'pdb_data_files/AF-Q6IEV9-F1-model_v2.pdb'


### Arguments parsing part
parser = argparse.ArgumentParser()
parser.add_argument(
    '-file_name',
     help='Complete path of the PDB file to be converted to a grid',
     type=str
    )
args = parser.parse_args()

def extract_PDB_coordinates_atoms(pdb_file_name=PDB_FILE_NAME, grid_size=500, shift=30):
    """This function extracts 3D coordinates of each atom from the PDB file and
    populates the 4D grid, where grid[0] is the carbon layer, 

    Args:
        pdb_file_name (str): The name of the PDB file
        grid_size (int, optional): The default size of the grid. Defaults to 700.

    Returns:
        protein_grid: a 4D grid shaped like a cube, where protein_grid[0] is the carbon layer,
        protein_grid[1] is the oxygen layer, protein_grid[2] is the nitrogen layer,
        and protein_grid[3] is the sulphur layer. 
    """

    log.info('extract_PDB_coordinate_atoms was called with ' + pdb_file_name)
    with open(pdb_file_name, 'r') as pdb_file:
        
        pdb_lines = pdb_file.readlines()
        protein_grid = np.zeros((4, grid_size, grid_size, grid_size))
        shift_y, shift_x, shift_z = 0,0,0
        for line in pdb_lines:
            line = line.strip()
            if line[0:4] == 'ATOM' and line[-1] != 'H':
                numerical_data = line[27:67].split(' ')
                numerical_data = [float(x) for x in numerical_data if x]
                atom_type = line[-1]
                entry = np.array(
                    [ATOM_DICT[atom_type], # atom type (C, O, N, S)
                    numerical_data[0], # x coordinate
                    numerical_data[1], # y_coordinate
                    numerical_data[2], # z_coordinate
                    numerical_data[4]] # confidence score
                )
                shift_x = min(shift_x, entry[1])
                shift_y = min(shift_y, entry[2])
                shift_z = min(shift_z, entry[3])
                #log.debug('The following point entry has been added to the grid: ' + str(entry))
                protein_grid[int(entry[0])] \
                [int((entry[1] - shift_x)*1000//200)] \
                [int((entry[2] - shift_y)*1000//200)] \
                [int((entry[3] - shift_z)*1000//200)] = entry[4]
    log.info('Extraction completed! ')
    
    return protein_grid

def get_protein_name_from_path(file_path):
    left_index = file_path.rfind('/') + 1
    right_index = file_path.find('-model')
    return file_path[left_index : right_index]

def visualize_grid(protein_grid, atom_layer):
    """This function visualizes a single atom layer of the 4D grid and outputs it to a PNG file.

    Args:
        protein_grid (np.array): 4D grid of atom locations
        atom_layer (char): specifies whether the carbon, oxygen, nitrogen, or sulphur layer
            should be graphed.
        TODO: ensure that atom_layer is a correct value
    """
    fig = plt.figure(figsize=(8, 6), dpi=80)
    ax = fig.add_subplot(211, projection='3d')
    z, x, y = protein_grid[ATOM_DICT[atom_layer]].nonzero()
    ax.scatter(x, y, z, alpha=1, marker='.', s=2)
    ax.grid(True)
    ax.set_title('All ' + atom_layer + ' atoms')
    plt.savefig('visuals/aa_grid_' + atom_layer + '.png')


with open('pdb_data_files/file_names.txt') as file_names:
    names = file_names.readlines()
    for name in names:
        complete_file_name = 'pdb_data_files/' + name[0:-1]
        coordinates_data = extract_PDB_coordinates_atoms(complete_file_name)
        protein_name = get_protein_name_from_path(complete_file_name)
        log.info('got here!')
        flat_grid = np.ndarray.ravel(coordinates_data)
        log.info('Your grid was flattenned to a 1D vector!')
        np.save('pdb_data_files/one_dim_grids/' + protein_name + '_grid', flat_grid)
        log.info('A flattened 1D version of your grid has been saved!')

# How to unload np array from the files:
#arr = np.load('pdb_data_files/one_dim_grids/AF-Q6IEV9-F1_grid.npy')
#log.info('1d array has been loaded.')

#        visualize_grid(coordinates_data, 'C')
#        visualize_grid(coordinates_data, 'O')
#        visualize_grid(coordinates_data, 'S')
#        visualize_grid(coordinates_data, 'N')
#        break