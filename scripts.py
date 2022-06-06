""" This module extracts the coordinates from a single PDB file and populates a 4D
grid containing 3D grids for the locations of all Carbon, Oxygen, Nitrogen, and Sulphur atoms.
"""

import numpy as np
import matplotlib.pyplot as plt
import logging as log
import log_config

"""
    TODO: Remove this part later.
    NOTES:
    -> Decided to use only numpy, without pandas because iteration over pandas
    is not a good practice.
    -> numpy can work only with integer indices
"""
ATOM_DICT = {'C':0, 'O':1, 'N':2, 'S':3}
# MAKE SURE TO USE THE COMPLETE PATH TO THE PDB FILE
PDB_FILE_NAME = '/home/users/tep18/AF-Q6IEV9-F1-model_v2.pdb'

def extract_PDB_coordinates_atoms(pdb_file_name, grid_size=700):
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
                #log.debug('The following point entry has been added to the grid: ' + str(entry))
                protein_grid[int(entry[0])] \
                [int(entry[1]*1000//100)] \
                [int(entry[2]*1000//100)] \
                [int(entry[3]*1000//100)] = entry[4]
    log.info('Extraction completed! ')
    
    return protein_grid

def create_populate_grid(pdb_data, grid_size=1500, granularity=0.05):
    """TODO: Figure out if you want to use pandas at all for this.
    The default size of the grid is 60 x 60 x 60 with granularity of 0.5 A.
    The order of layers is C - O - N - S.
    Args:
        pdb_df (pandas.DataFrame): A dataframe object containing the coordinates,
        atom types, and confidence scores for each atom in a PDB file.
    """
    pass
    #protein_grid = np.zeros((4, grid_size, grid_size, grid_size))
    #for atom in pdb_data:
    #    protein_grid[pdb_data[0]][pdb_data[1]][pdb_data[2]][pdb_data[3]] = pdb_data[4]
    # 
    #return protein_grid


def visualize_grid(protein_grid, atom_layer, grid_size=700):
    fig = plt.figure(figsize=(8, 6), dpi=80)
    ax = fig.add_subplot(211, projection='3d')
    # x_vals = np.arange(0, grid_size)
    # y_vals = np.arange(0, grid_size)
    z, x, y = protein_grid[ATOM_DICT[atom_layer]].nonzero()
    ax.scatter(x, y, z, alpha=1, marker='.', s=2)
    ax.grid(True)
    ax.set_title('All ' + atom_layer + ' atoms')
    plt.savefig('aa_grid_' + atom_layer + '.png')
        


pdb_file_name = '/home/users/tep18/AF-Q6IEV9-F1-model_v2.pdb'
coordinates_data = extract_PDB_coordinates_atoms(pdb_file_name)

print(coordinates_data)

#protein_grid = create_populate_grid(coordinates_data)

visualize_grid(coordinates_data, 'C')
visualize_grid(coordinates_data, 'O')
visualize_grid(coordinates_data, 'S')
visualize_grid(coordinates_data, 'N')

#x = np.linspace(0, 20, 100)
#plt.plot(x, np.sin(x))
#plt.savefig('test_fig.png')

#create_populate_grid(data_df)