import numpy as np
#import pandas as pd
import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D 

"""
    TODO: Remove this part later.
    NOTES:
    -> Decided to use only numpy, without pandas because iteration over pandas
    is not a good practice.
    -> numpy can work only with integer indices
"""
ATOM_DICT = {'C':0, 'O':1, 'N':2, 'S':3}
    
def extract_PDB_coordinates_atoms(pdb_file_name, grid_size=100, granularity=0.05):
    """column 0 - x
    column 1 - y
    column 2 - z
    column 3 - confidence score
    column 4 - atom type

    The function drops the hydrogen atoms right now.

    PDB files have 

    Args:
        pdb_file_name (string): The complete path of the PDB file to be analyzed.
    """
    
    with open(pdb_file_name, 'r') as pdb_file:
        pdb_lines = pdb_file.readlines()
        atoms_count_line = pdb_lines[-3].split(' ')
        atoms_count_line = [x for x in atoms_count_line if x]
        atoms_count = int(atoms_count_line[1])

        #data = np.zeros((atoms_count - 1, 5))
        protein_grid = np.zeros((4, grid_size, grid_size, grid_size))
    
        line_count = 0
        for line in pdb_lines:
            line = line.strip()
            if line[0:4] == 'ATOM' and line[-1] != 'H':
                numerical_data = line[27:67].split(' ')
                numerical_data = [float(x) for x in numerical_data if x]
                # numerical_data looks like that: ['28.826', '16.271', '2.280', '1.00', '39.86'] 
                atom_type = line[-1]
                entry = np.array(
                    [ATOM_DICT[atom_type], # atom type (C, O, N, S)
                    numerical_data[0], # x coordinate
                    numerical_data[1], # y_coordinate
                    numerical_data[2], # z_coordinate
                    numerical_data[4]] # confidence score
                )
                print(entry)
                #protein_grid[entry[0]*1000//50][entry[1]*1000//50][entry[2]*1000//50][entry[3]*1000//50] = entry[4]
                line_count += 1
    return protein_grid#data

def create_populate_grid(pdb_data, grid_size=1500, granularity=0.05):
    """TODO: Figure out if you want to use pandas at all for this.
    The default size of the grid is 60 x 60 x 60 with granularity of 0.5 A.
    The order of layers is C - O - N - S.
    Args:
        pdb_df (pandas.DataFrame): A dataframe object containing the coordinates,
        atom types, and confidence scores for each atom in a PDB file.
    """
    protein_grid = np.zeros((4, grid_size, grid_size, grid_size))
    for atom in pdb_data:
        protein_grid[pdb_data[0]][pdb_data[1]][pdb_data[2]][pdb_data[3]] = pdb_data[4]
    
    return protein_grid


def visualize_grid(protein_grid):
    fig = plt.figure()
    ax = fig.add_subplot(211, projection='3d')
    ax.plot(protein_grid[0][1], protein_grid[0][2], protein_grid[0][3], '*')
    ax.grid(True)
    ax.set_title('All carbon atoms')
        


pdb_file_name = '/home/users/tep18/AF-Q6IEV9-F1-model_v2.pdb'
coordinates_data = extract_PDB_coordinates_atoms(pdb_file_name)

#protein_grid = create_populate_grid(coordinates_data)

visualize_grid(coordinates_data)
#create_populate_grid(data_df)