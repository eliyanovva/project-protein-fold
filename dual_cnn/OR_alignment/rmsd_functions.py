import sys
import numpy as np
import plotly.graph_objects as go
import copy 
import pandas as pd 
import os


def adj_coordinates(uncentered_dir, df, out_dir):
    """
    Generate new pdb files with corrected coordinates stored in coord_file
    @coord_file: dataframe output from rmsd.dict_to_pdDataFrame
    @uncentered_dir: folder containing uncentered pdb files
    @out_dir: directory to store new .pdb files
    """

    os.chdir(uncentered_dir)
    for filename in os.listdir('.'):
        print(f'Working on {filename}')
        with open(filename, 'r') as f:
            x_column = None
            lines = f.readlines()
            
            if 'olfr' in filename.lower():
                pid = filename.split('_')[0]
                if '.pdb' in pid:
                    pid = pid.split('.pdb')[0]
            elif 'AF' in filename:
                pid = filename.split('-')[1]

            startind = (df.models.values == pid).argmax()
            idx = startind
            
            with open(f'{out_dir}/{filename}', 'w') as fout:
                for line in lines:
                    if line.startswith("TER") or line.startswith("END"):
                        break
                    if line.startswith("ATOM"):
                        # Might need to be fixed for lowercase entries
                        tokens = line.split()
                        if x_column is None:
                            try:
                                # look for x column
                                for i, x in enumerate(tokens):
                                    if "." in x and "." in tokens[i + 1] and "." in tokens[i + 2]:
                                        x_column = i
                                        break

                            except IndexError:
                                msg = f"error: Parsing coordinates for {line} in {filename}"
                                exit(msg)

                        # If statement is a safeguard
                        if df.iloc[[idx]].models.iloc[0] != pid:
                            print(f'Lines do not match length of pdb')
                            break
                        else:
                            tokens[x_column] = df.iloc[[idx]].x.iloc[0]
                            tokens[x_column + 1] = df.iloc[[idx]].y.iloc[0]
                            tokens[x_column + 2] = df.iloc[[idx]].z.iloc[0]
                        for i in range(len(tokens)):
                            if i == len(tokens) - 1:
                                fout.write(str(tokens[i]) + '\n')
                            else:
                                fout.write(str(tokens[i]) + '  ')
                        idx += 1


def get_coordinates_pdb(filename, is_gzip=False, return_atoms_as_int=False):
    """
    Get coordinates from the first chain in a pdb file
    and return a vectorset with all the coordinates.

    Parameters
    ----------
    filename : string
        Filename to read

    Returns
    -------
    atoms : list
        List of atomic types
    V : array
        (N,3) where N is number of atoms
    """

    # PDB files tend to be a bit of a mess. The x, y and z coordinates
    # are supposed to be in column 31-38, 39-46 and 47-54, but this is
    # not always the case.
    # Because of this the three first columns containing a decimal is used.
    # Since the format doesn't require a space between columns, we use the
    # above column indices as a fallback.
    x_column = None
    V = list()

    # Same with atoms and atom naming.
    # The most robust way to do this is probably
    # to assume that the atomtype is given in column 3.

    atoms = list()
    resid = list()

    if is_gzip:
        openfunc = gzip.open
        openarg = "rt"
    else:
        openfunc = open
        openarg = "r"

    with openfunc(filename, openarg) as f:
        lines = f.readlines()
        for line in lines:
            if line.startswith("TER") or line.startswith("END"):
                break
            if line.startswith("ATOM"):
                tokens = line.split()
                # Try to get the atomtype
                try:
                    atom = tokens[2]
                    atoms.append(atom)
                except ValueError:
                    msg = f"error: Parsing atomtype for the following line:" f" \n{line}"
                    exit(msg)

                if x_column is None:
                    try:
                        # look for x column
                        for i, x in enumerate(tokens):
                            if "." in x and "." in tokens[i + 1] and "." in tokens[i + 2]:
                                x_column = i
                                break

                    except IndexError:
                        msg = "error: Parsing coordinates " "for the following line:" f"\n{line}"
                        exit(msg)

                # Try to read the coordinates
                try:
                    V.append(np.asarray(tokens[x_column : x_column + 3], dtype=float))
                
                except ValueError:
                    # If that doesn't work, use hardcoded indices
                    try:
                        x = line[30:38]
                        y = line[38:46]
                        z = line[46:54]
                        V.append(np.asarray([x, y, z], dtype=float))
                    except ValueError:
                        msg = "error: Parsing input " "for the following line:" f"\n{line}"
                        exit(msg)
                # Try to read the resid number
                try:
                    resid.append(np.asarray(tokens[5], dtype=int))
                except ValueError:
                    msg = "error while reading resid - HL"
                    exit(msg)


#    if return_atoms_as_int:
#        atoms = [int_atom(atom) for atom in atoms]

    V = np.asarray(V)
    atoms = np.asarray(atoms)
    resid = np.asarray(resid)


    
    assert V.shape[0] == atoms.size 

    return atoms, V, resid


def rmsd(V, W):
    """
    Calculate Root-mean-square deviation from two sets of vectors V and W.

    Parameters
    ----------
    V : array
        (N,D) matrix, where N is points and D is dimension.
    W : array
        (N,D) matrix, where N is points and D is dimension.

    Returns
    -------
    rmsd : float
        Root-mean-square deviation between the two vectors
    """
    diff = np.array(V) - np.array(W)
    N = len(V)
    return np.sqrt((diff * diff).sum() / N)

def distance(V, W):
    """
    Calculate Root-mean-square deviation from two points V and W.

    rmsd : float
        Root-mean-square deviation between the two vectors
    """
    diff = np.array(V) - np.array(W)
    return np.sqrt((diff * diff).sum())


def centroid(X):
    """
    Centroid is the mean position of all the points in all of the coordinate
    directions, from a vectorset X.

    https://en.wikipedia.org/wiki/Centroid

    C = sum(X)/len(X)

    Parameters
    ----------
    X : array
        (N,D) matrix, where N is points and D is dimension.

    Returns
    -------
    C : float
        centroid
    """
    C = X.mean(axis=0)
    return C


def kabsch(P, Q):
    """
    Using the Kabsch algorithm with two sets of paired point P and Q, centered
    around the centroid. Each vector set is represented as an NxD
    matrix, where D is the the dimension of the space.
    The algorithm works in three steps:
    - a centroid translation of P and Q (assumed done before this function
      call)
    - the computation of a covariance matrix C
    - computation of the optimal rotation matrix U
    For more info see http://en.wikipedia.org/wiki/Kabsch_algorithm
    Parameters
    ----------
    P : array
        (N,D) matrix, where N is points and D is dimension.
    Q : array
        (N,D) matrix, where N is points and D is dimension.
    Returns
    -------
    U : matrix
        Rotation matrix (D,D)
    """

    # Computation of the covariance matrix
    C = np.dot(np.transpose(P), Q)

    # Computation of the optimal rotation matrix
    # This can be done using singular value decomposition (SVD)
    # Getting the sign of the det(V)*(W) to decide
    # whether we need to correct our rotation matrix to ensure a
    # right-handed coordinate system.
    # And finally calculating the optimal rotation matrix U
    # see http://en.wikipedia.org/wiki/Kabsch_algorithm
    V, S, W = np.linalg.svd(C)
    d = (np.linalg.det(V) * np.linalg.det(W)) < 0.0

    if d:
        S[-1] = -S[-1]
        V[:, -1] = -V[:, -1]

    # Create Rotation matrix U
    U = np.dot(V, W)

    return U

def angle(CA, CX, Center):
    """
    Calculate angle at CA pointing towards CA-center and CA-CX. 

    Parameters
    ----------
    A : array
        [x, y, z] coordinates.
    B : array
        [x, y, z] coordinates.
    Center : array
        [x, y, z] coordinates.


    Returns
    -------
    angle : float
        Degree angle at center between A and B 
    """
    
    v1 = [CX[0] - CA[0], CX[1] - CA[1], CX[2] - CA[2]]
    v2 = [Center[0] - CA[0], Center[1] - CA[1], Center[2] - CA[2]]
    v1mag = np.sqrt(v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2])
    v1norm = [v1[0] / v1mag, v1[1] / v1mag, v1[2] / v1mag]
    v2mag = np.sqrt(v2[0] * v2[0] + v2[1] * v2[1] + v2[2] * v2[2])
    v2norm = [v2[0] / v2mag, v2[1] / v2mag, v2[2] / v2mag]
    res = v1norm[0] * v2norm[0] + v1norm[1] * v2norm[1] + v1norm[2] * v2norm[2]
    
    angle = np.arccos(res)*(180/np.pi)
    
    return angle

def dict_plot( coordinate , mode='markers+lines', size = 3):
    x = []
    y = []
    z = []
    for i in coordinate: 
        x.append(i[0])
        y.append(i[1])
        z.append(i[2])

    fig = go.Figure()
    fig.add_trace(go.Scatter3d(
        x = x, 
        y = y, 
        z = z,
        mode = mode
    ))
    fig.update_traces( marker=dict(size=size, opacity = 0.6))
    
    return fig

def dict_plot_all( data , mode='markers+lines', size = 3):
    x = {}
    y = {}
    z = {}
    for i in data: 
        x[models] = []
        y[models] = []
        z[models] = []
        for i in data[models]['coord']: 
            x[models].append(i[0])
            y[models].append(i[1])
            z[models].append(i[2])
    fig = go.Figure()
    for models in x: 
        fig.add_trace(go.Scatter3d(
            x = x[models], 
            y = y[models], 
            z = z[models],
            mode = mode,
            name = models
    ))
    fig.update_traces( marker=dict(size=size, opacity = 0.4))
    # update_layout setting the axis visibility and background to False 
    fig.update_layout(width=600, height=600,
                  scene = dict(
                    xaxis = dict(
                         visible= False,
                         showbackground=False),
                    yaxis = dict(
                         visible= False,
                         showbackground=False),
                    zaxis = dict(
                         visible= False,
                         showbackground=False)),
                    margin=dict(r=10, l=10, b=10, t=10)
    )
    return fig

def df_plot_all(data_df, size=3, opacity=0.6, no_background=True):
    fig = go.Figure()
    for models in np.unique(data_df.models):
        fig.add_trace(go.Scatter3d(
            x = data_df.iloc[np.where(data_df.models == models)].x, 
            y = data_df.iloc[np.where(data_df.models == models)].y, 
            z = data_df.iloc[np.where(data_df.models == models)].z, 
            name = models, 
            mode = 'markers'
        ))
    fig.update_traces( marker=dict(size=size, opacity=opacity))

    # update_layout setting the axis visibility and background to False 
    if no_background:
        fig.update_layout(width=600, height=600,
                          scene = dict(
                            xaxis = dict(
                                 visible= False,
                                 showbackground=False),
                            yaxis = dict(
                                 visible= False,
                                 showbackground=False),
                            zaxis = dict(
                                 visible= False,
                                 showbackground=False)),
                            margin=dict(r=10, l=10, b=10, t=10)
    )
    fig.show()
    return


# remove "H" from models, since some models doesn't have it and it may throw off alignment 
def remove_H(data):
    for models in data:     
        H_atom = ["H" in i for i in data[models]['atom']]
        non_H_atom = np.invert(H_atom)
        data[models]['atom'] = data[models]['atom'][non_H_atom]
        data[models]['coord'] = data[models]['coord'][non_H_atom]
        data[models]['resid'] = data[models]['resid'][non_H_atom]
        print('Removed H atoms from', models)
    return data
        
# reorder the atoms, coord and resid so that atom orders are ('C' 'CA' 'CB' 'N' 'O' ...)
# sorted by alphabatical order for simplicity 
def reorder_atoms(data):
    for models in data: 
        orig_index = []
        reordered_index = []
        for resid in np.unique(data[models]['resid']):
            (orig_index,) = np.where(data[models]['resid'] == resid)
            reordered_index = np.argsort(data[models]['atom'][orig_index]) + np.min(orig_index)
            data[models]['atom'][orig_index] = copy.deepcopy(data[models]['atom'][reordered_index])
            data[models]['coord'][orig_index] = copy.deepcopy(data[models]['coord'][reordered_index])
        print('Reordered',models,'atoms :', data[models]['atom'][0:10])
    return data

# receives data in dictionary form and returns only N, CA, C and O atoms
def get_backbone(data, reorder=True):
    if reorder:
        data = reorder_atoms(data)
    for models in data: 
        N_atom = ["N" == i for i in data[models]['atom']]
        CA_atom = ["CA" == i for i in data[models]['atom']]
        C_atom = ["C" == i for i in data[models]['atom']]
        O_atom = ["O" == i for i in data[models]['atom']]
        # intersect all the True index 
        # extracts only the indexes of ["N","CA","CB","C","O"] for backbone data
        temp = np.logical_or(N_atom,CA_atom)
        temp = np.logical_or(temp,C_atom)
        temp = np.logical_or(temp,O_atom)
        for i in data[models]:
            data[models][i] = data[models][i][temp]
        print('Backbone atoms extracted for ',models)
    return data

# Manually subset the resid numbers for RMSD calculation, since calculating RMSD requires the atom numbers to be same 
# typically omit the ends (resid 25 - 285)
def trim_resid(data, start = 25, end = 285):
    backbone = copy.deepcopy(data)
    for models in backbone: 
        keep_resid = np.where((backbone[models]['resid'] <= end) & (backbone[models]['resid'] >= start ))
        backbone[models]['atom'] = backbone[models]['atom'][keep_resid]
        backbone[models]['coord'] = backbone[models]['coord'][keep_resid]
        backbone[models]['resid'] = backbone[models]['resid'][keep_resid]
        print('Resid trimed from',start, 'to', end, 'for', models)
    return backbone

# rmsd_data takes data in and returns rmsd_data dictionary for rmsd at every resid for all models
# all models in data needs to have the same length of resid 
def get_rmsd(data):
    align = list(data.keys())[0]
    rmsd_data = {}
    rmsd_resid = np.unique(data[align]['resid'])
    for models in data: 
        if models == align: 
            continue 
        rmsd_data[models] = []
        for resid in np.unique(data[models]['resid']):
            A_resid_index = np.where(data[models]['resid'] == resid)
            B_resid_index = np.where(data[align]['resid'] == resid)
            temp = rmsd(data[models]['coord'][A_resid_index], data[align]['coord'][B_resid_index])
            rmsd_data[models].append(temp)
        print('rmsd calculated by aligning', models, 'to', align)
    return rmsd_data

# rmsd_table takes in sequence, resid array, rmsd_data for the corresponding array and concatenate into a pd data table
def rmsd_table(data, sequence):
    rmsd_list = get_rmsd(data)
    resid = np.unique(data[list(data.keys())[0]]['resid'])
    aa = []
    for i in resid:
        aa.append(list(sequence)[i-1])
    table = pd.DataFrame(rmsd_list)
    table['resid'] = resid
    table['aa'] = aa
    print('Rmsd table contructed')
    return table

# Get rotation matrix, can only compare x to y. The first key in dict{} is used 
def align_data(data, align=None, rotate_data = None):
    remove = []
    if align == None :
        align = list(data.keys())[0]
    for models in data: 
        if models == align: 
            continue 
        if len(data[models]['coord']) == len(data[align]['coord']):
            rotation = kabsch(data[models]['coord'], data[align]['coord'])
        else: 
            print(models + ' NOT calculated, removed from data') 
            remove.append(models)
            continue
        if rotate_data != None: 
            rotate_data[models]['coord'] = np.dot(rotate_data[models]['coord'], rotation)
            print('Rotation calculated for', align,'and applied to ROTATE DATA ',models)
        else: 
            data[models]['coord'] = np.dot(data[models]['coord'], rotation)
            print('Rotation calculated for', align,'and applied onto',models)
                  
    if rotate_data != None: 
        if remove != []: 
            for i in remove: 
                  rotate_data.pop(i)
                  print('REMOVED '+i)
        print('RETURN new rotated data')
        return rotate_data
    print('RETURN updated rotated data')
    return data

# Center models to center by substracting coordinates by centroid 
def center_data(data):
    for models in data: 
        center = centroid(data[models]['coord'])
        data[models]['coord'] -= center
        print('Data centered for',models)
    return data

# Generates the combined sequence. IF consensus and h OR doesn't match on the aa, "X" is used 
def intersect_sequence(sequence_A, sequence_B):
    mutation_count = 0
    intersect_sequence = ""
    for aa in range(len(sequence_A)):
        if sequence_A[aa] == sequence_B[aa]:
            intersect_sequence += sequence_A[aa]
        else: 
            intersect_sequence += "X"
            mutation_count += 1
    print(mutation_count, "amino acid replaced with 'X' ")
    return intersect_sequence

# Recreate backbone_df from data, so that it is not truncated but full length 
# create a backbone_data for measuring center, RMSD and rotation 
def get_full_backbone(data, align_start=25, align_end=285, full_data=False):
    backbone = copy.deepcopy(data)
    backbone = get_backbone(backbone, reorder=True)
    backbone = center_data(backbone)
    # backbone = trim_resid(backbone)
    # backbone = align_data(backbone)

    # # create truncated backbone to get rotation (since you cannot apply kabsch on different length)
    short_backbone = trim_resid(backbone, start=align_start, end=align_end)
    align = list(short_backbone.keys())[0]
    for models in data: 
        if models == align: 
            continue 
        rotation = kabsch(short_backbone[models]['coord'], short_backbone[align]['coord'])
#   if full_data is True return rotation aligned full data instead of backbone
        if full_data:
            data[models]['coord'] = np.dot(data[models]['coord'], rotation)
        else:
            backbone[models]['coord'] = np.dot(backbone[models]['coord'], rotation)
            
    if full_data:
        return data
    return backbone

# for reproducibility and simplicity transform backbone/data dict{} into pd.DataFrame 
def dict_to_pdDataFrame(data):
    data_df = pd.DataFrame()
    x = []
    y = []
    z = []
    models = []
    resids = []
    atoms = []
    for model in data:
        for coord in data[model]['coord']:
            x.append(coord[0])
            y.append(coord[1])
            z.append(coord[2])
        for resid in data[model]['resid']:
            resids.append(resid)
        for atom in data[model]['atom']:
            atoms.append(atom)
            models.append(model)
    data_df['x'] = x
    data_df['y'] = y
    data_df['z'] = z
    data_df['models'] = models 
    data_df['resid'] = resids
    data_df['atom'] = atoms
    
    return data_df

#  find difference between distance of(CA-centroid) - (last C atom - centroid) 
#  if positive pointing inward else outward 
def CA_CX_distance(data_df, center):
    distance_df = pd.DataFrame({
        'resid': list(np.unique(data_df.resid))
    })
    for model in np.unique(data_df.models):
        distance_df[model] = None
        for resid in np.unique(data_df[data_df.models == model].resid): 
            temp = data_df[data_df.resid == resid][['x', 'y', 'z','atom','models']]
            temp = temp.loc[temp.models == model]
            temp = temp[temp.atom.str.contains('C')]
            CA = temp.iloc[1][['x', 'y', 'z']] # extracts x,y,z coordinates of CA 
            CX = temp.iloc[len(temp)-1][['x', 'y', 'z']] # extracts x,y,z coordinates of last C
            # calculates direction by subtracting distance between (CA-centroid) and (CX-centroid)
            direction =  distance(CA, center[model]) - distance(CX, center[model])
            distance_df.loc[distance_df.resid == resid, model] = direction
    # removes rows with null (None) value 
#     distance_df = distance_df.iloc[np.where(distance_df[distance_df.keys()[1]].notnull())]    
    return distance_df


# CA_CX_angle finds the angle between center-CA-CX.  
# The smaller the angle the more direct it's pointing towards the center.
def CA_CX_angle(data_df, center):
    angle_df = pd.DataFrame({
        'resid': list(np.unique(data_df.resid))
    })
    for model in np.unique(data_df.models):
        angle_df[model] = None
        for resid in np.unique(data_df[data_df.models == model].resid): 
            temp = data_df[data_df.resid == resid][['x', 'y', 'z','atom','models']]
            temp = temp.loc[temp.models == model]
            temp = temp.iloc[np.where(temp.atom.str.contains('C'))]
            CA = temp.iloc[1][['x', 'y', 'z']] # extracts x,y,z coordinates of CA 
            CX = temp.iloc[len(temp)-1][['x', 'y', 'z']] # extracts x,y,z coordinates of last C(CX)
            # calculates direction by subtracting distance between (CA-centroid) and (CX-centroid)
            degree = angle(CA, CX, center[model])
            angle_df.loc[angle_df.resid == resid, model] = degree
    # removes rows with null (None) value 
#     angle_df = angle_df.iloc[np.where(angle_df[angle_df.keys()[1]].notnull())]    
    return angle_df


# calculate the angle/distance square dfference of models to cryoEM
def CA_CX_difference(angle_df, distance_df):
    difference_df = pd.DataFrame({
        'resid': list(np.unique(angle_df.resid))
    })
    align = list(angle_df.keys())[1]
    models = list(angle_df.keys())[1:len(angle_df.keys())]
    for model in models:
        if model == align: 
            continue
        difference_df[model+'_angle_to'+align] = None 
        difference_df[model+'_distance_to'+align] = None 
        difference_df[model+'_angle_to'+align] = np.sqrt(list((angle_df[model] - angle_df[align])**2))
        difference_df[model+'_distance_to'+align] = np.sqrt(list((distance_df[model] - distance_df[align])**2))
    # drops Nan rows AND absolute value everything 
#     difference_df = difference_df.dropna().abs()
    return difference_df
# CA_CX_diff_tocryoEM.to_csv('./OR51_CA_CX_diff.csv')

# Returns distance between longest C atom CX to provided center point. 
def CX_center_distance(data_df, center):
    distance_df = pd.DataFrame({
        'resid': list(np.unique(data_df.resid))
    })
    
    for model in np.unique(data_df.models):
        distance_df[model] = None
        for resid in np.unique(data_df[data_df.models == model].resid): 
            temp = data_df[data_df.resid == resid][['x', 'y', 'z','atom','models']]
            temp = temp.loc[temp.models == model]
            temp = temp[temp.atom.str.contains('C')]
            CX = temp.iloc[1][['x', 'y', 'z']] # extracts x,y,z coordinates of last C
            distance_df.loc[distance_df.resid == resid, model] = distance(CX, center[model]) 
    # removes rows with null (None) value 
#     distance_df = distance_df.iloc[np.where(distance_df[distance_df.keys()[1]].notnull())]    
    return distance_df



