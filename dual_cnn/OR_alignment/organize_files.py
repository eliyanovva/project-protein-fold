import shutil
import os
import re
import pickle
import pandas as pd

os.chdir('/home/users/jvs15/project-protein-fold/dual_cnn/OR_alignment')
pickle_file = open('existing_symb_to_unip.p', 'rb')

# Keys = gene name | Values = UNIPROT
d = pickle.load(pickle_file)

# Find unique names in gene list and prepare to populate new matching file
os.chdir('/home/users/jvs15/project-protein-fold/dual_cnn')
df = pd.read_csv('odorants_ORs_paper.csv')
names = pd.DataFrame(df.Gene.unique(), columns=['Gene'])
names['File Names'] = 'NaN'

deorphan_dir = '/usr/project/csplus2/dietrich/datafiles/share/deorphanOR_human'
orphan_dir = '/usr/project/csplus2/dietrich/datafiles/share/orphanOR_human'


# Loop through each filename in orphan and deorphan directories and 
# if a uniprot name is a substring in any of those file names, add that
# filename to the resulting dictionary
for k in d.keys():
    if any([d[k] in x for x in os.listdir(deorphan_dir)]):
        d[k] = [x for x in os.listdir(deorphan_dir) if any([d[k] in x])][0]
        row = names.loc[names['Gene'] == f'h{k}'].index[0]
        names.at[row, 'File Names'] = d[k]
        shutil.copy(f'{deorphan_dir}/{d[k]}', f'/home/users/jvs15/project-protein-fold/dual_cnn/OR_alignment/ORs_uncentered/{d[k]}')
    elif any([d[k] in x for x in os.listdir(orphan_dir)]):
        d[k] = [x for x in os.listdir(orphan_dir) if any([d[k] in x])][0]
        row = names.loc[names['Gene'] == f'h{k}'].index[0]
        names.at[row, 'File Names'] = d[k]
        shutil.copy(f'{orphan_dir}/{d[k]}', f'/home/users/jvs15/project-protein-fold/dual_cnn/OR_alignment/ORs_uncentered/{d[k]}')
    else:
        print(f'Could not find pdb file {d[k]}')
        continue

names.to_csv('Gene_FileName_map.csv')
