import os
from os.path import exists
import shutil

alphafold_path = '/usr/xtmp/derf/alphafold/output/'

print(exists('/usr/xtmp/derf/alphafold/output/q9h208_h43r_h207r_k258t'))

for filename in os.listdir('./new_fasta'):
    src = alphafold_path + filename[:-6] + '/ranked_0.pdb'
    if (not exists(src)) :
        continue
    print(filename[:-6])
    new_name = 'AF-' + filename[:-6].upper() + '-F1-model_v2.pdb'
    shutil.copyfile(src, './all_pdb/' + new_name)
