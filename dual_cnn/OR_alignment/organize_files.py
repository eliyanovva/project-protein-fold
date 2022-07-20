import shutil
import os
import re
import pickle

print(os.chdir('/home/users/jvs15/project-protein-fold/dual_cnn/OR_alignment'))
pickle_file = open('existing_symb_to_unip.p', 'rb')

# Keys = gene name | Values = UNIPROT
d = pickle.load(pickle_file)

