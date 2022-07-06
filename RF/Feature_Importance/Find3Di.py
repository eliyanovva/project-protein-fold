#This script visualizes the most important structure elements in the proteins
#Script must be run while pwd is in the RF directory

import sys
sys.path.append('../../project-protein-fold/RF/')
import Globals
import os

Di_dict = Globals.initialize_3Di_dict()
indices = Globals.initialize_indices()

directory = 'pdb_data_files'

#Find the location of the 3Di kmer
def resinumber(protein, threeDi, domain): #domain is 1, 2, 3, or 4 corresponding to 3, 5, 6, or 7
    start = indices[protein][domain*2-2]
    seq = Di_dict[protein][domain-1]
    spot_in_seq = []
    current_find = -1
    while True:
        current_find = seq.find(threeDi, current_find+1)
        if current_find == -1:
            break
        spot_in_seq.append(current_find)
    ret = []
    for entry in spot_in_seq:
        ret.append(start + entry)
    return ret

print(resinumber('Q8VET4', 'LQ', 1))


"""
for filename in os.listdir(directory):
    f = os.path.join(directory, filename)
    if os.path.isfile(f):
        location_of_protein = directory + "/" + filename
        cmd.load(location_of_protein) #Load in pdb file
        name_of_protein = filename.replace("-model_v2.pdb", "")
        cmd.color("magenta", name_of_protein)

    """