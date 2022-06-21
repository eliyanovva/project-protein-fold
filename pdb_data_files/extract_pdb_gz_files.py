"""This short script extracts the uniprot ids of the genes used in the olfr_de folder by
Hiro's lab, and copies locally the datafiles from the alphafold output folder with the 
corresponding uniprot id. They are archived.
"""

import shutil
import os

all_pdb_file_names = os.listdir('/usr/project/csplus2/dietrich/datafiles/mouse')

with open('./uniprot-ensemble-map.tab', 'r') as read_file:
    lines = read_file.readlines()
    for line in lines:
        line_list = line.split()
        uniprot_id = line_list[1]
        for pdb_file_name in all_pdb_file_names:
            if pdb_file_name == 'AF-' + uniprot_id + '-F1-model_v2.fasta':
                #print(pdb_file_name)
                shutil.copy2('/usr/project/csplus2/dietrich/datafiles/mouse/' + pdb_file_name, '/home/users/tep18/new_ppp/project-protein-fold/AminoAcidSequences/new_mapping_mouse')
            #print(pdb_file_name)
            #break

