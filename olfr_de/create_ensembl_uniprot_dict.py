import json

with open('MGIBatchReport_20220531_153700.txt', 'r') as reference_file:
    reference_lines = reference_file.readlines()
    dict = {}
    for line in reference_lines:
        line_list = line.split()
        dict[line_list[-2]] = line_list[-1]

with open('mouse_genes_ensembl_uniprot.json', 'w') as dict_file:
    json.dump(dict, dict_file)    
