#Utilizes the dictionary, "mouse_genes_ensembl_uniprot.json" (accession:ENSEMBL) to create a csv 
#which matches accession number and ligand to logFC. This is stored in the file "uniprot_ligand_logfc.csv".

import logging as log
#from .. import constants

#FIXME: REMOVE THE UGLY PATH AND FIX THE CONSTANTS
with open('/home/users/tep18/new_ppp/project-protein-fold/graph_cnn/data_prep/pdb_adjacency_data/file_names.txt', 'r') as protein_txt_list:
    protein_list = protein_txt_list.readlines()
    for i in range(len(protein_list)):
        protein_list[i] = protein_list[i][:-13]
    #    print(protein_list[i])


with open('/home/users/tep18/new_ppp/project-protein-fold/graph_cnn/data_prep/mol_adjacency_data/filenames.txt', 'r') as ligand_txt_list:
    ligand_list = ligand_txt_list.readlines()
    for i in range(len(ligand_list)):
        left_index = ligand_list[i].find('_')
        if left_index > 0:
            ligand_list[i] = ligand_list[i][left_index + 1: -13]
            

with open ('uniprot-ensemble-map.tab', 'r') as fp:
    lines = fp.readlines()[1:]
    GENE_DICT = {}
    for line in lines:
        line_list = line.split()
        GENE_DICT[line_list[0]] = line_list[1]


with open('file_names.txt', 'r') as filenames:
    file_names_lines = filenames.readlines()
    for filename in file_names_lines:
        left_index = filename.rfind('_')
        print(filename[left_index + 1 : -5])
        if filename[left_index + 1: -5] in ligand_list:
            with open(filename[:-1], 'r') as fp: # ignore \n chars in the end
                with open('uniprot_ligand_logfc_pvalue.csv', 'a') as write_file:    
                    file_lines = fp.readlines()
                    for line in file_lines[1:]:
                        line_list = line.split(',')
                        
                        ensembl_gene_id = line_list[1]
                        logFC = line_list[2]
                        pValue = float(line_list[5])
                        FDR = float(line_list[6])

                        if FDR <= 0.06 and pValue <= 0.06:
                            try:
                                uniprot_ID = GENE_DICT[ensembl_gene_id[1:-1]]    
                                if uniprot_ID in protein_list:
                                    csv_line = str(uniprot_ID) + ',' + filename[:-5] + ',' + str(logFC) + ',' + str(pValue) + '\n'
                                    write_file.write(csv_line)
                            except:
                                log.info(ensembl_gene_id + ' is not in the dictionary')


