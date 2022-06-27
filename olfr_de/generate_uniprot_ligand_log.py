#Utilizes the dictionary, "mouse_genes_ensembl_uniprot.json" (accession:ENSEMBL) to create a csv 
#which matches accession number and ligand to logFC. This is stored in the file "uniprot_ligand_logfc.csv".

import logging as log
#from .. import constants


with open ('uniprot-ensemble-map.tab', 'r') as fp:
    lines = fp.readlines()[1:]
    GENE_DICT = {}
    for line in lines:
        line_list = line.split()
        GENE_DICT[line_list[0]] = line_list[1]


with open('file_names.txt', 'r') as filenames:
    file_names_lines = filenames.readlines()
    for filename in file_names_lines:
        with open(filename[:-1], 'r') as fp: # ignore \n chars in the end
            with open('uniprot_ligand_logfc_pvalue.csv', 'a') as write_file:    
                file_lines = fp.readlines()
                for line in file_lines[1:]:
                    line_list = line.split(',')
                    
                    ensembl_gene_id = line_list[1]
                    logFC = line_list[2]
                    pValue = float(line_list[5])

                    if pValue <= 0.1:
                        try:
                            uniprot_ID = GENE_DICT[ensembl_gene_id[1:-1]]    
                            csv_line = str(uniprot_ID) + ',' + filename[:-5] + ',' + str(logFC) + ',' + str(pValue) + '\n'
                            write_file.write(csv_line)
                        except:
                            log.info(ensembl_gene_id + ' is not in the dictionary')
