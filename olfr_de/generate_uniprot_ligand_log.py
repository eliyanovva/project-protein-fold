import json
import logging as log
import log_config

with open ('mouse_genes_ensembl_uniprot.json', 'r') as fp:
    GENE_DICT = json.load(fp)

with open('file_names.txt', 'r') as filenames:
    file_names_lines = filenames.readlines()
    for filename in file_names_lines:
        with open(filename[:-1], 'r') as fp:
            with open('uniprot_ligand_logfc.csv', 'a') as write_file:    
                file_lines = fp.readlines()
                for line in file_lines[1:]:
                    line_list = line.split(',')
                    ensembl_gene_id = line_list[1]
                    logFC = line_list[2]
                    try:
                        uniprot_ID = GENE_DICT[ensembl_gene_id[1:-1]]    
                        csv_line = str(uniprot_ID) + ',' + filename + ',' + str(logFC) + '\n'
                        write_file.write(csv_line)
                    except:
                        log.debug(ensembl_gene_id + ' is not in the dictionary')
