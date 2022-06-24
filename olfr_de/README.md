## Storing Data from the Matsunami Lab

This folder contains the entire repository of data from the Matsunami Lab. 

The csv files starting with pS6 contain data from gene expression experiments performed in mice exposed to the ligand listed in the file name and mice exposed to no odor. The change in gene expression for each receptor is described by logFC, and a p-value less than 0.05 indicates significant difference of expression, which is expected to correlate with whether the ligand has bound to the receptor.

This information is further stored in the uniprot_ligand_logfc.csv file, which matches protein and ligand pairs to their respective logFC values. 