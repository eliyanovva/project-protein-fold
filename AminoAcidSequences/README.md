This folder contains amino acid sequences for all proteins in the dataset. 

The folder, "new_mapping_mouse", contains separate fasta files for each protein, and these files are concatenated in allsequences.fasta. 

The script, "categorize_proteins.py", reduces the dimensions of the data by converting each amino acid into a letter corresponding to 1 of 7 categories. These categories correspond to the categories used in ["String-Based Models for Predicting RNA-Protein Interaction"](https://dl.acm.org/doi/10.1145/3107411.3107508). 

This script was run on the "new_mapping_mouse" folder, using the "new_accessions.txt" file, to create the "new_mapping_mouse_categorized" folder, which contains these new categorized sequences. These files are concatenated in "fully_categorized.fasta", which is used as input to our algorithm. 

Donald Adjeroh, Maen Allaga, Jun Tan, Jie Lin, Yue Jiang, Ahmed Abbasi, and Xiaobo Zhou. 2017. String-Based Models for Predicting RNA-Protein Interaction. In Proceedings of the 8th ACM International Conference on Bioinformatics, Computational Biology,and Health Informatics (ACM-BCB '17). Association for Computing Machinery, New York, NY, USA, 661–666. https://doi.org/10.1145/3107411.3107508

