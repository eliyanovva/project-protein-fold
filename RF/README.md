# SRF: String-Based Random Forest Model

SRF is a novel machine learning architecture for the fast, accurate, and interpretable prediction of pairings between Olfactory Receptors and odors. SRF formats protein and ligand structure and sequence as strings, with features input into the model being k-mers of these strings. SRF is unique in its ability to account for protein structure without requiring the co-crystal structure of the protein and ligand bound together.

SRF has achieved significant success with an average Receiver Operator Characteristic - Area Under the Curve (ROC-AUC) of 
0.9975, average Precision-Recall AUC of 0.9995, and average Balanced Accuracy of 0.9459.

This folder contains functions to train and test the model, predict protein-ligand pairs, and analyze the important features associated with protein-ligand binding. 