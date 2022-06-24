# project-protein-fold
The repository contains three different machine learning models which predict olfactory protein - ligand binding without docking them first. The models included are a Random Forest Model, a Convolutional Neural Network, and a Graph Convolutional Neural Network. All models were trained with data from wet lab experiments performed by the [Matsunami Lab](https://mgm.duke.edu/matsunami-lab) at Duke University. 

## ProteinFoldRF
Creating a String-Based Random Forest Model Informed by Tertiary Structure of Proteins and Ligands to Predict Binding.

The model input includes amino acid sequence, protein tertiary structure (represented by 3Di sequence), and ligand sequence and structure (represented by SMILES string) and output is a classification prediction as to whether the protein and ligand will bind or not, informed by logFC and p-value. 

## ProteinFold CNN

## ProteinFold GCN

## References

Donald Adjeroh, Maen Allaga, Jun Tan, Jie Lin, Yue Jiang, Ahmed Abbasi, and Xiaobo Zhou. 2017. String-Based Models for Predicting RNA-Protein Interaction. In Proceedings of the 8th ACM International Conference on Bioinformatics, Computational Biology,and Health Informatics (ACM-BCB '17). Association for Computing Machinery, New York, NY, USA, 661–666. https://doi.org/10.1145/3107411.3107508

Kim, S., Chen, J., Cheng, T., Gindulyte, A., He, J., He, S., Li, Q., Shoemaker, B. A., Thiessen, P. A., Yu, B., Zaslavsky, L., Zhang, J., & Bolton, E. E. (2019). PubChem in 2021: new data content and improved web interfaces. Nucleic Acids Res., 49(D1), D1388–D1395. https://doi.org/10.1093/nar/gkaa971

van Kempen M, Kim S, Tumescheit C, Mirdita M, Söding J, and Steinegger M. Foldseek: fast and accurate protein structure search. bioRxiv, doi:10.1101/2022.02.07.479398 (2022)
