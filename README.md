# project-protein-fold
The repository contains three different machine learning models which predict olfactory protein - ligand binding without docking them first. The models included are a Random Forest Model, a Convolutional Neural Network, and a Graph Convolutional Neural Network. All models were trained with data from wet lab experiments performed by the [Matsunami Lab](https://mgm.duke.edu/matsunami-lab) at Duke University. 


https://user-images.githubusercontent.com/77795920/181356910-98c6a2b1-53e1-4a01-bd44-2685e93a39dd.mp4


## ProteinFoldRF
Creating a String-Based Random Forest Model Informed by Tertiary Structure of Proteins and Ligands to Predict Binding.

The model input includes amino acid sequence, protein tertiary structure (represented by 3Di sequence), and ligand sequence and structure (represented by SMILES string) and output is a classification prediction as to whether the protein and ligand will bind or not, informed by logFC and p-value. 

## ProteinFold CNN

## ProteinFold GCN
This model has been trained specifically on mouse olfactory receptor data, and outputs the binding coefficient between proteins and ligands based on experimentally measured logFC score. The model uses a multi-input graph neural network, which represents both the protein and the ligand as a graph with an adjacency and a feature matrix. 

## References

Bioinformatics, 34(21), 2018, 3666–3674, doi: 10.1093/bioinformatics/bty374, Advance Access Publication Date: 10 May 2018, Original Paper

Claire A. de March, Yiqun Yu, Mengjue J. Ni, Kaylin A. Adipietro, Hiroaki Matsunami, Minghong Ma, and Jérôme Golebiowski. Conserved Residues Control Activation of Mammalian G Protein-Coupled Odorant Receptors. Journal of the American Chemical Society 2015 137 (26), 8611-8616, DOI: 10.1021/jacs.5b04659

Donald Adjeroh, Maen Allaga, Jun Tan, Jie Lin, Yue Jiang, Ahmed Abbasi, and Xiaobo Zhou. 2017. String-Based Models for Predicting RNA-Protein Interaction. In Proceedings of the 8th ACM International Conference on Bioinformatics, Computational Biology, and Health Informatics (ACM-BCB '17). Association for Computing Machinery, New York, NY, USA, 661-  666. https://doi.org/10.1145/3107411.3107508

Jimenez RC, Casajuana-Martin N, García-Recio A, Alcántara L, Pardo L, Campillo M, Gonzalez A. The mutational landscape of human olfactory G protein-coupled receptors. BMC Biol. 2021 Feb 5;19(1):21. doi: 10.1186/s12915-021-00962-0. PMID: 33546694; PMCID: PMC7866472.

Jumper, J., Evans, R., Pritzel, A. et al. Highly accurate protein structure prediction with AlphaFold. Nature 596, 583–589 (2021). https://doi.org/10.1038/s41586-021-03819-2

KDEEP: Protein–Ligand Absolute Binding Affinity Prediction via 3D-Convolutional Neural Networks
José Jiménez, Miha Škalič, Gerard Martínez-Rosell, and Gianni De Fabritiis
Journal of Chemical Information and Modeling 2018 58 (2), 287-296
DOI: 10.1021/acs.jcim.7b00650

Kim, S., Chen, J., Cheng, T., Gindulyte, A., He, J., He, S., Li, Q., Shoemaker, B. A., Thiessen, P. A., Yu, B., Zaslavsky, L., Zhang, J., & Bolton, E. E. (2019). PubChem in 2021: new data content and improved web interfaces. Nucleic Acids Res., 49(D1), D1388–D1395. https://doi.org/10.1093/nar/gkaa971

Large-Scale G Protein-Coupled Olfactory Receptor–Ligand Pairing
Xiaojing Cong, Wenwen Ren, Jody Pacalon, Rui Xu, Lun Xu, Xuewen Li, Claire A. de March, Hiroaki Matsunami, Hongmeng Yu, Yiqun Yu, and Jérôme Golebiowski
ACS Central Science 2022 8 (3), 379-387
DOI: 10.1021/acscentsci.1c01495

LeNail, (2019). NN-SVG: Publication-Ready Neural Network Architecture Schematics.
Journal of Open Source Software, 4(33), 747, https://doi.org/10.21105/joss.00747

National Center for Biotechnology Information (2022). PubChem Compound Summary for CID 7793, (-)-Citronellol. Retrieved July 26, 2022 from https://pubchem.ncbi.nlm.nih.gov/compound/citronellol.

Regalis, O. (2015). Beta-2-adrenergic-receptor. Wikimedia. Wikimedia. Retrieved July 25, 2022, from https://commons.wikimedia.org/wiki/File:Beta-2-adrenergic-receptor.png. 

Son J, Kim D (2021) Development of a graph convolutional neural network model for efficient prediction of protein-ligand binding affinities. PLoS ONE 16(4): e0249404. https://doi.org/10.1371/journal.pone.0249404 

Sriram K, Insel PA. G Protein-Coupled Receptors as Targets for Approved Drugs: How Many Targets and How Many Drugs? Mol Pharmacol. 2018 Apr;93(4):251-258. doi: 10.1124/mol.117.111062. Epub 2018 Jan 3. PMID: 29298813; PMCID: PMC5820538.

Stefan Doerr, Matthew J. Harvey, Frank Noé, and Gianni De Fabritiis. HTMD: High-throughput molecular dynamics for molecular discovery. Journal of Chemical Theory and Computation, 2016, 12 (4), pp 1845–1852. doi:10.1021/acs.jctc.6b00049

The UniProt Consortium, UniProt: the universal protein knowledgebase in 2021, Nucleic Acids Research, Volume 49, Issue D1, 8 January 2021, Pages D480–D489, https://doi.org/10.1093/nar/gkaa1100

UCSF ChimeraX: Structure visualization for researchers, educators, and developers. Pettersen EF, Goddard TD, Huang CC, Meng EC, Couch GS, Croll TI, Morris JH, Ferrin TE. Protein Sci. 2021 Jan;30(1):70-82

van Kempen M, Kim S, Tumescheit C, Mirdita M, Söding J, and Steinegger M. Foldseek: fast and accurate protein structure search. bioRxiv, doi:10.1101/2022.02.07.479398 (2022)

Wikipedia contributors. (2022, January 21). P-Cresol. In Wikipedia, The Free Encyclopedia. Retrieved 15:46, July 25, 2022, from https://en.wikipedia.org/w/index.php?title=P-Cresol&oldid=1067043032

Wikipedia contributors. (2022, July 24). Androstenone. In Wikipedia, The Free Encyclopedia. Retrieved 15:47, July 25, 2022, from https://en.wikipedia.org/w/index.php?title=Androstenone&oldid=1100049715 

Wikipedia contributors. (2022, July 23). Carvone. In Wikipedia, The Free Encyclopedia. Retrieved 15:47, July 25, 2022, from https://en.wikipedia.org/w/index.php?title=Carvone&oldid=1099965674

Yasi EA, Eisen SL, Wang H, Sugianto W, Minniefield AR, Hoover KA, Branham PJ, Peralta-Yahya P. Rapid Deorphanization of Human Olfactory Receptors in Yeast. Biochemistry. 2019 Apr 23;58(16):2160-2166. doi: 10.1021/acs.biochem.8b01208. Epub 2019 Apr 12. PMID: 30977365; PMCID: PMC6482435.

Yuan, S., Dahoun, T., Brugarolas, M. et al. Computational modeling of the olfactory receptor Olfr73 suggests a molecular basis for low potency of olfactory receptor-activating compounds. Commun Biol 2, 141 (2019). https://doi.org/10.1038/s42003-019-0384-8]
