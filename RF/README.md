# Random Forest Model

This folder contains all scripts for processing the amino acid sequences, 3Di sequences, and SMILES strings into the input matrix for the random forest model. It also contains scripts to create the vector of classification values for the model. Finally, it contains the Random Forest Algorithm, Optimization scripts, and Metrics functions. 

Our Random Forest model balances the data through Oversampling and optimizes hyperparameters through cross validation. 

At the moment, the algorithm is performing at an acceptable level with impressive speed. The ROC AUC is consistently around 0.8. 

WARNING: In the FixedClassificationModel.py and FeatureImportance.py files, there is an option to edit the number of processors used in the undersampling method. Reduce processors if running on a personal machine. 