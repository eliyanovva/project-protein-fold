#!/bin/bash

for file in /home/users/tep18/new_ppp/project-protein-fold/pdb_data_files/*.pdb;
    do obabel $file -O ${file%.*}.bgf; 
done