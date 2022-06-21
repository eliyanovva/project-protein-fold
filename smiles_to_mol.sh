#!/bin/bash

for file in /home/users/tep18/new_ppp/project-protein-fold/mol_files/*;
    do obabel $file -O ${file%.*}.mol --gen3D; 
done