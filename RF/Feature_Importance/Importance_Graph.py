#This Script visualizes the relative importance of ligand, protein structure, and protein sequence information

import math as m
import matplotlib.pyplot as plt

with open('important_features.txt') as f: 
    lines = f.readlines()
i=0
values = {'ligand':0, 'sequence':0, 'structure':0}
for line in lines:
    if i%2 == 0:
        if 'TM' in line:
            if line[0].isupper():
                type = 'structure'
            else:
                type = 'sequence'
        else:
            type = 'ligand'
    else:
        values[type] += float(line)
    i+=1

print(values)

fig = plt.figure()
plt.bar(list(values.keys()), list(values.values()))

fig.savefig('Relative_Feature_Importance.png')