#This Script visualizes the relative importance of ligand, protein structure, and protein sequence information

import math as m
import matplotlib.pyplot as plt

with open('sulfur_importance.txt') as f: 
    lines = f.readlines()
i=0
values = {'ligand':0, 'sequence':0, 'structure':0}
domains = {'TM3':0, 'TM5':0, 'TM6':0, 'TM7':0}
for line in lines:
    if i%2 == 0:
        if 'TM' in line:
            domain = line[-4:-1]
            if line[0].isupper():
                type = 'structure'
            else:
                type = 'sequence'
        else:
            domain = None
            type = 'ligand'
    else:
        values[type] += float(line)
        if domain != None:
            domains[domain] += float(line)
    i+=1

print(values)
print(domains)

fig,ax = plt.subplots()
plt.bar(list(values.keys()), list(values.values()))
ax.set_title("Feature Importance by Type")
plt.ylabel('Mean Decrease in Impurity')
fig.savefig('Relative_Feature_Importance_Sulfur.png')

fig2, ax2 = plt.subplots()
plt.bar(list(domains.keys()), list(domains.values()))
ax2.set_title("Feature Importance by Transmembrane Domain")
plt.ylabel('Mean Decrease in Impurity')
fig2.savefig('Domain_Importance_Sulfur.png')