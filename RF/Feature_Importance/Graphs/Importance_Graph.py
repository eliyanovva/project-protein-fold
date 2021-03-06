#This Script visualizes the relative importance of ligand, protein structure, and protein sequence information

import math as m
import matplotlib.pyplot as plt


#Sum the Mean Impurity Decrease for each category: ligand, sequence, structure and TM3, TM5, TM6, TM7
with open('../important_features.txt') as f: 
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

#Plot relative importance of ligand features, protein sequence, and protein structure
fig,ax = plt.subplots()
plt.bar(list(values.keys()), list(values.values()))
plt.xticks(size = 16)
plt.yticks(size = 12)
fig.suptitle("Feature Importance by Type", fontsize = 20)
plt.ylabel('Mean Decrease in Impurity', fontsize = 16)
fig.savefig('Relative_Feature_Importance.jpg', dpi = 400)

#Plot relative importance of the TM domains
fig2, ax2 = plt.subplots()
plt.bar(list(domains.keys()), list(domains.values()))
ax2.set_title("Feature Importance by Transmembrane Domain")
plt.ylabel('Mean Decrease in Impurity')
fig2.savefig('Domain_Importance.png')