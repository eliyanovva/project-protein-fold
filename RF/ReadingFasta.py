# Import featurize function
from Kmerizing import *

#Creating sequence class
class Seq:
    def __init__(self, name, sequence, dictionary):
        self.name = name
        self.sequence = sequence
        self.dictionary = dictionary
    def __repr__(self):
        return self.name

#Initialize sequence array
seqs = []

#Initialize Set of Features
Globals.initialize()

#Read fasta file
fasta = open("/home/users/sml96/bin/ProteinFoldRF/outputDb_ss.fasta")
i = 0
j = 0
for line in fasta:
    if j%2 == 0:
        name = line.replace('\n','')
        seqs.append(Seq(name, '', featurize('')))
        i += 1
    else:
        sequence = line.replace('\n','')
        seqs[i-1].sequence = sequence
        seqs[i-1].dictionary = featurize(sequence)
    j += 1

#This prints all of the k-mers identified in the sequences
#print(Globals.features)

print(len(Globals.features))
#TODO: For each sequence, assign a value to each feature
