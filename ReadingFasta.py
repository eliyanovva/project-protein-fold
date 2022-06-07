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
    j +=1
    

        # At the moment, the code prints different dictionaries for each fasta sequence. 
        #TODO: Retain different dictionaries for each sequence in separate variables that contain the same keyset.
        # I will probably need to create a class so I can have sequence objects with dictionaries that get updated by the function.
        # The challenge will be keeping the keyset as small as it can be while making sure all of the sequences can use the same keyset.
 
