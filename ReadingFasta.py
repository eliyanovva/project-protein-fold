# Import featurize function
from Kmerizing import *

fasta = open("./outputDb_ss.fasta")
for line in fasta:
    line = line.replace('\n','')
    if line[0] != '>':
        featurize(line) # At the moment, the code prints different dictionaries for each fasta sequence. 
        #TODO: Retain different dictionaries for each sequence in separate variables that contain the same keyset.
        # I will probably need to create a class so I can have sequence objects with dictionaries that get updated by the function.
        # The challenge will be keeping the keyset as small as it can be while making sure all of the sequences can use the same keyset.
        