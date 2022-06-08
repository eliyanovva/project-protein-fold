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
fasta = open("/home/users/sml96/bin/project-protein-fold/AminoAcidSequences/AAsequences_categorized.fasta")
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

#This prints the number of k-mers identified in all of the sequences
#print(len(Globals.features))
#TODO: For each sequence, assign a value to each feature

#Create input matrix to RF model
matrix = []
for kmer in Globals.features:
    newkmer = []
    for seq in seqs:
        if kmer not in seq.dictionary:
               seq.dictionary[kmer] = 0
        newkmer.append(seq.dictionary.get(kmer))
    matrix.append(newkmer)

#Prints a matrix with columns corresponding to k-mers and rows corresponding to proteins
#print(matrix)

