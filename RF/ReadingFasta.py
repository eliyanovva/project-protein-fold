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

#Initialize Set of Features
Globals.initialize()

#Read fasta file
fasta = open("/home/users/sml96/bin/project-protein-fold/AminoAcidSequences/categorized.fasta")
i = 0
j = 0
for line in fasta:
    if line[0] == '>':
        name = line.replace('\n','')
        Globals.seqs.append(Seq(name, '', featurize('',5)))
        i += 1
    else:
        sequence = line.replace('\n','')
        Globals.seqs[i-1].sequence = Globals.seqs[i-1].sequence + sequence
        Globals.seqs[i-1].dictionary = featurize(Globals.seqs[i-1].sequence,5)
    j += 1

#This prints all of the k-mers identified in the sequences
#print(Globals.features)

#This prints the number of k-mers identified in all of the sequences
#print(len(Globals.features))

#This prints number of sequences
#print(len(Globals.seqs))

for seq in Globals.seqs:
    newseq = []
    for kmer in Globals.features:
        if kmer not in seq.dictionary:
               seq.dictionary[kmer] = 0
        newseq.append(seq.dictionary.get(kmer))
    Globals.matrix.append(newseq)

#Prints a matrix with columns corresponding to k-mers and rows corresponding to proteins
#print(Globals.matrix)

#Prints protein names
#print(Globals.seqs)
