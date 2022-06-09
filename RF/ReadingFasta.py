# Import featurize function
from Kmerizing import *
import Globals

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

def makematrix(fasta, seqvar, feat, mat):
    i = 0
    j = 0
    for line in fasta:
        if line[0] == '>':
            name = line.replace('\n','')
            seqvar.append(Seq(name, '', featurize('',5,feat)))
            i += 1
        else:
            sequence = line.replace('\n','')
            seqvar[i-1].sequence = seqvar[i-1].sequence + sequence
            seqvar[i-1].dictionary = featurize(seqvar[i-1].sequence,5,feat)
        j += 1

    #This prints all of the k-mers identified in the sequences
    #print(Globals.categorized_features)

    #This prints the number of k-mers identified in all of the sequences
    #print(len(Globals.categorized_features))

    #This prints number of sequences
    #print(len(Globals.seq))

    for seq in seqvar:
        newseq = []
        for kmer in feat:
            if kmer not in seq.dictionary:
                seq.dictionary[kmer] = 0
            newseq.append(seq.dictionary.get(kmer))
        mat.append(newseq)

    #Prints a matrix with columns corresponding to k-mers and rows corresponding to proteins
    #print(Globals.categorized_matrix)

    #Prints protein names
    #print(Globals.categorized_seqs)

#Creating output for categorized amino acids
# Read fasta file
fasta1 = open("/home/users/sml96/bin/project-protein-fold/AminoAcidSequences/categorized.fasta")
# Make the matrix
makematrix(fasta1, Globals.categorized_seqs, Globals.categorized_features, Globals.categorized_matrix)
# View output
#print(Globals.categorized_seqs)
#print(Globals.categorized_features)
#print(Globals.categorized_matrix)

#Creating output for 3Di sequences
# Read fasta file
fasta2 = open("/home/users/sml96/bin/project-protein-fold/AminoAcidSequences/mouse.fasta") #TODO insert path to 3Di sequences here
# Make the matrix
makematrix(fasta2, Globals.di_seqs, Globals.di_features, Globals.di_matrix)
# View output
#print(Globals.di_seqs)
#print(Globals.di_features)
#print(Globals.di_matrix)

def import_variables():
    global sequence_seqs
    sequence_seqs = Globals.categorized_seqs
    global sequence_features 
    sequence_features = Globals.categorized_features
    global sequence_matrix
    sequence_matrix = Globals.categorized_matrix
    global structure_seqs
    structure_seqs = Globals.di_seqs
    global structure_features 
    structure_features = Globals.di_features
    global structure_matrix
    structure_matrix = Globals.di_matrix


