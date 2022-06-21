# Import featurize function
from Kmerizing import *
import Globals
import numpy as np
import Filtering
import labels

#Creating sequence class
class Seq:
    def __init__(self, name, sequence, dictionary):
        self.name = name
        self.sequence = sequence
        self.dictionary = dictionary
    def __repr__(self):
        return self.name
    def __eq__(self, other):
        return self.name == other.name


#Initialize Set of Features
Globals.initialize()

def make_seqvar(fasta, seqvar, feat):
    i = 0
    j = 0
    for line in fasta:
        if line[0] == '>':
            name = line.replace('\n','')
            name = name.replace('>', '')
            seqvar.append(Seq(name, '', featurize('',7,feat)))
            i += 1
        else:
            sequence = line.replace('\n','')
            seqvar[i-1].sequence = seqvar[i-1].sequence + sequence
            seqvar[i-1].dictionary = featurize(seqvar[i-1].sequence,7,feat)
        j += 1

    return seqvar, feat

def makematrix(seqvar, feat, mat):
    for seq in seqvar:
        newseq = []
        for kmer in feat:
            if kmer not in seq.dictionary:
                seq.dictionary[kmer] = 0
            newseq.append(seq.dictionary.get(kmer))
        mat.append(np.array(newseq))
    return mat

"""
def makematrix(fasta, seqvar, feat, mat):
    i = 0
    j = 0
    for line in fasta:
        if line[0] == '>':
            name = line.replace('\n','')
            name = name.replace('>', '')
            seqvar.append(Seq(name, '', featurize('',7,feat)))
            i += 1
        else:
            sequence = line.replace('\n','')
            seqvar[i-1].sequence = seqvar[i-1].sequence + sequence
            seqvar[i-1].dictionary = featurize(seqvar[i-1].sequence,7,feat)
        j += 1

#seq.dictionary is the kmer frequency dict for each protein

    for seq in seqvar:
        newseq = []
        for kmer in feat:
            if kmer not in seq.dictionary:
                seq.dictionary[kmer] = 0
            newseq.append(seq.dictionary.get(kmer))
        mat.append(np.array(newseq))
"""

#Creating output for categorized amino acids
# Read fasta file
#fasta1 = open("/home/users/sml96/bin/project-protein-fold/AminoAcidSequences/categorized.fasta")
#Remove insignificant kmers

# Make the matrix
#makematrix(fasta1, Globals.categorized_seqs, Globals.categorized_features, Globals.categorized_matrix)
# View output
#print(Globals.categorized_seqs)
#print(Globals.categorized_features)
#print(Globals.categorized_matrix)

#Creating output for 3Di sequences
# Read fasta file
#fasta2 = open("/home/users/sml96/bin/project-protein-fold/foldseek-master/foldseek/foldseek/outputDb_ss.fasta")
# Make the matrix
#makematrix(fasta2, Globals.di_seqs, Globals.di_features, Globals.di_matrix)
# View output
#print(Globals.di_seqs)
#print(Globals.di_features)
#print(np.array(Globals.di_matrix))

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

#Check whether categorized and 3di sequences contain all the same proteins in the same order
"""i = 0
for element in Globals.categorized_seqs:
    if element != Globals.di_seqs[i]:
        print("not same")
        print(element)
        i+=1
    i +=1"""
