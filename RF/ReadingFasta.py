# Import featurize function
from Kmerizing import *
import numpy as np

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

        #newseqs = []
        #newseq_counts = list(np.repeat(0, len(seqvar)))

        for kmer in feat:
            if kmer not in seq.dictionary:
                seq.dictionary[kmer] = 0
            newseq.append(seq.dictionary.get(kmer))
        #if newseqs.count(newseq) == 0:
        #    newseqs.append(newseq)
        #ind = newseqs.index(newseq)
        #newseq_counts[ind] += 1
        mat.append(np.array(newseq))
    return mat
