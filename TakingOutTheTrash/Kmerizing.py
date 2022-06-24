#This script contains a function to extract kmers and kmer-frequency from protein sequences

def featurize(seq,k,feat):
     dict = {}
     for i in range(0, len(seq) - k + 1):
          kmer = ""
          for j in range(k):
               kmer += seq[i + j]
          if kmer not in dict:
               dict[kmer] = 0
               feat.add(kmer)
          dict[kmer] += 1
     return dict