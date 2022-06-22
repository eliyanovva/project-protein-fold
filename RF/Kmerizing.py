#from skbio import Sequence
"""
def featurize(seq,k,feat):
     s = Sequence(seq)
     dict = {}
     for kmer in s.iter_kmers(k, overlap=True):
          key = str(kmer)
          if key not in dict:
               dict[key] = 0
               feat.add(key)
          dict[key] += 1
     return dict
"""
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