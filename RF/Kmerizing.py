from skbio import Sequence
import Globals

Globals.initialize()

def featurize(seq,k):
     s = Sequence(seq)
     dict = {}
     for kmer in s.iter_kmers(k, overlap=True):
          key = str(kmer)
          if key not in dict:
               dict[key] = 0
               Globals.features.add(key)
          dict[key] += 1
     return dict

