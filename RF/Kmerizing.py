from skbio import Sequence

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

