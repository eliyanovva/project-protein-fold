from skbio import Sequence

def featurize(seq):
     s = Sequence(seq)
     dict = {}
     for kmer in s.iter_kmers(5, overlap=True):
          key = str(kmer)
          if key not in dict:
               dict[key] = 0
          dict[key] += 1
     return dict
     
