#Imports
# Need iter_kmers function
from skbio import Sequence
# Need createDictionary function
from AllKmers import *

#Input 3Di sequence here
s = Sequence('DDDPLDPAWWFQEQADPDDLVVLVVVLVVLVVLLVLLCVLLVVLLVLCVPDPVCPDLLSLLSNLLSVLSNQLSVLQSVLLSCCSVPSPDIDTPVSQLVSQLSNQLSVLLNLLSLLVSLVLLLCCQVPVVCSCVVPDPVVSVVSSVVSNVVSNVLSVVLSVLLVPFRFDDRSYQHASGNHSVSSLVRTPDDCVVNVVSVVVSCCCRLVVSLVSSVVSVVSSVVSLVPDPDPVSNVVSCLSCVLSVVLSCLVNVLVCCLVPPDPPPDDPSVSNVSVVCSSSVSSSCNSVSVQVSDPVNVVSVVVVVVVVVVD')

#Initialize 3Di dictionary
set1 = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S' 'T', 'X']
k = 5
createDictionary(set1, k)
print(dict['DDDPL'])

for kmer in s.iter_kmers(5, overlap=True):
     key = str(kmer)
     dict[key] += 1
print(dict)
     
