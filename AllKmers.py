# Python 3 program to create dictionary of k-mers

# Wrapper over
# recursive function createDictionaryRec()
dict = {}

def createDictionary(set, k):
    n = len(set)
    dict = createDictionaryRec(set, "", n, k)
 
# The main recursive method
# to print all possible
# strings of length k
def createDictionaryRec(set, prefix, n, k):
    global dict
    # Base case: k is 0,
    # print prefix
    if (k == 0) :
        dict[prefix] = 0
        return 
 
    # One by one add all characters
    # from set and recursively
    # call for k equals to k-1
    for i in range(n):
    
            # Next character of input added
            newPrefix = prefix + set[i]
            
            # k is decreased, because
            # we have added a new character
            createDictionaryRec(set, newPrefix, n, k - 1)
 
# Driver Code
if __name__ == "__main__":
     
    set1 = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S' 'T', 'X']
    k = 5
    createDictionary(set1, k)
    print(dict['LDPAW'])
 
# This code is adapted from code by
# ChitraNayal