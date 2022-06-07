# Python 3 program to print all
# possible strings of length k
     
# The method that prints all
# possible strings of length k.
# It is mainly a wrapper over
# recursive function printAllKLengthRec()
def printAllKLength(set, k):
 
    n = len(set)
    printAllKLengthRec(set, "", n, k)
 
# The main recursive method
# to print all possible
# strings of length k
def printAllKLengthRec(set, prefix, n, k):
     
    # Base case: k is 0,
    # print prefix
    if (k == 0) :
        print(prefix)
        return
 
    # One by one add all characters
    # from set and recursively
    # call for k equals to k-1
    for i in range(n):
    
            # Next character of input added
            newPrefix = prefix + set[i]
            
            # k is decreased, because
            # we have added a new character
            printAllKLengthRec(set, newPrefix, n, k - 1)
 
# Driver Code
if __name__ == "__main__":
     
    set1 = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S' 'T', 'X']
    k = 5
    printAllKLength(set1, k)
 
# This code is adapted from code by
# ChitraNayal