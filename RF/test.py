import numpy as np
array = np.array([1,2,3])
while array.size > 0:
        index = np.argmax(array)
        print(array[index])
        array = np.delete(array, index)