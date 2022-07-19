import pandas as pd

with open('Residue_Location.csv', 'w') as f:
    print('gene,start,sequence', file = f)
    with open('Aligned_Location.txt') as w:
        lines = w.readlines()
        for line in lines:
            line = line.replace('\n', '')
            line = line[0:9] + ',' + line[9:12] + ',' + line[12:]
            line = line.replace(" ", "")
            print(line, file = f)
