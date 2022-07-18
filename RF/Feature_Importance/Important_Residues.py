import pandas as pd

with open('Residue_Location.csv', 'w') as f:
    df = pd.read_csv('../../data_files/TMdomains/TM.csv')
    print('accession,start,stop,sequence', file = f)
    for i in range(len(df)):
        accession = df.loc[i, 'protein']
        sequence = df.loc[i, 'TM6']
        start = int(df.loc[i, 's6']) + 18
        stop = int(df.loc[i, 's6']) + 27
        sequence = sequence[18:28]
        print(accession + ',' + str(start) + ',' + str(stop) + ',' + sequence, file = f)
