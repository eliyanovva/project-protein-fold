with open('TM.txt') as f:
    lines = f.readlines()

with open('TM.csv', 'w') as f:
    printstatement = 'protein,TM3,s3,e3,TM5,s5,e5,TM6,s6,e6,TM7,s7,e7'
    for line in lines:
        if line[0] == '>':
            print('stop')
            print(printstatement)
            print(printstatement, file = f)
            print('start')
            line = line.replace('>', '')
            printstatement = line.replace('\n', '')
            print(printstatement)
        else:
            printstatement += ',' + line.replace('\n', '')
