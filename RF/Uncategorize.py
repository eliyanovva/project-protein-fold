#This script must be interpreted by pymol and displays the top 5 AA sequences on the pdb file

#Amino Acid Codes:
#a: {Ala, Gly, Val} => {A, G, V}
#b: {Ile, Leu, Phe, Pro} => {I, L, F, P}
#c: {Tyr, Met, Thr, Ser} => {Y, M, T, S}
#d: {His, Asn, Gln, Trp} => {H, N, Q, W}
#e: {Arg, Lys} => {R, K}
#f: {Asp, Glu} => {D, E}
#g: {Cys} => {C}

cmd.load("AF-A0A1L1SQ02-F1-model_v2.pdb") #Load in pdb file
cmd.color("magenta", "AF-A0A1L1SQ02-F1")

def uncat(feature):
    ret = ['']
    for letter in feature:
        length = len(ret)
        if letter == 'a':
            for i in range(0, length):
                ret.append(ret[i] + 'G')
                ret.append(ret[i] + 'V')
                ret[i] += 'A'
        if letter == 'b':
            for i in range(0, length):
                ret.append(ret[i] + 'L')
                ret.append(ret[i] + 'F')
                ret.append(ret[i] + 'P')
                ret[i] += 'I'
        if letter == 'c':
            for i in range(0, length):
                ret.append(ret[i] + 'M')
                ret.append(ret[i] + 'T')
                ret.append(ret[i] + 'S')
                ret[i] += 'Y'
        if letter == 'd':
            for i in range(0, length):
                ret.append(ret[i] + 'W')
                ret.append(ret[i] + 'N')
                ret.append(ret[i] + 'Q')
                ret[i] += 'H'
        if letter == 'e':
            for i in range(0, length):
                ret.append(ret[i] + 'K')
                ret[i] += 'R'
        if letter == 'f':
            for i in range(0, length):
                ret.append(ret[i] + 'E')
                ret[i] += 'D'
        if letter == 'g':
            for i in range(0, length):
                ret[i] += 'C'
    return ret


printstatement1 = ""
i = 0
for sequence in uncat('cgbbb'):
    if i > 0:
        printstatement1 += ' or '
    printstatement1 += 'pepseq ' + sequence
    i = 1
print(printstatement1)

printstatement2 = ""
i = 0
for sequence in uncat('ebcae'):
    if i > 0:
        printstatement2 += ' or '
    printstatement2 += 'pepseq ' + sequence
    i = 1
print(printstatement2)

printstatement3 = ""
i = 0
for sequence in uncat('gaecd'):
    if i > 0:
        printstatement3 += ' or '
    printstatement3 += 'pepseq ' + sequence
    i = 1
print(printstatement3)

printstatement4 = ""
i = 0
for sequence in uncat('cfggb'):
    if i > 0:
        printstatement4 += ' or '
    printstatement4 += 'pepseq ' + sequence
    i = 1
print(printstatement4)

printstatement5 = ""
i = 0
for sequence in uncat('dabfd'):
    if i > 0:
        printstatement5 += ' or '
    printstatement5 += 'pepseq ' + sequence
    i = 1
print(printstatement5)


cmd.select("AA1", printstatement1)
cmd.color("red", "AA1")

cmd.select("AA2", printstatement2)
cmd.color("orange", "AA2")

cmd.select("AA3", printstatement3)
cmd.color("yellow", "AA3")

cmd.select("AA4", printstatement4)
cmd.color("green", "AA4")

cmd.select("AA5", printstatement5)
cmd.color("blue", "AA5")