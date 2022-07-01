#This script must be interpreted by pymol and displays a list of sequences as spheres on the pdb file

#Amino Acid Codes:
#a: {Ala, Gly, Val} => {A, G, V}
#b: {Ile, Leu, Phe, Pro} => {I, L, F, P}
#c: {Tyr, Met, Thr, Ser} => {Y, M, T, S}
#d: {His, Asn, Gln, Trp} => {H, N, Q, W}
#e: {Arg, Lys} => {R, K}
#f: {Asp, Glu} => {D, E}
#g: {Cys} => {C}

cmd.load("AF-A0A1L1SQ02-F1-model_v2.pdb") #Load in pdb file

for entry in ['ccbbc', 'aaeae', 'aaege', 'cbcbd', 'ababa']:
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

    printstatement = ""
    i = 0
    for sequence in uncat(entry):
        if i > 0:
            printstatement += ' or '
        printstatement += 'pepseq ' + sequence
        i = 1
    print(printstatement)

    cmd.select("AAs", printstatement)

    cmd.show("spheres", "AAs")