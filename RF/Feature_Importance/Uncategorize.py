#This script must be interpreted by pymol and displays the top 5 AA sequences on the pdb file

#Amino Acid Codes:
#a: {Ala, Gly, Val} => {A, G, V}
#b: {Ile, Leu, Phe, Pro} => {I, L, F, P}
#c: {Tyr, Met, Thr, Ser} => {Y, M, T, S}
#d: {His, Asn, Gln, Trp} => {H, N, Q, W}
#e: {Arg, Lys} => {R, K}
#f: {Asp, Glu} => {D, E}
#g: {Cys} => {C}

#Imports
import os

#Method to determine all possible combinations of amino acids for each code
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

directory = '../../../../data_files/pdb_data_files'

#Iterate through each pdb
for filename in os.listdir(directory):
    f = os.path.join(directory, filename)
    if os.path.isfile(f):
        location_of_protein = directory + "/" + filename
        cmd.load(location_of_protein) #Load in pdb file
        name_of_protein = filename.replace("-model_v2.pdb", "")
        cmd.hide("everything")
        cmd.show("cartoon", name_of_protein)
        cmd.color("magenta", name_of_protein)
        cmd.set("cartoon_transparency",  "0.75", name_of_protein)

        #Select important Amino Acids
        printstatement1 = ""
        i = 0
        for sequence in uncat('eccbc'):
            if i > 0:
                printstatement1 += ' or '
            printstatement1 += 'pepseq ' + sequence
            i = 1
        print(printstatement1)

        printstatement2 = ""
        i = 0
        for sequence in uncat('acecc'):
            if i > 0:
                printstatement2 += ' or '
            printstatement2 += 'pepseq ' + sequence
            i = 1
        print(printstatement2)

        printstatement3 = ""
        i = 0
        for sequence in uncat('gcace'):
            if i > 0:
                printstatement3 += ' or '
            printstatement3 += 'pepseq ' + sequence
            i = 1
        print(printstatement3)

        printstatement4 = ""
        i = 0
        for sequence in uncat('cacec'):
            if i > 0:
                printstatement4 += ' or '
            printstatement4 += 'pepseq ' + sequence
            i = 1
        print(printstatement4)

        printstatement5 = ""
        i = 0
        for sequence in uncat('faacd'):
            if i > 0:
                printstatement5 += ' or '
            printstatement5 += 'pepseq ' + sequence
            i = 1
        print(printstatement5)

        printstatement6 = ""
        i = 0
        for sequence in uncat('bdbbd'):
            if i > 0:
                printstatement6 += ' or '
            printstatement6 += 'pepseq ' + sequence
            i = 1
        print(printstatement6)

        printstatement7 = ""
        i = 0
        for sequence in uncat('aacfc'):
            if i > 0:
                printstatement7 += ' or '
            printstatement7 += 'pepseq ' + sequence
            i = 1
        print(printstatement7)

        printstatement8 = ""
        i = 0
        for sequence in uncat('dbbdc'):
            if i > 0:
                printstatement8 += ' or '
            printstatement8 += 'pepseq ' + sequence
            i = 1
        print(printstatement8)

        printstatement9 = ""
        i = 0
        for sequence in uncat('bdaaa'):
            if i > 0:
                printstatement9 += ' or '
            printstatement9 += 'pepseq ' + sequence
            i = 1
        print(printstatement9)

        printstatement10 = ""
        i = 0
        for sequence in uncat('daaaa'):
            if i > 0:
                printstatement10 += ' or '
            printstatement10 += 'pepseq ' + sequence
            i = 1
        print(printstatement10)
        

        TM3 = "Res 70-170"
        TM5 = "Res 150-280"
        TM6 = "Res 210-310"
        TM7 = "Res 240-350"
        
        #Color important Amino Acids
        cmd.select("AA1", printstatement1 + " and " + TM6)
        cmd.create("obj1","AA1")
        cmd.color("red", "obj1")
        cmd.set("cartoon_transparency",  "0", "obj1")

        cmd.select("AA2", printstatement2 + " and " + TM6)
        cmd.create("obj2","AA2")
        cmd.color("orange", "obj2")
        cmd.set("cartoon_transparency",  "0", "obj2")

        cmd.select("AA3", printstatement3 + " and " + TM6)
        cmd.create("obj3","AA3")
        cmd.color("yellow", "obj3")
        cmd.set("cartoon_transparency",  "0", "obj3")

        cmd.select("AA4", printstatement4 + " and " + TM6)
        cmd.create("obj4","AA4")
        cmd.color("green", "obj4")
        cmd.set("cartoon_transparency",  "0", "obj4")

        cmd.select("AA5", printstatement5 + " and " + TM7)
        cmd.create("obj5","AA5")
        cmd.color("blue", "obj5")
        cmd.set("cartoon_transparency",  "0", "obj5")

        cmd.select("AA6", printstatement6 + " and " + TM6)
        cmd.create("obj6","AA6")
        cmd.color("pink", "obj6")
        cmd.set("cartoon_transparency",  "0", "obj6")

        cmd.select("AA7", printstatement7 + " and " + TM5)
        cmd.create("obj7","AA7")
        cmd.color("cyan", "obj7")
        cmd.set("cartoon_transparency",  "0", "obj7")

        cmd.select("AA8", printstatement8 + " and " + TM6)
        cmd.create("obj8","AA8")
        cmd.color("white", "obj8")
        cmd.set("cartoon_transparency",  "0", "obj8")

        cmd.select("AA9", printstatement9 + " and " + TM5)
        cmd.create("obj9","AA9")
        cmd.color("brown", "obj9")
        cmd.set("cartoon_transparency",  "0", "obj9")

        cmd.select("AA10", printstatement10 + " and " + TM5)
        cmd.create("obj10","AA10")
        cmd.color("gray40", "obj10")
        cmd.set("cartoon_transparency",  "0", "obj10")

        #Save as png
        cmd.deselect()
        savelocation = "../../Images/" + name_of_protein + ".png"
        cmd.png(savelocation)
        
        #Hide protein
        cmd.hide("everything", name_of_protein)