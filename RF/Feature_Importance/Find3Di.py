#This script visualizes the most important structure elements in the proteins
#Script must be run while pwd is in the RF directory

import sys
sys.path.append('../../project-protein-fold/RF/')
import Globals
import os

Di_dict = Globals.initialize_3Di_dict()
indices = Globals.initialize_indices()

directory = '../data_files/pdb_data_files'

#Find the location of the 3Di kmer
def resinumber(protein, threeDi, domain): #domain is 1, 2, 3, or 4 corresponding to 3, 5, 6, or 7
    start = indices[protein][domain*2-2]
    seq = Di_dict[protein][domain-1]
    spot_in_seq = []
    current_find = -1
    while True:
        current_find = seq.find(threeDi, current_find+1)
        if current_find == -1:
            break
        spot_in_seq.append(current_find)
    ret = []
    for entry in spot_in_seq:
        ret.append(start + entry)
    return ret

for filename in os.listdir(directory):
    f = os.path.join(directory, filename)
    if os.path.isfile(f):
        location_of_protein = directory + "/" + filename
        cmd.load(location_of_protein) #Load in pdb file
        name_of_protein = filename.replace("-model_v2.pdb", "")
        cmd.color("magenta", name_of_protein)
        
        name = name_of_protein.replace("AF-", "")
        name = name.replace("-F1", "")

        if name in Di_dict:
            #Visualize the 3Di features
            printstatement1 = ""
            i = 0
            for structure in resinumber(name,'QVVCV', 3):
                if i > 0:
                    printstatement1 += ' or '
                printstatement1 += 'Resi ' + str(structure) + '- ' + str(structure+5)
                i = 1
            print(printstatement1)

            printstatement2 = ""
            i = 0
            for structure in resinumber(name,'PPNVS', 4):
                if i > 0:
                    printstatement2 += ' or '
                printstatement2 += 'Resi ' + str(structure) + '- ' + str(structure+5)
                i = 1
            print(printstatement2)

            printstatement3 = ""
            i = 0
            for structure in resinumber(name,'PNPSS', 4):
                if i > 0:
                    printstatement3 += ' or '
                printstatement3 += 'Resi ' + str(structure) + '- ' + str(structure+5)
                i = 1
            print(printstatement3)

            printstatement4 = ""
            i = 0
            for structure in resinumber(name,'VNPLV', 3):
                if i > 0:
                    printstatement4 += ' or '
                printstatement4 += 'Resi ' + str(structure) + '- ' + str(structure+5)
                i = 1
            print(printstatement4)

            printstatement5 = ""
            i = 0
            for structure in resinumber(name,'VVCCV', 4):
                if i > 0:
                    printstatement5 += ' or '
                printstatement5 += 'Resi ' + str(structure) + '- ' + str(structure+5)
                i = 1
            print(printstatement5)

            #Color important structures
            cmd.select("Di1", printstatement1)
            cmd.color("red", "Di1")

            cmd.select("Di2", printstatement2)
            cmd.color("orange", "Di2")

            cmd.select("Di3", printstatement3)
            cmd.color("yellow", "Di3")

            cmd.select("Di4", printstatement4)
            cmd.color("green", "Di4")

            cmd.select("Di5", printstatement5)
            cmd.color("blue", "Di5")

        #Save as png
        savelocation = "/Feature_Importance/Images/" + "Di_" + name_of_protein + ".png"
        cmd.png(savelocation)
        
        #Hide protein
        cmd.hide("everything", name_of_protein)
