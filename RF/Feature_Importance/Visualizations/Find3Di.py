#This script must be interpreted in pymol and visualizes the most important structure elements in the proteins
#Script must be run while pwd is in the RF directory

import sys
sys.path.append('../../../project-protein-fold/RF/')
import Globals
import CombineLigandsProteins
import os

plist = Globals.initialize_protein_list()
Di_dict = Globals.initialize_3Di_dict(list(plist))
indices = Globals.initialize_indices(list(plist))
print(indices)

directory = '../data_files/pdb_data_files'

#Find the location of the 3Di kmer
def resinumber(protein, threeDi, domain): #domain is 1, 2, 3, or 4 corresponding to 3, 5, 6, or 7
    print(indices[protein])
    print(domain*2-2)
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
        cmd.hide("everything")
        cmd.show("cartoon", name_of_protein)
        cmd.color("magenta", name_of_protein)
        cmd.set("cartoon_transparency",  "0.75", name_of_protein)
        
        name = name_of_protein.replace("AF-", "")
        name = name.replace("-F1", "")

        save = False 

        if name in Di_dict:
            #Visualize the 3Di features
            if len(resinumber(name,'VVLVV', 2)) > 0:
                printstatement1 = ""
                i = 0
                for structure in resinumber(name, 'VVLVV', 2):
                    if i > 0:
                        printstatement1 += ' or '
                    printstatement1 += 'Resi ' + str(structure) + '-' + str(structure+5)
                    i = 1
                cmd.select("Di1", printstatement1)
                cmd.create("obj1", "Di1")
                cmd.color("red", "obj1")
                cmd.set("cartoon_transparency",  "0", "obj1")
                save = True
            
            if len(resinumber(name,'VCCVV', 3)) > 0:
                printstatement2 = ""
                i = 0
                for structure in resinumber(name,'VCCVV', 3):
                    if i > 0:
                        printstatement2 += ' or '
                    printstatement2 += 'Resi ' + str(structure) + '-' + str(structure+5)
                    i = 1
                cmd.select("Di2", printstatement2)
                cmd.create("obj2","Di2")
                cmd.color("orange", "obj2")
                cmd.set("cartoon_transparency",  "0", "obj2")
                save = True
            
            if len(resinumber(name,'SSNNV', 4)) > 0:
                printstatement3 = ""
                i = 0
                for structure in resinumber(name,'SSNNV', 4):
                    if i > 0:
                        printstatement3 += ' or '
                    printstatement3 += 'Resi ' + str(structure) + '-' + str(structure+5)
                    i = 1
                cmd.select("Di3", printstatement3)
                cmd.create("obj3", "Di3")
                cmd.color("yellow", "obj3")
                cmd.set("cartoon_transparency",  "0", "obj3")
                save = True
            
            if len(resinumber(name,'VVSVC', 2)) > 0:
                printstatement4 = ""
                i = 0
                for structure in resinumber(name,'VVSVC', 2):
                    if i > 0:
                        printstatement4 += ' or '
                    printstatement4 += 'Resi ' + str(structure) + '-' + str(structure+5)
                    i = 1
                cmd.select("Di4", printstatement4)
                cmd.create("obj4", "Di4")
                cmd.color("green", "obj4")
                cmd.set("cartoon_transparency",  "0", "obj4")
                save = True

            if len(resinumber(name,'SPQVP', 4)) > 0:
                printstatement5 = ""
                i = 0
                for structure in resinumber(name,'SPQVP', 4):
                    if i > 0:
                        printstatement5 += ' or '
                    printstatement5 += 'Resi ' + str(structure) + '-' + str(structure+5)
                    i = 1
                cmd.select("Di5", printstatement5)
                cmd.create("obj5", "Di5")
                cmd.color("blue", "obj5")
                cmd.set("cartoon_transparency",  "0", "obj5")
                save = True

        #Save as png
        if save:
            cmd.deselect()
            savelocation = "Feature_Importance/Visualizations/Sulfur_Images/" + "Di_" + name_of_protein + ".png"
            cmd.png(savelocation)
        
        #Hide protein
        cmd.hide("everything", name_of_protein)
