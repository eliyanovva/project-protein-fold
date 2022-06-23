#Code categorizes amino acids into 7 different categories to reduce dimensionality

#Open file containing the protein accession numbers
fr = open('new_accessions.txt', "r")
protein_list = fr.readlines()
for i in range(len(protein_list)):
    protein_list[i] = protein_list[i][:-1]

#Open each file containing Amino Acid sequence
for id in protein_list:
    fr = open("new_mapping_mouse/AF-"+id+"-F1-model_v2.fasta", "r")
    header = fr.readline()
    data = fr.readlines()
    fr.close()

    data_str = ''
    for item in data:
        data_str += item

#Amino Acid Codes:
#a: {Ala, Gly, Val} => {A, G, V}
#b: {Ile, Leu, Phe, Pro} => {I, L, F, P}
#c: {Tyr, Met, Thr, Ser} => {Y, M, T, S}
#d: {His, Asn, Gln, Trp} => {H, N, Q, W}
#e: {Arg, Lys} => {R, K}
#f: {Asp, Glu} => {D, E}
#g: {Cys} => {C}

    #Convert the amino acids to their corresponding codes
    data1 = data_str.replace("A", "a").replace("G", "a").replace("V", "a")
    data2 = data1.replace("I", "b").replace("L", "b").replace("F", "b").replace("P", "b")
    data3 = data2.replace("Y", "c").replace("M", "c").replace("T", "c").replace("S", "c")
    data4 = data3.replace("H", "d").replace("N", "d").replace("Q", "d").replace("W", "d")
    data5 = data4.replace("R", "e").replace("K", "e")
    data6 = data5.replace("D", "f").replace("E", "f")
    data7 = data6.replace("C", "g")

    #Save the new categorized sequences to the folder: "new_mapping_mouse_categorized"
    fw = open("new_mapping_mouse_categorized/AF-"+id+"-F1-model_v2_categorized.fasta", "w")
    fw.write(header)
    fw.write(data7)
    fw.close()