import pandas as pd

ligand_dict = {}
df = pd.read_csv('../Ligands_withSMILE/ligand_SMILEs.csv')

files = df['ligand file'].tolist()
smiles = df['SMILE'].tolist()
for i in range(len(files)):
    ligand_dict[files[i]] = smiles[i]

def initialize_ligand_dict():
    files = df['ligand file'].tolist()
    smiles = df['SMILE'].tolist()
    for i in range(len(files)):
        ligand_dict[files[i]] = smiles[i]

    return ligand_dict

def initialize_ligand_list():
    ligands = ['pS6_DE_1p_citronellol.csv', 'pS6_DE_1p_isoamylAcetate.csv', 'pS6_DE_1p_ethylTiglate.csv',
               'pS6_DE_1p_bIonone.csv', 'pS6_DE_1p_butyricAcid.csv',
               'pS6_DE_1p_paraCresol.csv', 'pS6_DE_1p_bCaryophyllene.csv', 'pS6_DE_p1_isovalericAcid.csv',
               'pS6_DE_1p_Octanal.csv', 'pS6_DE_1p_heptanal.csv',
               'pS6_DE_1p_tbm.csv', 'pS6_DE_1p_bDamascone.csv', 'pS6_DE_1p_pyridine.csv', 'pS6_DE_1p_propionicAcid.csv',
               'pS6_DE_1p_methylSalicylate.csv',
               'pS6_DE_p01_e2butene1thiol.csv', 'pS6_DE_1p_3methyl1butanethiol.csv', 'pS6_DE_1p_ethylButyrate.csv',
               'pS6_DE_1p_hexylTiglate.csv', 'pS6_DE_1p_indole.csv',
               'pS6_DE_500mM_2propylthietane.csv', 'pS6_DE_1p_2heptanone.csv', 'pS6_DE_p01_cyclopentanethiol.csv',
               'pS6_DE_1p_dimethyltrisulfide.csv',
               'pS6_DE_1p_guaiacol.csv', 'pS6_DE_1p_Benzaldehyde.csv', 'pS6_DE_p01_citral.csv',
               'pS6_DE_3mM_androstenone.csv', 'pS6_DE_100p_ebFarnesene.csv',
               'pS6_DE_1p_acetophenone.csv', 'pS6_DE_1p_transCinnamaldehyde.csv', 'pS6_DE_1p_linalool.csv',
               'pS6_DE_1p_2hexanone.csv', 'pS6_DE_1p_isopropylTiglate.csv',
               'pS6_DE_1p_aPinene.csv', 'pS6_DE_1p_diacetyl.csv', 'pS6_DE_1p_geranoil.csv',
               'pS6_DE_1p_heptanoicAcid.csv']
    return ligands

def initialize_protein_list():
    """
    acc_ids = ['Q8VET1', 'Q8VFR3', 'Q8VGZ7', 'Q8VG59', 'Q8VF91', 'Q8VGD7', 'Q8VGC8', 'Q8VFM9', 'A2AVB5', 'L7N1X3', 'Q8VEZ3', 'A0A1L1SQF6', 'Q8VGC7', 'Q8VGC6',
               'A0A1D5RLR5', 'Q7TQV4', 'Q7TQV7', 'Q8VH21', 'Q8VFZ6', 'Q8VFZ2', 'Q9EQ90', 'Q8VGJ3', 'E9Q840', 'Q8VGJ7', 'Q8VGJ5', 'Q9R0K2', 'Q60881', 'Q8VFU6',
               'E9PYP4', 'Q8VFU9', 'Q8VG27', 'E9Q1P0', 'Q7TS48', 'Q8VH10', 'Q8VGM2', 'L7N1Y6', 'Q7TQU8', 'Q8VES9', 'Q8VGS0', 'Q8VGS3', 'Q7TRH8', 'Q8VGS7',
               'E9Q5F1', 'Q8VGK2', 'Q8VFE7', 'E9Q2B9', 'L7MU75', 'Q8VGK8', 'Q8VF14', 'Q7TRM0', 'Q8VF13', 'L7N1Z5', 'Q8VGA1', 'A3KPP7', 'Q8VGX8', 'Q8VFT4',
               'Q8VFT6', 'Q8VGX5', 'Q8VGA9', 'Q8VGX7', 'L7N1Z8', 'Q8VGY9', 'Q8VGX2', 'Q8VGW8', 'Q8VH01', 'Q8VH00', 'Q7TRX7', 'Q7TR67', 'E9Q542', 'Q8VGP9',
               'A0A1L1SUC9', 'Q7TQT8', 'Q8VFE3', 'L7N210', 'Q0VEV7', 'Q8VEZ5', 'Q8VET5', 'Q8VFD2', 'Q7TRX0', 'Q7TRX3', 'Q8VF05', 'Q8VF06', 'Q8VGH1', 'Q8VGH2',
               'E9Q5P8', 'Q8VFW8', 'Q8VF18', 'Q3KPA8', 'Q8VFW5', 'Q8VGY6', 'Q8VGN0', 'Q8VFW3', 'Q8VGQ3', 'Q8VGR9', 'Q8VGQ7', 'Q8VGQ6', 'Q8VF17', 'Q8VGY4',
               'Q7TR71', 'Q8VFE8', 'E9Q0M4', 'Q8VG01', 'Q8VEU1', 'L7N205', 'Q8VEU7', 'Q8VEU6', 'Q7TRN7', 'Q8VFG9', 'Q8VH20', 'Q8VGI3', 'Q8VGI2', 'Q8VFG1', 'Q8VGI4',
               'Q7TQY1', 'Q8VFA9', 'Q8VFX1', 'Q7TR49', 'Q8VEV7', 'A2RT31', 'Q9EQ84', 'Q8VF64', 'Q8VGG9', 'A2ATE5', 'Q8VGV6', 'Q8VGV7', 'Q8VGV1', 'A0A1L1SQT2',
               'Q8VG72', 'Q8VG04', 'Q8VG06', 'Q8VG07', 'Q8VEY3', 'Q924H8', 'Q7TRV2', 'Q8VG09', 'Q7TQR3', 'Q7TQR2', 'Q8VG66', 'Q8VGL9', 'Q8VG64', 'Q8VGL0', 'Q8VGL3',
               'Q8VGL4', 'Q8VGZ6', 'Q8VGJ1', 'Q60895', 'Q8VEW1', 'A2ATG2', 'Q8VGD8', 'Q8VFQ4', 'Q8VFQ3', 'Q8VFQ1', 'Q8VF52', 'Q7TRD9', 'E9Q438', 'Q8VGW6', 'Q8VGW1',
               'Q8VF57', 'Q920Z2', 'Q60882', 'Q8VF58', 'Q8VG90', 'Q8VGW9', 'Q8VG96', 'Q7TRD6', 'Q8VG94', 'Q60893', 'Q8VG16', 'Q8VG14', 'Q7TQQ4', 'Q7TQQ7', 'Q9EQQ5',
               'Q8VG77', 'Q8VGP6', 'Q9QWU6', 'Q7TS37', 'Q8VFC4', 'Q8VFC6', 'Q8VFK1', 'Q7TRE5', 'Q7TRE4', 'Q60885', 'E9Q0Y7', 'Q7TS38', 'Q8VFP4', 'Q60889', 'Q5NC59',
               'Q8VGE3', 'Q8VGE1', 'E9Q546', 'Q8VFX4', 'Q8VGT2', 'E9Q545', 'Q8VGT4', 'Q8VFX2', 'Q7TRA7', 'K7N609', 'Q9JHB2', 'E9Q549', 'Q7TRT9', 'Q9EQA6', 'Q9EPG2',
               'E9PWU0', 'E9Q985', 'Q2M2Q2', 'Q8VG42', 'Q8VG49', 'Q8VF43', 'Q7TRE6', 'Q8VFG4', 'Q7TRB0', 'Q8VGB9', 'Q7TR59', 'E9Q413', 'Q8VGB3', 'Q8VFJ5', 'Q7TRB9',
               'Q8VGB4', 'Q7TRU6', 'Q8VF78', 'Q8VGU9', 'Q8VGU8', 'Q8VGU7', 'E9Q3K2', 'Q8VGU3', 'Q8VF72', 'Q7TS51', 'E9Q0Q2', 'Q8VGM3', 'Q7TRJ1', 'Q0VEL5']
    """
    acc_ids = []
    fr = open("../AminoAcidSequences/allsequences.fasta", "r")
    lines = fr.readlines()
    #fw = open("../AminoAcidSequences/new_accessions.txt", "w")
    for line in lines:
        if line[0] == ">":
            acc_ids.append(line[1:-1])
            #fw.write(str(line[1:-1]) + "\n")
    fr.close()
    #fw.close()

    return acc_ids

