import numpy as np 

def initialize():
    #categorized variables
    global categorized_features
    categorized_features = set()
    global categorized_seqs
    categorized_seqs = []
    global categorized_matrix
    categorized_matrix = []
    #3Di variables
    global di_features
    di_features = set()
    global di_seqs
    di_seqs = []
    global di_matrix
    di_matrix = []
    #ligand variables
    global ligand_dict
    ligand_dict = {'2hepatanone': 'CCCCCC(=O)C','2hexanone': 'CCCCC(=O)C','3methyl1butanethiol': 'CC(C)CCS',
'acetophenone': 'CC(=O)C1=CC=CC=C1','aPinene': 'CC1=CCC2CC1C2(C)C','bCaryophyllene': 'CC1=CCCC(=C)C2CC(C2CC1)(C)C',
'bDamascone': 'CC=CC(=O)C1=C(CCCC1(C)C)C','Benzaldehyde': 'C1=CC=C(C=C1)C=O','bIonone': 'CC1=C(C(CCC1)(C)C)C=CC(=O)C',
'butyricAcid': 'CCCC(=O)O','citronellol': 'CC(CCC=C(C)C)CCO','diacetyl': 'CC(=O)C(=O)C','dimethylSulfide': 'CSC',
'dimethyltrisulfide': 'CSSSC','ethylButyrate': 'CCCC(=O)OCC','ethylTiglate': 'CCOC(=O)C(=CC)C','geranoil': 'CC(=CCCC(=CCO)C)C',
'guaiacol': 'COC1=CC=CC=C1O','heptanal': 'CCCCCCC=O','heptanoicAcid': 'CCCCCCC(=O)O','hexylTiglate': 'CCCCCCOC(=O)C(=CC)C',
'indole': 'C1=CC=C2C(=C1)C=CN2','isoamyl acetate': 'CC(C)CCOC(=O)C','isopropyl tiglate': 'CC=C(C)C(=O)OC(C)C',
'linalool': 'CC(=CCCC(C)(C=C)O)C','methylSalicylate': 'COC(=O)C1=CC=CC=C1O','Octanal': 'CCCCCCCC=O','paraCresol': 'CC1=CC=C(C=C1)O',
'propionicAcid': 'CCC(=O)O','pyridine': 'C1=CC=NC=C1','tbm': 'CC(C)(C)S','transCinnamaldehyde': 'C1=CC=C(C=C1)C=CC=O',
'androstenone': 'CC12CCC3C(C1CC=C2)CCC4C3(CCC(=O)C4)C','ebFarnesene': 'CC(=CCCC(=CCCC(=C)C=C)C)C','2proplythietane': 'CCCC1CCS1',
'citral': 'CC(=CCCC(=CC=O)C)C','cyclopentanethiol': 'C1CCC(C1)S','e2butene1thiol': 'CC=CCS','isovalericAcid': 'CC(C)CC(=O)O'}

