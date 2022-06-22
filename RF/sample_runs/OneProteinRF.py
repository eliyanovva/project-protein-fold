#Protein is A0A1L1SQF6
import PreparingMatrix
import CombineLigandsProteins
import ReadingFasta
import RandomForest
import SmileKmer

ligand_dict = {"pS6_DE_1p_citronellol.csv":'CC(CCC=C(C)C)CCO', "pS6_DE_1p_isoamylAcetate.csv":'CC(C)CCOC(=O)C', "pS6_DE_1p_ethylTiglate.csv":'CCOC(=O)C(=CC)C',
"pS6_DE_1p_bIonone.csv":'CC1=C(C(CCC1)(C)C)C=CC(=O)C', "pS6_DE_1p_butyricAcid.csv":'CCCC(=O)O', "pS6_DE_1p_paraCresol.csv":'CC1=CC=C(C=C1)O', "pS6_DE_1p_bCaryophyllene.csv":'CC1=CCCC(=C)C2CC(C2CC1)(C)C',
"pS6_DE_p1_isovalericAcid.csv":'CC(C)CC(=O)O', "pS6_DE_1p_Octanal.csv":'CCCCCCCC=O', "pS6_DE_1p_heptanal.csv":'CCCCCCC=O', "pS6_DE_1p_tbm.csv":'CC(C)(C)S', "pS6_DE_1p_bDamascone.csv":'CC=CC(=O)C1=C(CCCC1(C)C)C',
"pS6_DE_1p_pyridine.csv":'C1=CC=NC=C1', "pS6_DE_1p_propionicAcid.csv":'CCC(=O)O', "pS6_DE_1p_methylSalicylate.csv":'COC(=O)C1=CC=CC=C1O', "pS6_DE_p01_e2butene1thiol.csv":'CC=CCS',
"pS6_DE_1p_3methyl1butanethiol.csv":'CC(C)CCS', "pS6_DE_1p_ethylButyrate.csv":'CCCC(=O)OCC', "pS6_DE_1p_hexylTiglate.csv":'CCCCCCOC(=O)C(=CC)C', "pS6_DE_1p_indole.csv":'C1=CC=C2C(=C1)C=CN2',
"pS6_DE_500mM_2propylthietane.csv":'CCCC1CCS1', "pS6_DE_1p_2heptanone.csv":'CCCCCC(=O)C', "pS6_DE_p01_cyclopentanethiol.csv":'C1CCC(C1)S',
"pS6_DE_1p_dimethyltrisulfide.csv":'CSSSC', "pS6_DE_1p_guaiacol.csv":'COC1=CC=CC=C1O', "pS6_DE_1p_Benzaldehyde.csv":'C1=CC=C(C=C1)C=O', "pS6_DE_p01_citral.csv":'CC(=CCCC(=CC=O)C)C',
"pS6_DE_3mM_androstenone.csv":'CC12CCC3C(C1CC=C2)CCC4C3(CCC(=O)C4)C', "pS6_DE_100p_ebFarnesene.csv":'CC(=CCCC(=CCCC(=C)C=C)C)C', "pS6_DE_1p_acetophenone.csv":'CC(=O)C1=CC=CC=C1',
"pS6_DE_1p_transCinnamaldehyde.csv":'C1=CC=C(C=C1)C=CC=O', "pS6_DE_1p_linalool.csv":'CC(=CCCC(C)(C=C)O)C', "pS6_DE_1p_2hexanone.csv":'CCCCC(=O)C', "pS6_DE_1p_isopropylTiglate.csv":'CC=C(C)C(=O)OC(C)C',
"pS6_DE_1p_aPinene.csv":'CC1=CCC2CC1C2(C)C', "pS6_DE_1p_diacetyl.csv":'CC(=O)C(=O)C', "pS6_DE_1p_geranoil.csv":'CC(=CCCC(=CCO)C)C', "pS6_DE_1p_heptanoicAcid.csv":'CCCCCCC(=O)O'}

SmileKmer.importmatrix(ligand_dict, 5, 1)
ligand_matrix = SmileKmer.ligmat

CombineLigandsProteins.exportdicts()
logFC = CombineLigandsProteins.logdic
p = CombineLigandsProteins.pdic

logmat = []
pmat = []
for ligand in list(logFC['A0A1L1SQF6'].keys()):
    logmat.append(logFC['A0A1L1SQF6'][ligand])
    pmat.append(p['A0A1L1SQF6'][ligand])

RandomForest.train(ligand_matrix, logmat)
RandomForest.train(ligand_matrix, pmat)
