import numpy as np
import Duplicates
import CombineLigandsProteins
import Globals
import labels
import ReadingFasta
import SmileKmer

#Import dictionary matching ligands to SMILES String
ligand_dict = Globals.initialize_ligand_dict()

CombineLigandsProteins.import_final()
logFC = CombineLigandsProteins.logFC_data
FDR = CombineLigandsProteins.FDR_data
AA_filter_TM3 = CombineLigandsProteins.AA3_kmers
AA_filter_TM5 = CombineLigandsProteins.AA5_kmers
AA_filter_TM6 = CombineLigandsProteins.AA6_kmers
AA_filter_TM7 = CombineLigandsProteins.AA7_kmers
Di_filter_TM3 = CombineLigandsProteins.Di3_kmers
Di_filter_TM5 = CombineLigandsProteins.Di5_kmers
Di_filter_TM6 = CombineLigandsProteins.Di6_kmers
Di_filter_TM7 = CombineLigandsProteins.Di7_kmers
filter_kmers = CombineLigandsProteins.filter_kmers

neutral_ligands = ['pS6_DE_1p_2hac.csv', 'pS6_DE_1p_3methyl1butanethiol.csv', 'pS6_DE_1p_Benzaldehyde.csv', 'pS6_DE_1p_aPinene.csv',
             'pS6_DE_1p_bDamascone.csv', 'pS6_DE_1p_dimethyltrisulfide.csv', 'pS6_DE_1p_heptanal.csv', 'pS6_DE_1p_indole.csv',
             'pS6_DE_1p_isoamylAcetate.csv', 'pS6_DE_1p_isopropylTiglate.csv', 'pS6_DE_1p_methylSalicylate.csv',
             'pS6_DE_1p_nMenthol.csv', 'pS6_DE_1p_pCarvone.csv', 'pS6_DE_1p_sbt.csv', 'pS6_DE_1p_tmt.csv',
             'pS6_DE_3mM_androstenone.csv', 'pS6_DE_500mM_2propylthietane.csv',	'pS6_DE_p01_e2butene1thiol.csv']
neutral_ligands.sort()

possible = labels.extract_highvals(logFC, FDR, neutral_ligands)
n_proteins = []

f = open('possible_pos_pairs.txt', 'w')

for id in possible:
    if len(possible[id]) >= 5:
        n_proteins.append(id)
        f.write(id)
        for lig in possible[id]:
            f.write("," + lig)
        f.write("\n")
f.close()
n_proteins.sort()

nAA_dict = Globals.initialize_AA_dict(n_proteins)
nDi_dict = Globals.initialize_3Di_dict(n_proteins)

nAA_seqvar_TM3, ignoreAA3 = ReadingFasta.make_seqvar_TMS(nAA_dict, 0, 5, {}, set())
nAA_seqvar_TM5, ignoreAA5 = ReadingFasta.make_seqvar_TMS(nAA_dict, 1, 5, {}, set())
nAA_seqvar_TM6, ignoreAA6 = ReadingFasta.make_seqvar_TMS(nAA_dict, 2, 5, {}, set())
nAA_seqvar_TM7, ignoreAA7 = ReadingFasta.make_seqvar_TMS(nAA_dict, 3, 5, {}, set())

nDi_seqvar_TM3, ignoreDi3 = ReadingFasta.make_seqvar_TMS(nDi_dict, 0, 5, {}, set())
nDi_seqvar_TM5, ignoreDi5 = ReadingFasta.make_seqvar_TMS(nDi_dict, 1, 5, {}, set())
nDi_seqvar_TM6, ignoreDi6 = ReadingFasta.make_seqvar_TMS(nDi_dict, 2, 5, {}, set())
nDi_seqvar_TM7, ignoreDi7 = ReadingFasta.make_seqvar_TMS(nDi_dict, 3, 5, {}, set())

nAA_seqvar = [nAA_seqvar_TM3, nAA_seqvar_TM5, nAA_seqvar_TM6, nAA_seqvar_TM7]
nAA_feat = [AA_filter_TM3, AA_filter_TM5, AA_filter_TM6, AA_filter_TM7]
nDi_seqvar = [nDi_seqvar_TM3, nDi_seqvar_TM5, nDi_seqvar_TM6, nDi_seqvar_TM7]
nDi_feat = [Di_filter_TM3, Di_filter_TM5, Di_filter_TM6, Di_filter_TM7]
n_unip = Duplicates.n_remove_proteins(nAA_seqvar, nAA_feat, nDi_seqvar, nDi_feat, n_proteins)

n_lig_counts = SmileKmer.n_ligand_matrix(ligand_dict, 5, neutral_ligands, filter_kmers)
n_uni_lig = Duplicates.n_remove_ligands(n_lig_counts)

num_ligands = len(n_uni_lig)

nlig_mat = []
for lig in n_uni_lig:
    nlig_mat.append(np.array(list(n_lig_counts[lig].values())))

nAA_mat_TM3 = ReadingFasta.make_nmatrix(nAA_seqvar_TM3, AA_filter_TM3, [], n_unip, num_ligands)
nAA_mat_TM5 = ReadingFasta.make_nmatrix(nAA_seqvar_TM5, AA_filter_TM5, [], n_unip, num_ligands)
nAA_mat_TM6 = ReadingFasta.make_nmatrix(nAA_seqvar_TM6, AA_filter_TM6, [], n_unip, num_ligands)
nAA_mat_TM7 = ReadingFasta.make_nmatrix(nAA_seqvar_TM7, AA_filter_TM7, [], n_unip, num_ligands)

nDi_mat_TM3 = ReadingFasta.make_nmatrix(nDi_seqvar_TM3, Di_filter_TM3, [], n_unip, num_ligands)
nDi_mat_TM5 = ReadingFasta.make_nmatrix(nDi_seqvar_TM5, Di_filter_TM5, [], n_unip, num_ligands)
nDi_mat_TM6 = ReadingFasta.make_nmatrix(nDi_seqvar_TM6, Di_filter_TM6, [], n_unip, num_ligands)
nDi_mat_TM7 = ReadingFasta.make_nmatrix(nDi_seqvar_TM7, Di_filter_TM7, [], n_unip, num_ligands)

nAA_mat = np.concatenate((np.array(nAA_mat_TM3, dtype= np.uint8), np.array(nAA_mat_TM5, dtype= np.uint8),
                          np.array(nAA_mat_TM6, dtype= np.uint8), np.array(nAA_mat_TM7, dtype= np.uint8)), axis = 1)

nDi_mat = np.concatenate((np.array(nDi_mat_TM3, dtype= np.uint8), np.array(nDi_mat_TM5, dtype= np.uint8),
                          np.array(nDi_mat_TM6, dtype= np.uint8), np.array(nDi_mat_TM7, dtype= np.uint8)), axis = 1)

n_intermed = np.concatenate((np.array(nAA_mat, dtype = np.uint8), np.array(nDi_mat, dtype = np.uint8)) , axis = 1)
n_final_lig = np.repeat(nlig_mat, len(n_unip), axis = 0)

n_final_mat = np.concatenate((n_intermed, n_final_lig), axis=1)

"""
print(len(n_unip))                  #31
print(len(neutral_ligands))         #18
print(len(filter_kmers))            #81
print(len(n_intermed))              #558
print(len(n_intermed[0]))           #4856
print(len(n_uni_lig))               #18
print(len(n_final_lig))             #558
print(len(n_final_lig[0]))          #81
"""


def import_final():
    global nMat
    nMat = n_final_mat
    global nproteins
    nproteins = n_unip
    global nligands
    nligands = n_uni_lig