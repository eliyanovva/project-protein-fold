import numpy as np
import Duplicates
import CombineLigandsProteins
import Globals
import labels
import ReadingFasta
import SmileKmer
import random

#Import dictionary matching ligands to SMILES String
ligand_dict = Globals.initialize_ligand_dict()

CombineLigandsProteins.import_final()
logFC = CombineLigandsProteins.logFC_data
FDR = CombineLigandsProteins.FDR_data
AA_seqvar_TM3 = CombineLigandsProteins.AA3_seqs
AA_seqvar_TM5 = CombineLigandsProteins.AA5_seqs
AA_seqvar_TM6 = CombineLigandsProteins.AA6_seqs
AA_seqvar_TM7 = CombineLigandsProteins.AA7_seqs
Di_seqvar_TM3 = CombineLigandsProteins.Di3_seqs
Di_seqvar_TM5 = CombineLigandsProteins.Di5_seqs
Di_seqvar_TM6 = CombineLigandsProteins.Di6_seqs
Di_seqvar_TM7 = CombineLigandsProteins.Di7_seqs
AA_filter_TM3 = CombineLigandsProteins.AA3_kmers
AA_filter_TM5 = CombineLigandsProteins.AA5_kmers
AA_filter_TM6 = CombineLigandsProteins.AA6_kmers
AA_filter_TM7 = CombineLigandsProteins.AA7_kmers
Di_filter_TM3 = CombineLigandsProteins.Di3_kmers
Di_filter_TM5 = CombineLigandsProteins.Di5_kmers
Di_filter_TM6 = CombineLigandsProteins.Di6_kmers
Di_filter_TM7 = CombineLigandsProteins.Di7_kmers
filter_kmers = CombineLigandsProteins.filter_kmers
unique_proteins = list(CombineLigandsProteins.uni_prot)
unique_ligands = CombineLigandsProteins.uni_lig
lig_counts_filter = CombineLigandsProteins.lig_counts_filter

print(len(unique_proteins))         #313
print(len(unique_ligands))          #26

#https://www.geeksforgeeks.org/randomly-select-n-elements-from-list-in-python/
unip = []
for i in range(20):
    id = random.choice(unique_proteins)
    unip.append(id)
    unique_proteins.remove(id)

new_combos = labels.extract_new_combos(FDR, unip, unique_ligands)

AA_mat_TM3 = ReadingFasta.make_combomatrix(AA_seqvar_TM3, AA_filter_TM3, [], new_combos)
AA_mat_TM5 = ReadingFasta.make_combomatrix(AA_seqvar_TM5, AA_filter_TM5, [], new_combos)
AA_mat_TM6 = ReadingFasta.make_combomatrix(AA_seqvar_TM6, AA_filter_TM6, [], new_combos)
AA_mat_TM7 = ReadingFasta.make_combomatrix(AA_seqvar_TM7, AA_filter_TM7, [], new_combos)

Di_mat_TM3 = ReadingFasta.make_combomatrix(Di_seqvar_TM3, Di_filter_TM3, [], new_combos)
Di_mat_TM5 = ReadingFasta.make_combomatrix(Di_seqvar_TM5, Di_filter_TM5, [], new_combos)
Di_mat_TM6 = ReadingFasta.make_combomatrix(Di_seqvar_TM6, Di_filter_TM6, [], new_combos)
Di_mat_TM7 = ReadingFasta.make_combomatrix(Di_seqvar_TM7, Di_filter_TM7, [], new_combos)

AA_matrix = np.concatenate((np.array(AA_mat_TM3, dtype = np.uint8), np.array(AA_mat_TM5, dtype = np.uint8),
                            np.array(AA_mat_TM6, dtype = np.uint8), np.array(AA_mat_TM7, dtype = np.uint8)) , axis = 1)

Di_matrix = np.concatenate((np.array(Di_mat_TM3, dtype = np.uint8), np.array(Di_mat_TM5, dtype = np.uint8),
                            np.array(Di_mat_TM6, dtype = np.uint8), np.array(Di_mat_TM7, dtype = np.uint8)) , axis = 1)

lig_mat = []            #matrix to store ligand features
for id in unip:
    for lig in new_combos[id]:
            lig_mat.append(np.array(list(lig_counts_filter[lig].values())))

print(len(AA_matrix))           #499
print(len(Di_matrix))           #499
print(len(lig_mat))             #499

intermed_matrix = np.concatenate((np.array(AA_matrix, dtype = np.uint8), np.array(Di_matrix, dtype = np.uint8)) , axis = 1)
final_matrix = np.concatenate((intermed_matrix, np.array(lig_mat, dtype = np.uint8)), axis = 1)

def import_final():
    global X
    X = final_matrix
    global combos
    combos = new_combos