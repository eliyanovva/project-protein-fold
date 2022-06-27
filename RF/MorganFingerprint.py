from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import DataStructs
import numpy as np
import Globals
import SmileKmer
import numpy as np
import ReadingFasta
import labels
import Globals
import Filtering
import FixedClassificationModel

def Smiles2Finger(smiles):
    m1 = Chem.MolFromSmiles(smiles)
    fp1 = AllChem.GetHashedMorganFingerprint(m1, 2, nBits=1024)
    array = np.zeros((0,), dtype=np.int8)
    DataStructs.ConvertToNumpyArray(fp1, array)
    return array

def ligandmatrix(smiledict, num_proteins):
    fingerprints = []
    for smile in list(smiledict.keys()):
        for i in range(num_proteins):
            fingerprints.append(Smiles2Finger(smiledict[smile])) 
    return np.array(fingerprints)

#Create classification dictionary
logFC, pVal = labels.labels()
classified, pos_counts, neg_counts = labels.classified_logFC_pVal(logFC, pVal)

#Initialize Variables
#categorized variables
categorized_features = set()
categorized_seqs = []
categorized_matrix = []
#3Di variables
di_features = set()
di_seqs = []
di_matrix = []

#Creating output for categorized amino acids
#Read fasta file
fasta1 = open("../AminoAcidSequences/fully_categorized.fasta")
#Create kmer frequency dictionary
seqvar1, features1 = ReadingFasta.make_seqvar(fasta1, categorized_seqs, categorized_features)
#Remove insignificant kmers
filter_feat = Filtering.richness_protein(features1, seqvar1, pos_counts, neg_counts)
# Make the matrix
AA_mat = ReadingFasta.makematrix(seqvar1, filter_feat, categorized_matrix)

#Creating output for 3Di sequences
# Read fasta file
fasta2 = open("../3DiSequences/fullset_ss.fasta")
#Create kmer frequency dictionary
seqvar2, features2 = ReadingFasta.make_seqvar(fasta2, di_seqs, di_features)
#Remove insignificant kmers
filter_feat2 = Filtering.richness_protein(features2, seqvar2, pos_counts, neg_counts)
# Make the matrix
Di_mat = ReadingFasta.makematrix(seqvar2, filter_feat2, di_matrix)

#Concatenate AA and 3Di matrices
intermed_matrix = np.concatenate((np.array(AA_mat, dtype = np.uint8), np.array(Di_mat, dtype = np.uint8)) , axis = 1)
#Expand the protein matrix to account for the ligands
ligand_count = 55
proteins_matrix = np.repeat(intermed_matrix, repeats = ligand_count, axis = 0)

#Import dictionary matching ligands to SMILES String
ligand_dict = Globals.initialize_ligand_dict()
#Create ligands matrix
ligand_matrix = ligandmatrix(ligand_dict, 1084)


#Concatenate protein and ligand matrices
final_matrix = np.concatenate((proteins_matrix, np.array(ligand_matrix, dtype = np.uint8)), axis = 1)

#Create Classification Vector
proteins = seqvar1
logFCmat = []
for protein in proteins:
    for ligand in list(ligand_dict.keys()):
        logFCmat.append(float(classified[str(protein.name)][ligand]))

FixedClassificationModel.train(final_matrix, logFCmat)