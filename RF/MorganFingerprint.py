from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import DataStructs
import numpy as np
import Globals

smiledict = Globals.initialize_ligand_dict()

def Smiles2Finger(smiles):
    m1 = Chem.MolFromSmiles(smiles)
    fp1 = AllChem.GetHashedMorganFingerprint(m1, 2, nBits=1024)
    array = np.zeros((0,), dtype=np.int8)
    DataStructs.ConvertToNumpyArray(fp1, array)
    return array

fingerprints = []

for smile in list(smiledict.keys()):
    fingerprints.append(Smiles2Finger(smiledict[smile])) 

print(np.array(fingerprints))
