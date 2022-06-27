from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import DataStructs
import numpy as np

def Smiles2Finger(smiles):
    m1 = Chem.MolFromSmiles(smiles)
    fp1 = AllChem.GetHashedMorganFingerprint(m1, 2, nBits=1024)
    array = np.zeros((0,), dtype=np.int8)
    DataStructs.ConvertToNumpyArray(fp1, array)
    return array.nonzero()