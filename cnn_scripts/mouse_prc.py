from moleculekit.molecule import Molecule
from moleculekit.tools.voxeldescriptors import getVoxelDescriptors, viewVoxelFeatures
from moleculekit.tools.atomtyper import prepareProteinForAtomtyping
from moleculekit.smallmol.smallmol import SmallMol
from moleculekit.home import home

def dict_to_fasta(d, name):
    """
    Writes a fasta file from a dictionary where keys are sequence names and
    values are sequences
    :param d: dictionary
    :param dir: directory for fasta file to be stored
    :param name: filename to write
    :return:
    """
    with open(f'{name}.fa', 'w') as f:
        for k, v in d.items():
            f.write(f'>{k}\n')
            f.write(f'{v}\n')

mut_dict = dict()

Olfr124 = Molecule('/home/users/bmp40/project-protein-fold/cnn_scripts/Olfr124.pdb')
Olfr124_seq = list(Olfr124.sequence().values())[0]

Olfr1362 = Molecule('/home/users/bmp40/project-protein-fold/cnn_scripts/Olfr1362.pdb')
Olfr1362_seq = list(Olfr1362.sequence().values())[0]
print(Olfr1362_seq[101])

mut_dict.update({'Olfr124_A183I': Olfr124_seq[0:182] + 'I' + Olfr124_seq[183:]})
mut_dict.update({'Olfr124_H176A': Olfr124_seq[0:175] + 'A' + Olfr124_seq[176:]})
mut_dict.update({'Olfr124_P182A': Olfr124_seq[0:181] + 'A' + Olfr124_seq[182:]})
mut_dict.update({'Olfr1362_F102Y': Olfr1362_seq[0:101] + 'Y' + Olfr1362_seq[102:]})
mut_dict.update({'Olfr1362_I107L': Olfr1362_seq[0:106] + 'L' + Olfr1362_seq[107:]})
mut_dict.update({'Olfr1362_I204A': Olfr1362_seq[0:203] + 'A' + Olfr1362_seq[204:]})
mut_dict.update({'Olfr1362_S106A': Olfr1362_seq[0:105] + 'A' + Olfr1362_seq[106:]})
mut_dict.update({'Olfr1362_V207L': Olfr1362_seq[0:206] + 'L' + Olfr1362_seq[207:]})

dict_to_fasta(mut_dict, 'mut_mouse_fasta')