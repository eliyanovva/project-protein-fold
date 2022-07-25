import os
from os.path import exists

human_path1 = '/home/users/bmp40/share/deorphanOR_human/'
human_path2 = '/home/users/bmp40/share/orphanOR_human/'

for filename in os.listdir('./new_fasta'):
    #print(filename[:-6].upper())
    if exists(os.path.join('./new_mouse_with_mut', filename)):
        print('exists')
        os.remove(os.path.join('./new_fasta', filename))
    elif exists(os.path.join(human_path1,'AF-' + filename[:-6].upper() + '-F1-model_v2.pdb')):
        print('exists')
        os.remove(os.path.join('./new_fasta', filename))
    elif exists(os.path.join(human_path2,'AF-' + filename[:-6].upper() + '-F1-model_v2.pdb')):
        print('exists')
        os.remove(os.path.join('./new_fasta', filename))
    else:
        print('not exists')