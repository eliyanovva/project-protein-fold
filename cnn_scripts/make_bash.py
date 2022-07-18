import os

command = '#!/bin/bash\n'

for filename in os.listdir('./new_mouse_with_mut/'):
    command += '\ncat monomer_test9.sh | sed \'s/AER304C/'+filename+'/\' > to_run_'+filename[:-6]+' ; sbatch to_run_'+filename[:-6]


bash_script = open("alphafold.sh", "w")
bash_script.write(command)
bash_script.close()