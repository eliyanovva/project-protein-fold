#This script writes all of the names of the ligands to the compound_names.txt file. 
#These names are the Matsunami Lab's shortened versions of the chemical name.

with open('file_names.txt', 'r') as file_names:
    names = file_names.readlines()
    with open('compound_names.txt', 'w') as write_file:
        for name in names:
            left_index = name.rfind('_')
            right_index = name.find('.csv')
            write_file.write(name[left_index + 1 : right_index] + '\n')