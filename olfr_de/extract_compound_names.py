with open('file_names.txt', 'r') as file_names:
    names = file_names.readlines()
    with open('compound_names.txt', 'w') as write_file:
        for name in names:
            left_index = name.rfind('_')
            right_index = name.find('.csv')
            write_file.write(name[left_index + 1 : right_index] + '\n')