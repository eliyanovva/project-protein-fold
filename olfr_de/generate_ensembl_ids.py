with open('./pS6_DE_1p_2e3mp.csv', 'r') as file:
    with open('./ensemble_names.csv', 'w') as write_file:
        lines = file.readlines()
        for line in lines:
            line_list = line.split(',')
            write_file.write(line_list[1][1:-1] + ',')