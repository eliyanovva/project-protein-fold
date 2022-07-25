import os

for filename in os.listdir('./'):
    f = os.path.join('./', filename)
    # checking if it is a file
    if os.path.isfile(f):
        if 'to_run' in filename or 'monomer' in filename:
            print(filename)
            os.remove(f)