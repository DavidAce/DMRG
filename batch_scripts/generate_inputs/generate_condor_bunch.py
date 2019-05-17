import numpy as np
import os
import random

parentpath = os.path.abspath('..')
os.chdir(parentpath)

src_directory='input'
tgt_directory='bunch'

realizations = 1000
bunch_list = []
seed_counter = 0

if not os.path.exists(tgt_directory):
    os.makedirs(tgt_directory)

def chunks(l, n):
    # For item i in a range that is a length of l,
    for i in range(0, len(l), n):
        # Create an index range for l of n items:
        yield l[i:i+n]

complete_list = []
for dirName, subdirList, fileList in os.walk(src_directory):
    for src_filename in fileList:
        for realization in range(realizations):
            complete_list.append(dirName + '/' + src_filename + ' ' + str(seed_counter))
            seed_counter = seed_counter + 1

random.shuffle(complete_list)
filepath = tgt_directory + '/bunch.condor'
with open(filepath, 'w') as file_handler:
    for i, item in enumerate(complete_list):
        file_handler.write("{}\n".format(item))

print("Seeds        :   ", seed_counter)
print("Sim types    :   ", seed_counter/realizations)
print("Realizations :   ", realizations)