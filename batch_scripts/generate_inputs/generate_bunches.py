import numpy as np
import os
import random

parentpath = os.path.abspath('..')
os.chdir(parentpath)

src_directory='input'
tgt_directory='bunch'
bunch_size = 1000
bunch_list = []
if not os.path.exists(tgt_directory):
    os.makedirs(tgt_directory)

def chunks(l, n):
    # For item i in a range that is a length of l,
    for i in range(0, len(l), n):
        # Create an index range for l of n items:
        yield l[i:i+n]

current_bunch = []
for dirName, subdirList, fileList in os.walk(src_directory):
    for src_filename in fileList:
        current_bunch.append(dirName + '/' +src_filename)

#current_bunch.sort()
random.shuffle(current_bunch)
bunch_list = chunks(current_bunch,bunch_size)


for i, bunch in enumerate(bunch_list):
    filepath = tgt_directory + '/bunch_' + str(i) + '.txt'
    with open(filepath, 'w') as file_handler:
        for item in bunch:
            file_handler.write("{}\n".format(item))




