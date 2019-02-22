import numpy as np
import os


parentpath = os.path.abspath('..')
os.chdir(parentpath)

src_directory='input'
tgt_directory='bunch'
bunch_size = 100
bunch_list = []

for dirName, subdirList, fileList in os.walk(src_directory):
    current_bunch = []
    for src_filename in fileList:
        current_bunch.append(dirName + '/' +src_filename)
        if len(current_bunch) >= bunch_size:
            bunch_list.append(current_bunch)
            current_bunch = []

if not os.path.exists(tgt_directory):
    os.makedirs(tgt_directory)

for i, bunch in enumerate(bunch_list):
    filepath = tgt_directory + '/bunch_' + str(i) + '.txt'
    with open(filepath, 'w') as file_handler:
        for item in bunch:
            file_handler.write("{}\n".format(item))
