import numpy as np
import os



def replace_value(line,pos,val):
    old_val = line.split()[pos]
    index_start = line.find(old_val)
    index_end   = index_start + len(old_val)
    len_diff = len(old_val) - len(val)
    return line[:index_start] + val + ' '*len_diff + line[index_end:]


def generate_input_file(settings,input_filename,template_filename):
    with open(template_filename, 'r') as template:
        with open(input_filename, 'w') as file:
            for line in template:
                for var,val in settings.items():
                    if line.find(var) >= 0:
                        line = replace_value(line, 2, val)
                file.write(line)

