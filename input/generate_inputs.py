import numpy as np


# Using the input_template.cfg, make num_copies new files enumerated as 0....num_copies,
# while replacing the fields stated in find_replace.
# Use '@' as a wildcard to be replaced by the the copy number.

template_filename = 'input_template.cfg'
output_basename   = 'test_'
num_copies = 5


find_replace = {
    "idmrg::max_steps"      : '3000',
    "hdf5::output_filename": "testoutput_@.h5",
    "model::seed": '@'
    }



# This function replaces entries found in the string "line":
#
#   line.split()[pos] -> val
#
# and wildcards "@" that may be presen in "val":
#
#   @ -> num
#
# while keeping the same effective width of the entry.
def replace_value(line,pos,val,num=0):
    val = val.replace('@', str(num))
    old_val = line.split()[pos]
    len_diff = len(val) - len(old_val)
    return line.replace(old_val + ' ' * len_diff, val.ljust(len(old_val)))


for i in range(num_copies):
    with open(template_filename, 'r') as input:
        with open(output_basename + str(i) + '.cfg', 'w') as new_input:
            for line in input:
                for var,val in find_replace.items():
                    if var in line.split():
                        line = replace_value(line, 2, val, i)
                new_input.write(line)
