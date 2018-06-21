import numpy as np


# Using the input_template.cfg, make num_copies new files enumerated as 0....num_copies,
# while replacing the fields stated in find_replace.
# Use '@' as a wildcard to be replaced by the the copy number.

template_filename = 'input_template.cfg'
input_basename   = 'fes_'
num_copies = 25


find_replace = {
    "idmrg::on"                 : 'true',
    "idmrg::max_steps"          : '500000',
    "idmrg::chi_max"            : '@',
    "precision::eigThreshold"   : '1e-10',
    "hdf5::output_filename"     : "fes_@.h5"
}



# This function replaces entries found in the string "line":
#
#   line.split()[pos] -> val
#
# and wildcards "@" that may be present in "val":
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
        with open(input_basename + str(i) + '.cfg', 'w') as new_input:
            for line in input:
                for var,val in find_replace.items():
                    if var in line.split():
                        if(var == "idmrg::chi_max"):
                            line = replace_value(line, 2, val, 4*i+8)
                        else:
                            line = replace_value(line, 2, val, i)
                new_input.write(line)
