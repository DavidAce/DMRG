import numpy as np


# Using the input_template.cfg, make num_copies new files enumerated as 0....num_copies,
# while replacing the fields stated in find_replace.
# Use '@' as a wildcard to be replaced by the the copy number.

template_filename = 'input_template.cfg'
input_basename   = 'mbl_sg_'
num_copies = 500


find_replace = {
    "model::selfdual_tf_rf_ising::J_mu"     : '2     ',
    "model::selfdual_tf_rf_ising::h_mu"     : '0.01  ',
    "model::selfdual_tf_rf_ising::J_sigma"  : '1     ',
    "model::selfdual_tf_rf_ising::h_sigma"  : '1     ',
    "model::selfdual_tf_rf_ising::lambda"   : '0.2   ',
    "model::selfdual_tf_rf_ising::d"        : '2     ',
    "xdmrg::on"                             : 'true  ',
    "xdmrg::max_length"                     : '24    ',
    "xdmrg::max_sweeps"                     : '50    ',
    "xdmrg::chi_max"                        : '16    ',
    "xdmrg::chi_grow"                       : 'false ',
    "xdmrg::seed"                           : '@     ',
    "xdmrg::print_freq"                     : '1     ',
    "xdmrg::store_freq"                     : '1     ',
    "hdf5::output_filename"                 : input_basename + '@.h5'
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
                        line = replace_value(line, 2, val, i)
                new_input.write(line)
