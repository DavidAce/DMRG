import numpy as np
import os


# Using the input_template.cfg, make num_copies new files enumerated as 0....num_copies,
# while replacing the fields stated in find_replace.
# Use '@' as a wildcard to be replaced by the the copy number.

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


def generate_input_files(settings,input_filenames,template_filename, location='-1'):
    assert(len(settings) == len(input_filenames))
    num_files = len(input_filenames)
    if location == '-1':
        location = ''
    else:
        location = location + '/'
    print("Writing ", num_files, " files to ", location)
    for i in range (0,num_files):
        canonical_filename = os.path.normpath(location + input_filenames[i])
        if location != '':
            os.makedirs(os.path.dirname(canonical_filename), exist_ok=True)
        with open(template_filename, 'r') as template:
            with open(canonical_filename, 'w') as file:
                for line in template:
                    for var,val in settings[i].items():
                        if var in line.split():
                            line = replace_value(line, 2, val, i)
                    file.write(line)

