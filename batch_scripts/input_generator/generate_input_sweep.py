from generate_inputs import *
import numpy as np


# Using the input_template.cfg, make num_copies new files enumerated as 0....num_copies,
# while replacing the fields stated in find_replace.
# Use '@' as a wildcard to be replaced by the the copy number.

template_filename = 'input_template.cfg'
basename    = 'mbl_'
location    = "../input/"


copies     = 500    #Number of copies for each point on the sweep
J_log_mean = {1}
h_log_mean = range(-5,5,1)

settings = []
input_filenames = []
num_total  = 0
num_batch  = 0
for J in J_log_mean:
    for h in h_log_mean:
        for copy in range(0,copies):
            input_filenames.append(basename + str(num_batch) + '_' + str(copy) + '.cfg')
            settings.append( {
                "model::selfdual_tf_rf_ising::J_log_mean"     : str(J),
                "model::selfdual_tf_rf_ising::h_log_mean"     : str(h),
                "xdmrg::seed"                                 : str(num_total),
                "hdf5::output_folder"                         : 'output_' + basename + str(num_batch),
                "hdf5::output_filename"                       : basename + str(copy) + '.h5'
            })
            num_total  = num_total + 1
        num_batch = num_batch + 1
print(settings)

generate_input_files(settings, input_filenames, template_filename,location)



