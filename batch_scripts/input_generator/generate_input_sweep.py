from generate_inputs import *
import numpy as np


# Using the input_template.cfg, make num_copies new files enumerated as 0....num_copies,
# while replacing the fields stated in find_replace.
# Use '@' as a wildcard to be replaced by the the copy number.

template_filename = 'input_template.cfg'
basename    = 'mbl_'
location    = "../input/"
settings = []
input_filenames = []
realizations = np.arange(0,5,1)  # Number of copies for each point on the sweep
lengths = [24] #np.arange(12,24,2)
lambdas = np.arange(0, 0.2, 0.1)
J_log_mean = [1]
h_log_mean = np.flipud(np.arange(-3, 3.5, 0.5)+1)
num_total = 0

for num_L in range(len(lengths)):
    for num_l in range(len(lambdas)):
        for num_j in range(len(J_log_mean)):
            for num_h in range(len(h_log_mean)):
                for num_r in range(len(realizations)):
                    # input_filenames.append(basename + 'l' + str(num_l) + '_J'+ str(num_j) + '_h'+ str(num_h) + '_' + str(num_r) + '.cfg')
                    input_filenames.append('L_'+ str(lengths[num_L]) + '/' +  basename + str(num_total) + '_l' + str(num_l) + '_J'+ str(num_j) + '_h'+ str(num_h) + '_' + str(num_r) + '.cfg')
                    settings.append({
                        "model::selfdual_tf_rf_ising::J_log_mean"     : "{:.2f}".format(J_log_mean[num_j]),
                        "model::selfdual_tf_rf_ising::h_log_mean"     : "{:.2f}".format(h_log_mean[num_h]),
                        "model::selfdual_tf_rf_ising::lambda"         : "{:.2f}".format(lambdas[num_l]),
                        "model::selfdual_tf_rf_ising::J_sigma"        : "0.5",
                        "model::selfdual_tf_rf_ising::h_sigma"        : "0.5",
			            "xdmrg::max_length"                           : "24",
			            "xdmrg::chi_max"                              : "32",
                        "xdmrg::seed"                                 : str(num_total),
                        "hdf5::output_folder"                         : 'output/L_'+ str(lengths[num_L]) + '/l_'+str(num_l) + '/J_' +str(num_j) + '/h_'+ str(num_h),
                        "hdf5::output_filename"                       : basename + str(num_r) + '.h5'
                    })
                    num_total  = num_total + 1

generate_input_files(settings, input_filenames, template_filename,location)



