from src.generate_inputs import *
import numpy as np


# Using the input_template.cfg, make num_copies new files enumerated as 0....num_copies,
# while replacing the fields stated in find_replace.
# Use '@' as a wildcard to be replaced by the the copy number.

template_filename = 'input_template.cfg'
basename    = 'mbl_'
location    = "../input"

realizations = np.arange(0,30,1)  # Number of copies for each point on the sweep
sites        = np.arange(12,36,4)
lambdas      = np.arange(0, 0.1, 0.1)
deltas       = np.array([0]) #np.arange(-1.0, 2.0, 1.0)
J_log_mean   = np.array([1])
h_log_mean   = J_log_mean - deltas
# h_log_mean   = np.flipud(np.arange(-1, 1.5, 0.5)+1)
num_total = 0
settings = []
input_filenames = []

for num_L in sites:
    for num_l in range(len(lambdas)):
        for num_j in range(len(J_log_mean)):
            for num_h in range(len(h_log_mean)):
                for num_r in realizations:
                    input_filenames.append('L_'+ str(num_L) + '/' +  basename + str(num_total) + '_l' + str(num_l) + '_J'+ str(num_j) + '_h'+ str(num_h) + '_' + str(num_r) + '.cfg')
                    settings.append({
                        "model::selfdual_tf_rf_ising::J_log_mean"     : "{:.2f}".format(J_log_mean[num_j]),
                        "model::selfdual_tf_rf_ising::h_log_mean"     : "{:.2f}".format(h_log_mean[num_h]),
                        "model::selfdual_tf_rf_ising::lambda"         : "{:.2f}".format(lambdas[num_l]),
                        "model::selfdual_tf_rf_ising::J_sigma"        : "1.0",
                        "model::selfdual_tf_rf_ising::h_sigma"        : "1.0",
                        "model::seed"                                 : str(num_total),
                        "xdmrg::num_sites"                            : str(num_L),
			            "xdmrg::chi_max"                              : "64",
                        "fdmrg::num_sites"                            : str(num_L),
                        "fdmrg::chi_max"                              : "64",
                        "hdf5::output_folder"                         : 'output/L_'+ str(num_L) + '/l_'+str(num_l) + '/J_' +str(num_j) + '/h_'+ str(num_h),
                        "hdf5::output_filename"                       : basename + str(num_r) + '.h5'
                    })
                    num_total  = num_total + 1

generate_input_files(settings, input_filenames, template_filename,location)



