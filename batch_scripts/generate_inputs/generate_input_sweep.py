from src.generate_inputs import *
import numpy as np


# Using the input_template.cfg, make num_copies new files enumerated as 0....num_copies,
# while replacing the fields stated in find_replace.

template_filename = 'input_template.cfg'
basename    = 'mbl'
location    = "../input"
os.makedirs(location, exist_ok=True)


sites        = np.linspace(16,28,4, dtype=int)
lambdas      = [0] # np.linspace(0,0.2,3)
deltas       = [0] # np.linspace(-1.0,1.0,5)
J_log_mean   = np.array([1])
h_log_mean   = J_log_mean - deltas
# h_log_mean   = np.flipud(np.arange(-1, 1.5, 0.5)+1)
num_total = 0
settings = []
input_filenames = []
print("Generating", len(sites) * len(lambdas) * len(J_log_mean) * len(h_log_mean) , "input files")
for num_L in sites:
    for num_l in range(len(lambdas)):
        for num_j in range(len(J_log_mean)):
            for num_h in range(len(h_log_mean)):
                input_filename_dir = location + '/' + 'L_'+ str(num_L)
                os.makedirs(input_filename_dir, exist_ok=True)
                input_filename = input_filename_dir + '/' + basename + '_l' + str(num_l) + '_J'+ str(num_j) + '_h'+ str(num_h) + '.cfg'
                settings = {
                    "model::selfdual_tf_rf_ising::J_log_mean"     : "{:.2f}".format(J_log_mean[num_j]),
                    "model::selfdual_tf_rf_ising::h_log_mean"     : "{:.2f}".format(h_log_mean[num_h]),
                    "model::selfdual_tf_rf_ising::lambda"         : "{:.2f}".format(lambdas[num_l]),
                    "model::selfdual_tf_rf_ising::J_sigma"        : "1.0",
                    "model::selfdual_tf_rf_ising::h_sigma"        : "1.0",
                    "model::seed"                                 : str(num_total),
                    "xdmrg::num_sites"                            : str(num_L),
                    "xdmrg::chi_max"                              : "128",
                    "hdf5::output_filename"                       : 'output/L_'+ str(num_L) + '/l_'+str(num_l) + '/J_' +str(num_j) + '/h_'+ str(num_h)+ '/' + basename + '.h5'
                }
                num_total = num_total + 1
                generate_input_file(settings, input_filename, template_filename)

# generate_input_files(settings, input_filenames, template_filename,location)



