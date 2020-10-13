from generate_inputs.src.generate_inputs import *
import numpy as np


# Using the input_template.cfg, make num_copies new files enumerated as 0....num_copies,
# while replacing the fields stated in find_replace.

template_filename = 'input_template.cfg'
basename    = 'mbl'
location    = "input"
os.makedirs(location, exist_ok=True)


sites        = np.array([40,48,56,64]) #np.linspace(16,36,6, dtype=int)
# sites        = np.linspace(16,36,6, dtype=int)
lambdas      = [0.01,0.02,0.03] # np.linspace(0,0.2,3)
deltas       = [0] # np.linspace(-1.0,1.0,5)
J_mean   = np.array([1])
h_mean   = J_mean - deltas
# h_log_mean   = np.flipud(np.arange(-1, 1.5, 0.5)+1)
num_total = 0
settings = []
input_filenames = []
print("Generating", len(sites) * len(lambdas) * len(J_mean) * len(h_mean) , "input files")
for num_L,val_L in enumerate(sites):
    for num_l,val_l in enumerate(lambdas):
        for num_j,val_j in enumerate(J_mean):
            for num_h,val_h in enumerate(h_mean):
                os.makedirs(location, exist_ok=True)
                str_L = str(val_L)
                str_l = "{:.3f}".format(val_l)
                str_j = "{:.3f}".format(val_j)
                str_h = "{:.3f}".format(val_h)

                input_filename = location + '/' + basename + '_L'+ str_L + '_l' + str_l + '_J'+ str_j + '_h'+ str_h + '.cfg'
                settings = {
                    "output::output_filepath"            : 'output/L_'+ str_L + '/l_'+ str_l + '/J_' + str_j + '/h_'+ str_h + '/' + basename + '.h5',
                    "threading::num_threads"             : "1",
                    "model::model_size"                  : str_L,
                    "model::ising_sdual::J_mean"         : str_j,
                    "model::ising_sdual::h_mean"         : str_h,
                    "model::ising_sdual::lambda"         : str_l,
                    "model::ising_sdual::J_stdv"         : "1.0",
                    "model::ising_sdual::h_stdv"         : "1.0",
                    "xdmrg::chi_lim_max"                 : "320",
                    "xdmrg::max_states "                 : "4",
                }
                num_total = num_total + 1
                generate_input_file(settings, input_filename, template_filename)




