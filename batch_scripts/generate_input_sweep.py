from generate_inputs.src.generate_inputs import *
import numpy as np


# Using the input_template.cfg, make num_copies new files enumerated as 0....num_copies,
# while replacing the fields stated in find_replace.

template_filename = 'input_template.cfg'
basename    = 'mbl'
location    = "input"
os.makedirs(location, exist_ok=True)


sites        = np.array([16,20,24]) #np.linspace(16,36,6, dtype=int)
# sites        = np.linspace(16,36,6, dtype=int)
lambdas      = [0.01] # np.linspace(0,0.2,3)
deltas       = [0] # np.linspace(-1.0,1.0,5)
J_mean   = np.array([1])
h_mean   = J_mean - deltas
# h_log_mean   = np.flipud(np.arange(-1, 1.5, 0.5)+1)
num_total = 0
settings = []
input_filenames = []
print("Generating", len(sites) * len(lambdas) * len(J_mean) * len(h_mean) , "input files")
for num_L in sites:
    for num_l in range(len(lambdas)):
        for num_j in range(len(J_mean)):
            for num_h in range(len(h_mean)):
                os.makedirs(location, exist_ok=True)
                input_filename = location + '/' + basename + '_L'+ str(num_L) + '_l' + str(num_l) + '_J'+ str(num_j) + '_h'+ str(num_h) + '.cfg'
                settings = {
                    "output::output_filepath"            : 'output/L_'+ str(num_L) + '/l_'+str(num_l) + '/J_' +str(num_j) + '/h_'+ str(num_h)+ '/' + basename + '.h5',
                    "threading::num_threads"             : "1",
                    "model::model_size"                  : str(num_L),
                    "model::ising_sdual::J_mean"         : "{:.2f}".format(J_mean[num_j]),
                    "model::ising_sdual::h_mean"         : "{:.2f}".format(h_mean[num_h]),
                    "model::ising_sdual::lambda"         : "{:.2f}".format(lambdas[num_l]),
                    "model::ising_sdual::J_stdv"         : "1.0",
                    "model::ising_sdual::h_stdv"         : "1.0",
                    "xdmrg::chi_lim_max"                 : "512",
                    "xdmrg::max_states "                 : "2",
                }
                num_total = num_total + 1
                generate_input_file(settings, input_filename, template_filename)

# generate_input_files(settings, input_filenames, template_filename,location)



