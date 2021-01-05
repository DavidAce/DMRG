from generate_inputs.src.generate_inputs import *
import numpy as np
from fractions import Fraction
from itertools import product

# Using the input_template.cfg, make num_copies new files enumerated as 0....num_copies,
# while replacing the fields stated in find_replace.

template_filename = 'input_flbit_template.cfg'
basename    = 'mbl'
location    = "input"


sites               = [16,24,32]
J                   = [[0.000, 0.500, 0.000]]
w                   = [[0.500, 0.500, 0.500]]
f                   = [0.1]
initial_state       = ["PRODUCT_STATE_NEEL"]
output_prefix       = "output"

num_total = 0
settings = []
input_filenames = []

print("Generating", len(sites) * len(J) * len(w) * len(f) * len(initial_state), "input files")

for val_L,val_J,val_w, val_f, init, in  product(sites,J,w, f, initial_state):
    str_L = str(val_L)
    str_J = "J[{:+.2f}_{:+.2f}_{:+.2f}]".format(val_J[0], val_J[1], val_J[2])
    str_w = "w[{:+.2f}_{:+.2f}_{:+.2f}]".format(val_w[0], val_w[1], val_w[2])
    str_f = "{:+.2f}".format(val_f)
    str_J1,str_J2, str_J3 = "{:+.2f}".format(val_J[0]),"{:+.2f}".format(val_J[1]),"{:+.2f}".format(val_J[2])
    str_w1,str_w2, str_w3 = "{:+.2f}".format(val_w[0]),"{:+.2f}".format(val_w[1]),"{:+.2f}".format(val_w[2])

    input_filename = location + '/' + basename + '_L'+ str_L + '_' + str_J + '_' + str_w + '_f' + str_f + '.cfg'

    settings = {
        "output::output_filepath"            : output_prefix + '/L_'+ str_L + '/' + str_J + '/' + str_w + '/f' + str_f + '/' + basename + '.h5',
        "threading::num_threads"             : "2",
        "console::verbosity"                 : "2",
        "strategy::initial_state"            : str(init),
        "model::model_size"                  : str_L,
        "model::lbit::J1"                    : str_J1,
        "model::lbit::J2"                    : str_J2,
        "model::lbit::J3"                    : str_J3,
        "model::lbit::w1"                    : str_w1,
        "model::lbit::w2"                    : str_w2,
        "model::lbit::w3"                    : str_w3,
        "flbit::chi_lim_max"                 : "256",
        "flbit::time_step_init_real"         : "1e-2",
        "flbit::time_step_init_imag"         : "0.0",
        "flbit::time_step_max_size"          : "1e4",
        "flbit::time_step_per_size"          : "10",
        "flbit::time_step_growth_factor"     : "10",
        "flbit::time_limit"                  : "1e6",
    }
    os.makedirs(location, exist_ok=True)
    num_total = num_total + 1
    print(input_filename, "L:", str_L,str_J, str_w, "f:", str_f)
    generate_input_file(settings, input_filename, template_filename)
