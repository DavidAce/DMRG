from generate_inputs.src.generate_inputs import *
import numpy as np
from fractions import Fraction
from itertools import product

# Using the input_template.cfg, make num_copies new files enumerated as 0....num_copies,
# while replacing the fields stated in find_replace.

template_filename = 'input_flbit_template.cfg'
basename    = 'mbl'
location    = "input"


sites               = [12]
J                   = [[0.000, 0.200, 0.000]]
w                   = [[1.000, 0.200, 0.040]]
f                   = [0.01, 0.20,0.30, 0.60]
u                   = [6]
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
    str_f = "{:.2f}".format(val_f)
    str_u = "{}".format(val_u)
    str_J1,str_J2, str_J3 = "{:+.2f}".format(val_J[0]),"{:+.2f}".format(val_J[1]),"{:+.2f}".format(val_J[2])
    str_w1,str_w2, str_w3 = "{:+.2f}".format(val_w[0]),"{:+.2f}".format(val_w[1]),"{:+.2f}".format(val_w[2])

    input_filename = "{}/{}/L_{}/{}/{}/f_{}/u_{}.cfg".format(location,basename,str_L,str_J,str_w, str_f, str_u)

    settings = {
        "output::output_filepath"            : "L_{}/{}/{}/f_{}/u_{}/{}.h5".format(output_prefix,str_L,str_J,str_w, str_f, str_u, basename),
        "console::verbosity"                 : "2",
        "strategy::initial_state"            : str(init),
        "model::model_size"                  : str_L,
        "model::lbit::J1_mean"               : str_J1,
        "model::lbit::J2_mean"               : str_J2,
        "model::lbit::J3_mean"               : str_J3,
        "model::lbit::J1_wdth"               : str_w1,
        "model::lbit::J2_wdth"               : str_w2,
        "model::lbit::J3_wdth"               : str_w3,
        "model::lbit::J3_wdth"               : str_w3,
        "model::lbit::f_mixer"               : str_f,
        "model::lbit::u_layer"               : str_u,
        "flbit::chi_lim_max"                 : "256",
        "flbit::time_start_real"             : "1e-1",
        "flbit::time_start_imag"             : "0",
        "flbit::time_final_real"             : "1e6",
        "flbit::time_final_imag"             : "0",
        "flbit::time_num_steps"              : "150",
    }
    os.makedirs(location, exist_ok=True)
    num_total = num_total + 1
    print(input_filename, "L:", str_L,str_J, str_w, "f:", str_f)
    generate_input_file(settings, input_filename, template_filename)
