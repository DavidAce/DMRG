from generate_inputs.src.generate_inputs import *
import numpy as np
from fractions import Fraction
from itertools import product

# Using the input_template.cfg, make num_copies new files enumerated as 0....num_copies,
# while replacing the fields stated in find_replace.

template_filename = 'input_template_flbit.cfg'
basename    = 'mbl'
location    = "input"


sites               = [6,8,10,12,14,16]
J                   = [[0.000, 0.000, 0.000]]
w                   = [[1.000, 1.000, 0.100]]
f                   = [0.20, 0.30, 0.40, 0.50]
u                   = [2,3,4]
b                   = [2,3,4,5,6]
r                   = [2,3,4,5,6,7,8]
initial_state       = ["PRODUCT_STATE_NEEL"]
output_prefix       = "output"

num_total = 0
settings = []
input_filenames = []

print("Generating", len(sites) * len(J) * len(w) * len(f) * len(initial_state), "input files")

for val_L,val_J,val_w, val_b, val_f,val_u, val_r, init, in  product(sites,J,w, b,f, u, r, initial_state):
    str_L = str(val_L)
    str_J = "J[{:+.2f}_{:+.2f}_{:+.2f}]".format(val_J[0], val_J[1], val_J[2])
    str_w = "w[{:.2f}_{:.2f}_{:.2f}]".format(val_w[0], val_w[1], val_w[2])
    str_b = "{:.2f}".format(val_b)
    str_f = "{:.2f}".format(val_f)
    str_u = "{}".format(val_u)
    str_r = "{}".format(val_r)
    str_J1,str_J2, str_J3 = "{:+.2f}".format(val_J[0]),"{:+.2f}".format(val_J[1]),"{:+.2f}".format(val_J[2])
    str_w1,str_w2, str_w3 = "{:+.2f}".format(val_w[0]),"{:+.2f}".format(val_w[1]),"{:+.2f}".format(val_w[2])

    input_filename = "{}/{}_L{}_{}_{}_b{}_f{}_u{}_r{}.cfg".format(location,basename,str_L,str_J,str_w, str_b,str_f, str_u,str_r)

    settings = {
        "output::output_filepath"            : "{}/L_{}/{}/{}/b_{}/f_{}/u_{}/r_{}/{}.h5".format(output_prefix,str_L,str_J,str_w,str_b, str_f, str_u,str_r, basename),
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
        "model::lbit::J2_base"               : str_b,
        "model::lbit::J2_span"               : str_r,
        "model::lbit::f_mixer"               : str_f,
        "model::lbit::u_layer"               : str_u,
        "flbit::chi_lim_max"                 : "1024",
        "flbit::time_start_real"             : "1e-1",
        "flbit::time_start_imag"             : "0",
        "flbit::time_final_real"             : "1e8",
        "flbit::time_final_imag"             : "0",
        "flbit::time_num_steps"              : "100",
    }
    os.makedirs(location, exist_ok=True)
    num_total = num_total + 1
    print(input_filename, "L:", str_L,str_J, str_w, "b:", str_b, "f:", str_f,"u:", str_u, "r:",str_r)
    generate_input_file(settings, input_filename, template_filename)
