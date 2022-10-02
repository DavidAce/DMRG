from generate_inputs.src.generate_inputs import *
import numpy as np
from fractions import Fraction
from itertools import product
import platform

# Using the input_template.cfg, make num_copies new files enumerated as 0....num_copies,
# while replacing the fields stated in find_replace.

template_filename = 'input_template_flbit.cfg'
basename    = 'mbl'
location    = "input"


sites               = [8,12,16,20]
J                   = [[0.000, 0.000, 0.000]]
w                   = [[1.000, 0.250, 0.100]] # for w2, nearest neighbors have this order of magnitude
f                   = [0.2500]
u                   = [4]
x                   = [0.8]
r                   = [-1]
initial_state       = ["PRODUCT_STATE_NEEL"]
output_prefix       = "output"
tmp_storage = "/tmp"
if "lith" in  platform.node():
    tmp_storage = "/scratch/local"


num_total = 0
settings = []
input_filenames = []

print("Generating", len(sites) * len(J) * len(w) * len(f) * len(initial_state), "input files")

for val_L,val_J,val_w, val_x, val_f,val_u, val_r, init, in  product(sites,J,w, x,f, u, r, initial_state):
    str_L = str(val_L)
    str_J = "J[{:+.3f}_{:+.3f}_{:+.3f}]".format(val_J[0], val_J[1], val_J[2])
    str_w = "w[{:.3f}_{:.3f}_{:.3f}]".format(val_w[0], val_w[1], val_w[2])
    str_x = "{:.3f}".format(val_x)
    str_f = "{:.3f}".format(val_f)
    str_u = "{}".format(val_u)
    str_r = "{}".format(val_r)
    str_J1,str_J2, str_J3 = "{:+.3f}".format(val_J[0]),"{:+.3f}".format(val_J[1]),"{:+.3f}".format(val_J[2])
    str_w1,str_w2, str_w3 = "{:+.3f}".format(val_w[0]),"{:+.3f}".format(val_w[1]),"{:+.3f}".format(val_w[2])

    input_filename = "{}/{}_L{}_{}_{}_x{}_f{}_u{}_r{}.cfg".format(location,basename,str_L,str_J,str_w, str_x,str_f, str_u,str_r)

    settings = {
        "storage::output_filepath"           : "{}/L_{}/{}/{}/x_{}/f_{}/u_{}/r_{}/{}.h5".format(output_prefix,str_L,str_J,str_w,str_x, str_f, str_u,str_r, basename),
        "storage::temp_dir"                  : tmp_storage,
        "console::loglevel"                  : "2",
        "strategy::initial_state"            : str(init),
        "model::model_size"                  : str_L,
        "model::lbit::J1_mean"               : str_J1,
        "model::lbit::J2_mean"               : str_J2,
        "model::lbit::J3_mean"               : str_J3,
        "model::lbit::J1_wdth"               : str_w1,
        "model::lbit::J2_wdth"               : str_w2,
        "model::lbit::J3_wdth"               : str_w3,
        "model::lbit::J3_wdth"               : str_w3,
        "model::lbit::J2_xcls"               : str_x,
        "model::lbit::J2_span"               : str_r,
        "model::lbit::f_mixer"               : str_f,
        "model::lbit::u_layer"               : str_u,
        "flbit::bond_max"                    : "2048",
        "flbit::time_start_real"             : "1e-1",
        "flbit::time_start_imag"             : "0",
        "flbit::time_final_real"             : "1e12",
        "flbit::time_final_imag"             : "0",
        "flbit::time_num_steps"              : "100",
    }
    os.makedirs(location, exist_ok=True)
    num_total = num_total + 1
    print(input_filename, "L:", str_L,str_J, str_w, "x:", str_x, "f:", str_f,"u:", str_u, "r:",str_r)
    generate_input_file(settings, input_filename, template_filename)
