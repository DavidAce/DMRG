from generate_inputs.src.generate_inputs import *
import numpy as np
from fractions import Fraction
from itertools import product
import platform

# Using the input_template.cfg, make num_copies new files enumerated as 0....num_copies,
# while replacing the fields stated in find_replace.

template_filename = 'input_template_flbit.cfg'
basename          = 'mbl'
input_prefix      = "input-utest"
output_prefix     = "output-utest"


sites               = [8,12,16,20,24,28,32]
J                   = [[0.00, 0.00, 0.00]]
w                   = [[1.00, 1.00, 1.00]] # for w2, nearest neighbors have this order of magnitude
x                   = [1.0]
r                   = [-1]
u_depth             = [8,12,16,20,24,28,32]
u_fmix              = [1.0]
u_tstd              = [1.0]
u_cstd              = [1.0]
u_tgw8              = ['IDENTITY']
u_cgw8              = ['EXPDECAY']
u_bond              = [64,128,192,256]
# initial_state       = ["PRODUCT_STATE_NEEL"]
initial_state       = ["RANDOM_PRODUCT_STATE_ZEROMAG"]
tmp_storage = "/tmp"
if "lith" in  platform.node():
    tmp_storage = "/scratch/local"


num_total = 0
settings = []
input_filenames = []
batch_size = len(sites) * len(J) * len(w) * len(x) * len(r) * len(u_depth) * len(u_fmix) * len(u_tstd) * len(u_cstd) * len(u_tgw8) * len(u_cgw8) * len(initial_state)
print(f"Generating {batch_size} input files")

for val_L,val_J,val_w, val_x, val_r, val_u, val_f, val_tstd, val_cstd, val_tgw8, val_cgw8, val_bond, init in product(sites,J,w, x, r, u_depth, u_fmix, u_tstd, u_cstd, u_tgw8, u_cgw8, u_bond, initial_state):
    str_L = str(val_L)
    str_J = "[{:+.2f}±{:.2f}_{:+.2f}±{:.2f}_{:+.2f}±{:.2f}]".format(val_J[0], val_w[0],val_J[1], val_w[1], val_J[2], val_w[2])
    # str_J1,str_J2, str_J3 = "{:+.2f}".format(val_J[0]),"{:+.2f}".format(val_J[1]),"{:+.2f}".format(val_J[2])
    # str_w1,str_w2, str_w3 = "{:+.2f}".format(val_w[0]),"{:+.2f}".format(val_w[1]),"{:+.2f}".format(val_w[2])
    str_J1,str_J2, str_J3 = ["{:+.2f}".format(val) for val in val_J]
    str_w1,str_w2, str_w3 = ["{:+.2f}".format(val) for val in val_w]
    str_x = "{:.2f}".format(val_x)
    str_r = "{}".format(val_r)
    str_rL = "{}".format(val_r) if val_r >= 0 else "L"
    str_u = "{}".format(val_u)
    str_f = "{:.2f}".format(val_f)
    str_tstd = "{:.2f}".format(val_tstd)
    str_cstd = "{:.2f}".format(val_cstd)
    str_tgw8 = "{}".format(val_tgw8)
    str_cgw8 = "{}".format(val_cgw8)
    str_bond = "{}".format(val_bond)


    # Generate a unique identifier for the config file and output directory
    # Hamiltonian part:
    input_filename =  f"{input_prefix}/{basename}_L{str_L}_J{str_J}_x{str_x}_r{str_rL}"
    output_filedir =   f"{output_prefix}/L{str_L}/J{str_J}/x{str_x}/r{str_rL}"
    # Unitary circuit string:
    str_circuit = f"d{str_u}_f{str_f}"
    str_circuit += (f"_tw{str_tstd}{str_tgw8[:2]}" if (len(u_tstd) > 1 or len(u_tgw8) > 1) else '')
    str_circuit += (f"_cw{str_cstd}{str_cgw8[:2]}" if (len(u_cstd) > 1 or len(u_cgw8) > 1) else '')
    str_circuit += (f"_bond{str_bond}" if (len(u_bond) > 1) else '')
    input_filename += f"_u[{str_circuit}].cfg"
    output_filedir += f"/u[{str_circuit}]"

    settings = {
        "storage::output_filepath"           : f"{output_filedir}/{basename}.h5",
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
        "model::lbit::J2_span"               : str_r,
        "model::lbit::xi_Jcls"               : str_x,
        "model::lbit::u_depth"               : str_u,
        "model::lbit::u_fmix"                : str_f,
        "model::lbit::u_tstd"                : str_tstd,
        "model::lbit::u_cstd"                : str_cstd,
        "model::lbit::u_tgw8"                : str_tgw8,
        "model::lbit::u_cgw8"                : str_cgw8,
        "flbit::bond_max"                    : "2048",
        "flbit::time_start_real"             : "1e-1",
        "flbit::time_start_imag"             : "0",
        "flbit::time_final_real"             : "1e10",
        "flbit::time_final_imag"             : "0",
        "flbit::time_num_steps"              : "500",
        "flbit::cls::mpo_circuit_svd_bondlim": str_bond,
    }
    os.makedirs(input_prefix, exist_ok=True)
    num_total = num_total + 1
    print(input_filename)
    generate_input_file(settings, input_filename, template_filename)
