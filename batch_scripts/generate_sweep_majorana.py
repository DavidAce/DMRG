from generate_inputs.src.generate_inputs import *
import numpy as np
from fractions import Fraction
from itertools import product
import platform

# Using the input_template.cfg, make num_copies new files enumerated as 0....num_copies,
# while replacing the fields stated in find_replace.

template_filename = 'input_template_majorana.cfg'
basename    = 'mbl'
location    = "input"


sites               = np.array([10])
gs                  = [0.000, 0.01, 0.5]
deltas              = [-8.0, -7.0, -6.0, -5.0, -4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]
initial_state       = ["RANDOM_PRODUCT_STATE"]
multisite_mps_size_def  = [1]
multisite_mps_size_max  = [2]
output_prefix       = "output"

tmp_storage = "/tmp"
if "lith" in  platform.node():
    tmp_storage = "/scratch/local"

num_total = 0
settings = []
input_filenames = []

print("Generating", len(sites) * len(gs) * len(deltas) * len(multisite_mps_size_def) * len(initial_state), "input files")

for val_L,val_g, val_d, init, multi in  product(sites,gs,deltas,initial_state,multisite_mps_size_max):
    str_L = str(val_L)
    str_d = "{:+.4f}".format(val_d)
    str_g = "{:.4f}".format(val_g)
    str_j = "{:.6f}".format(val_j)
    str_h = "{:.6f}".format(val_h)

    extra_prefix       = ""
    if len(initial_state) > 1:
        if init == "RANDOM_PRODUCT_STATE":
            extra_prefix = extra_prefix + "_rps"
        else:
            extra_prefix = extra_prefix + "_res"

    if len(multisite_mps_size_def) > 1:
            extra_prefix = extra_prefix + "_multi" + str(multi)

    input_filename = "{}/{}_L{}_g{}_d{}.cfg".format(location+extra_prefix,basename,str_L,str_g,str_d)
    settings = {
        "storage::output_filepath"            : "{}/L_{}/l_{}/d_{}/{}.h5".format(output_prefix+extra_prefix,str_L,str_g,str_d, basename),
        "storage::temp_dir"                  : tmp_storage,
        "threading::stl_threads"             : "2",
        "threading::omp_threads"             : "2",
        "console::loglevel"                  : "2",
        "model::model_size"                  : str_L,
        "model::ising_majorana::delta"       : str_d,
        "model::ising_majorana::g"           : str_g,
        "xdmrg::chi_lim_max"                 : "1024",
        "xdmrg::max_states"                  : "1",
        "strategy::multisite_mps_size_def"   : str(multisite_mps_size_def[0]),
        "strategy::multisite_mps_size_max"   : str(multi),
        "strategy::initial_state"            : str(init),
    }
    os.makedirs(location + extra_prefix, exist_ok=True)
    num_total = num_total + 1
    print(input_filename, "L:", str_L,"g:",str_g, "d:", str_d)
    generate_input_file(settings, input_filename, template_filename)
