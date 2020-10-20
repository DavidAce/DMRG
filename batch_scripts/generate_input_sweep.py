from generate_inputs.src.generate_inputs import *
import numpy as np
from fractions import Fraction

# Using the input_template.cfg, make num_copies new files enumerated as 0....num_copies,
# while replacing the fields stated in find_replace.

template_filename = 'input_template.cfg'
basename    = 'mbl'
location    = "input"
os.makedirs(location, exist_ok=True)


sites        = np.array([16,24])
lambdas      = [0.005,0.010, 0.015,0.020]
deltas       = [-0.1, -0.05, 0.0, 0.05, 0.1]

num_total = 0
settings = []
input_filenames = []

def delta(J,h):
    return np.log(J) - np.log(h)

def undelta(delta):
    exp_delta = np.exp(delta)
    fr = Fraction(exp_delta).limit_denominator()
    maxJh = np.max([fr.numerator, fr.denominator])
    J = round(fr.numerator/maxJh,6)
    h = round(fr.denominator/maxJh,6)
    return J, h


print("Generating", len(sites) * len(lambdas) * len(deltas), "input files")
for num_L,val_L in enumerate(sites):
    for num_l,val_l in enumerate(lambdas):
        for num_d,val_d in enumerate(deltas):
            val_j,val_h = undelta(val_d)
            os.makedirs(location, exist_ok=True)
            str_L = str(val_L)
            str_d = "{:+.4f}".format(val_d)
            str_l = "{:.4f}".format(val_l)
            str_j = "{:.6f}".format(val_j)
            str_h = "{:.6f}".format(val_h)

            input_filename = location + '/' + basename + '_L'+ str_L + '_l' + str_l + '_d'+ str_d + '.cfg'
            settings = {
                "output::output_filepath"            : 'output/L_'+ str_L + '/l_'+ str_l + '/d_' + str_d + '/' + basename + '.h5',
                "threading::num_threads"             : "2",
                "model::model_size"                  : str_L,
                "model::ising_sdual::lambda"         : str_l,
                "model::ising_sdual::J_mean"         : str_j,
                "model::ising_sdual::h_mean"         : str_h,
                "model::ising_sdual::J_stdv"         : "1.0",
                "model::ising_sdual::h_stdv"         : "1.0",
                "xdmrg::chi_lim_max"                 : "512",
                "xdmrg::max_states "                 : "2",
            }
            num_total = num_total + 1
            print(input_filename, "L:", str_L,"l:",str_l, "d:", "{:+6.4f}".format(val_d),"j/h:", undelta(val_d))
            generate_input_file(settings, input_filename, template_filename)




