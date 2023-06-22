from utils.generators import get_config_product, write_config_file, write_seed_files
from utils.flbit import get_max_time, get_output_filepath, get_config_filename
from seed_setup import get_seed_setup
import platform

config_paths = {
    'project_prefix'    : 'mbl',
    'config_template'   : 'template_configs/flbit.cfg',
    'config_dir'        : "config-L[28]",
    'output_dir'        : "output",
    'seed_dir'          : "config-L[24]", # Output one .json file per .cfg, describing the seed values to run
    'temp_dir'          : "/scratch/local" if "lith" in platform.node() else "/tmp"
}

config_ranges = {
    "storage::output_filepath": [get_output_filepath],
    "storage::temp_dir": [config_paths['temp_dir']],
    "console::loglevel": ['2'],
    "solver::svd_truncation_lim": ['1e-4'],
    "solver::svd_truncation_init": ['1e-4'],
    "solver::svd_switchsize_bdc": ['16'],
    "solver::svd_save_fail": ['false'],
    "strategy::initial_state": ["PRODUCT_STATE_NEEL"],
    "model::model_size": ['28'],
    "model::lbit::J1_mean": ['+0.00'],
    "model::lbit::J2_mean": ['+0.00'],
    "model::lbit::J3_mean": ['+0.00'],
    "model::lbit::J1_wdth": ['1.00'],
    "model::lbit::J2_wdth": ['1.00'],
    "model::lbit::J3_wdth": ['1.00'],
    "model::lbit::J2_span": ['-1'],
    "model::lbit::xi_Jcls": ['1.00'],
    "model::lbit::u_depth": ['16'],
    "model::lbit::u_fmix": ['0.20','0.30','0.40'],
    "model::lbit::u_tstd": ['1.0'],
    "model::lbit::u_cstd": ['1.0'],
    "model::lbit::u_tgw8": ['IDENTITY'],
    "model::lbit::u_cgw8": ['EXPDECAY'],
    "flbit::bond_max": ['8192'],
    "flbit::time_start_real": ['1e-1'],
    "flbit::time_start_imag": ['0.0'],
    "flbit::time_final_real": [get_max_time], #"{:.1e}".format(max_time),
    "flbit::time_final_imag": ['0.0'],
    "flbit::time_num_steps": ['100'],
    "flbit::cls::mpo_circuit_svd_bondlim": ['20'],
}

configs = get_config_product(config_ranges, config_paths)
for config in configs:
    # Write the config file
    config_filename = get_config_filename(config, config_ranges, config_paths)
    config_template = config_paths['config_template']
    write_config_file(config, config_template, config_filename)

seed_setup = get_seed_setup('lbit93-precision')
write_seed_files(seed_setup, config_paths)