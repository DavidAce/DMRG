from utils.generators import get_config_product, write_config_file, write_seed_files
from utils.flbit import get_max_time, get_output_filepath, get_config_filename
import platform


config_paths = {
    'config_template'   : 'template_configs/flbit.cfg',
    'config_dir'        : "config",
    'output_stem'       : 'mbl',
    'output_dir'        : "output",
    'temp_dir'          : "/scratch/local" if "lith" in platform.node() else "/tmp"
}


config_ranges = {
    "filename" : [''],
    "storage::output_filepath": [get_output_filepath],
    "storage::temp_dir": [config_paths['temp_dir']],
    "console::loglevel": ['2'],
    "solver::svd_truncation_lim": '1e-4',
    "solver::svd_truncation_init": '1e-4',
    "solver::svd_switchsize_bdc": '16',
    "solver::svd_save_fail": 'false',
    "strategy::initial_state": ["PRODUCT_STATE_NEEL"],
    "model::model_size": ['24'],
    "model::lbit::J1_mean": ['0.0'],
    "model::lbit::J2_mean": ['0.0'],
    "model::lbit::J3_mean": ['0.0'],
    "model::lbit::J1_wdth": ['1.0'],
    "model::lbit::J2_wdth": ['1.0'],
    "model::lbit::J3_wdth": ['1.0'],
    "model::lbit::J2_span": ['-1'],
    "model::lbit::xi_Jcls": ['1.0'],
    "model::lbit::u_depth": ['16'],
    "model::lbit::u_fmix": ['0.2','0.3','0.4'],
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

batch_setup = {
    'projectname': 'lbit93-precision',
    'seed_extent': [10000],      # Number of seeds to run
    'seed_offset': [200000000],  # Start the seed count here
    'seed_stride': 100000000,    # Leave this much room between configs
}

configs = get_config_product(config_ranges, config_paths)
for config in configs:
    # Write the config file
    config_filename = get_config_filename(config, config_ranges, config_paths)
    config_template = config_paths['config_template']
    write_config_file(config, config_template, config_filename)

write_seed_files(batch_setup, config_paths)
