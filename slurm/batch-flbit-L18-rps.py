from utils.generators import get_config_product, write_config_file, write_batch_files, move_directories
from utils.flbit import get_max_time, get_output_filepath, get_config_filename, update_batch_status
from batch_setup import get_batch_setup
import os
import platform

config_paths = {
    'config_template'   : 'template_configs/flbit.cfg',
    'output_prfx'       : "/mnt/WDB-AN1500/mbl_transition",
    'output_stem'       : 'mbl',
    'config_dir'        : "config-L18-rps",
    'output_dir'        : "output-rps",
    'status_dir'        : "status-rps",
    'temp_dir'          : "/scratch/local" if "lith" in platform.node() else (os.environ.get('PDC_TMP') if "PDC_TMP" in os.environ else "/tmp")
}

config_ranges = {
    "filename" : [''],
    "storage::output_filepath": [get_output_filepath],
    "storage::temp_dir": [config_paths['temp_dir']],
    "console::loglevel": ['2'],
    "solver::svd_truncation_lim": ['1e-5'],
    "solver::svd_truncation_init": ['1e-5'],
    "solver::svd_switchsize_bdc": ['16'],
    "solver::svd_save_fail": ['false'],
    "strategy::initial_state": ["PRODUCT_STATE_NEEL_SHUFFLED"],
    "model::model_size": ['18'],
    "model::lbit::J1_mean": ['+0.00'],
    "model::lbit::J2_mean": ['+0.00'],
    "model::lbit::J3_mean": ['+0.00'],
    "model::lbit::J1_wdth": ['1.00'],
    "model::lbit::J2_wdth": ['1.00'],
    "model::lbit::J3_wdth": ['1.00'],
    "model::lbit::J2_span": ['-1'],
    "model::lbit::xi_Jcls": ['1.00'],
    "model::lbit::u_depth": ['16'],
    "model::lbit::u_fmix": ['0.20','0.40'],
    "model::lbit::u_tstd": ['1.0'],
    "model::lbit::u_cstd": ['1.0'],
    "model::lbit::u_tgw8": ['EXPDECAY'],
    "model::lbit::u_cgw8": ['EXPDECAY'],
    "flbit::bond_max": ['8192'],
    "flbit::time_scale": ['LOGSPACED'],
    "flbit::time_start_real": ['1e-1'],
    "flbit::time_start_imag": ['0.0'],
    "flbit::time_final_real": [get_max_time], #"{:.1e}".format(max_time),
    "flbit::time_final_imag": ['0.0'],
    "flbit::time_num_steps": ['200'],
    "flbit::cls::mpo_circuit_svd_bondlim": ['20'],
    "flbit::cls::num_rnd_circuits": ['1']  # Enables l-bit correlation matrix O(i,j) calculation
}

configs = get_config_product(config_ranges, config_paths)
for config in configs:
    # Set up the config file
    config['filename'] = get_config_filename(config, config_ranges, config_paths)
    config['template'] = config_paths['config_template']

batch_setup = get_batch_setup('lbit101-rps')
write_batch_files(batch_setup=batch_setup, configs=configs, config_paths=config_paths)
update_batch_status(config_paths=config_paths)
move_directories(batch_setup=batch_setup, config_paths=config_paths)