from utils.generators import get_config_product, write_config_file, write_batch_files, move_directories
from utils.flbit import get_max_time, get_max_steps, get_output_filepath, get_config_filename, update_batch_status
from batch_setup import get_batch_setup
import os
import platform

config_paths = {
    'config_template'   : 'template_configs/flbit-v2.cfg',
    'output_prfx'       : "/mnt/WDB-AN1500/mbl_transition",
    'output_stem'       : 'mbl',
    'config_dir'        : "config-lbit123",
    'output_dir'        : "output-lbit123",
    'status_dir'        : "status-lbit123",
    'temp_dir'          : "/scratch/local" if "lith" in platform.node() else (os.environ.get('PDC_TMP') if "PDC_TMP" in os.environ else "/tmp")
}

config_ranges = {
    "filename" : [''],
    "storage::output_filepath": [get_output_filepath],
    "storage::temp_dir": [config_paths['temp_dir']],
    "storage::copy_from_temp_freq": ['5'],
    "console::loglevel": ['2'],
    "solver::svd_truncation_lim": ['1e-5'],
    "solver::svd_truncation_init": ['1e-5'],
    "solver::svd_switchsize_bdc": ['16'],
    "solver::svd_save_fail": ['false'],
    "strategy::initial_state": ["PRODUCT_STATE_NEEL"],
    "model::model_type": ['lbit'],
    "model::model_size": ['12', '16', '20'],
    "model::lbit::J1_mean": ['+0.00'],
    "model::lbit::J2_mean": ['+0.00'],
    "model::lbit::J3_mean": ['+0.00'],
    "model::lbit::J1_wdth": ['1.00'],
    "model::lbit::J2_wdth": ['1.00'],
    "model::lbit::J3_wdth": ['1.00'],
    "model::lbit::J2_span": ['-1'],
    "model::lbit::xi_Jcls": ['1.00'],
    "model::lbit::u_depth": ['8'],
    "model::lbit::u_fmix": ['0.20'],
    "model::lbit::u_lambda": ['0.1', '1.0', '5.0'],
    "model::lbit::u_wkind": ['EXPDECAY'],
    "model::lbit::u_mkind": ['MATRIX_V2'],
    "flbit::run_effective_model": ['true'],
    "flbit::time_scale": ['LOGSPACED'],
    "flbit::bond_max": ['8192'],
    "flbit::time_start_real": ['1e-1'],
    "flbit::time_start_imag": ['0.0'],
    "flbit::time_final_real": [get_max_time],  # "{:.1e}".format(max_time),
    "flbit::time_final_imag": ['0.0'],
    "flbit::time_num_steps": [get_max_steps],
    "flbit::cls::num_rnd_circuits": ['1'],
    "flbit::cls::mpo_circuit_svd_bondlim": ['20'],
    "flbit::opdm::num_rps": ['10'],
}

configs = get_config_product(config_ranges, config_paths)
for config in configs:
    # Write the config file
    config['filename'] = get_config_filename(config, config_ranges, config_paths)
    config_template = config_paths['config_template']
    write_config_file(config, config_template, config['filename'])

batch_setup = get_batch_setup('lbit123-fleff')
write_batch_files(batch_setup=batch_setup, configs=configs, config_paths=config_paths)
update_batch_status(config_paths=config_paths)
move_directories(batch_setup=batch_setup, config_paths=config_paths)