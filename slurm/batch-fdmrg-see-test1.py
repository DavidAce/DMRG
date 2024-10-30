from utils.generators import get_config_product, write_config_file, write_batch_files, move_directories
from utils.fdmrg import get_output_filepath, get_config_filename, update_batch_status
from batch_setup import get_batch_setup
import os
import platform

config_paths = {
    'config_template'   : 'template_configs/fdmrg-ising-majorana.cfg',
    'output_prfx'       : "/mnt/WDB-AN1500/mbl_transition",
    'output_stem'       : 'mbl',
    'config_dir'        : "config",
    'output_dir'        : "output",
    'status_dir'        : "status",
    'temp_dir'          : "/scratch/local" if "lith" in platform.node() else (os.environ.get('PDC_TMP') if "PDC_TMP" in os.environ else "/tmp")
}

config_ranges = {
    "filename" : [''],
    "storage::output_filepath": [get_output_filepath],
    "storage::temp_dir": [config_paths['temp_dir']],
    "console::loglevel": ['2'],
    "solver::svd_truncation_lim": ['1e-8'],
    "solver::svd_truncation_init": ['1e-8'],
    "solver::svd_switchsize_bdc": ['16'],
    "strategy::initial_state": ["PRODUCT_STATE_NEEL"],
    "model::model_type": ['ising_majorana'],
    "model::model_size": ['16'],
    "model::ising_majorana::g": ['0.50'],
    "model::ising_majorana::delta": ['-8.00', '0.00', '8.00'],
    "fdmrg::ritz": ['SR'],
    "fdmrg::max_iters": ['50'],
    "fdmrg::min_iters": ['1'],
    "fdmrg::warmup_iters": ['2'],
    "fdmrg::bond_max": ['2048'],
    "fdmrg::bond_init": ['32'],
    "fdmrg::print_freq": ['1'],
    "fdmrg::store_wavefn": ['true'],

}

configs = get_config_product(config_ranges, config_paths)
for config in configs:
    # Set up the config file
    config['filename'] = get_config_filename(config, config_ranges, config_paths)
    config['template'] = config_paths['config_template']

batch_setup = get_batch_setup('fdmrg-see-test1')
write_batch_files(batch_setup=batch_setup, configs=configs, config_paths=config_paths)
update_batch_status(config_paths=config_paths)
move_directories(batch_setup=batch_setup, config_paths=config_paths)