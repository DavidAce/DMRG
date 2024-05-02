from utils.generators import get_config_product, write_config_file, write_batch_files, move_directories
from utils.xdmrg import get_output_filepath, get_config_filename, update_batch_status
from xdmrg_batches import get_xdmrg_batch_setup
import os
import platform

config_paths = {
    'config_template'   : 'template_configs/xdmrg-ising-majorana.cfg',
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
    "storage::file_collision_policy": ['REVIVE'],
    "storage::temp_dir": [config_paths['temp_dir']],
    "storage::mps::state_emid::policy": ["ITER|FINISH|REPLACE"],
    "storage::table::kvornings_marker::policy": ["FINISH"],
    "storage::dataset::subsystem_entanglement_entropies::bond_lim": ["2048"],
    "storage::dataset::subsystem_entanglement_entropies::trnc_lim": ["1e-6"],
    "console::loglevel": ['2'],
    "solver::svd_truncation_lim": ['1e-8'],
    "solver::svd_truncation_init": ['1e-8'],
    "solver::svd_switchsize_bdc": ['16'],
    "strategy::initial_state": ["PRODUCT_STATE_NEEL"],
    "precision::use_parity_shifted_mpo": ["true"],
    "precision::use_parity_shifted_mpo_squared": ["true"],
    "model::model_type": ['ising_majorana'],
    "model::model_size": ['12'],
    "model::ising_majorana::g": ['0.00', '0.01', '0.02', '0.03'],
    "model::ising_majorana::delta": [f'{x:+.2f}' for x in range(-6,7)],
    "xdmrg::max_iters": ['500'],
    "xdmrg::min_iters": ['1'],
    "xdmrg::warmup_iters": ['4'],
    "xdmrg::bond_max": ['2048'],
    "xdmrg::bond_init": ['32'],
    "xdmrg::print_freq": ['1'],
}

configs = get_config_product(config_ranges, config_paths)
for config in configs:
    # Write the config file
    config['filename'] = get_config_filename(config, config_ranges, config_paths)
    config_template = config_paths['config_template']
    write_config_file(config, config_template, config['filename'])

batch_setup = get_xdmrg_batch_setup('xdmrg1-fse')
write_batch_files(batch_setup=batch_setup, configs=configs, config_paths=config_paths)
update_batch_status(config_paths=config_paths)
move_directories(batch_setup=batch_setup, config_paths=config_paths)