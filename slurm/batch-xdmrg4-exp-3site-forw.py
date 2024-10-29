from utils.generators import get_config_product, write_config_file, write_batch_files, move_directories
from utils.xdmrg import get_output_filepath, get_config_filename, update_batch_status
from batches_xdmrg import get_xdmrg_batch_setup
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
    "console::loglevel": ['2'],
    "storage::output_filepath": [get_output_filepath],
    "storage::resume_policy": ['IF_UNSUCCESSFUL'],
    "storage::file_collision_policy": ['REPLACE'],
    "storage::temp_dir": [config_paths['temp_dir']],
    "storage::mps::state_emid::policy": ["ITER|FINISH|REPLACE"],
    "storage::table::opdm::policy": ["FINISH|RBDS|RTES"],
    "storage::table::opdm_spectrum::policy": ["FINISH|RBDS|RTES"],
    "storage::dataset::subsystem_entanglement_entropies::bits_err": ["1e-6"],
    "storage::dataset::subsystem_entanglement_entropies::eig_size": ["8192"],
    "storage::dataset::subsystem_entanglement_entropies::bond_lim": ["4096"],
    "storage::dataset::subsystem_entanglement_entropies::trnc_lim": ["1e-6"],
    "strategy::iter_max_warmup": ['8'],
    "strategy::dmrg_blocksize_policy": ["MAX"],
    "strategy::dmrg_min_blocksize": ["2"],
    "strategy::dmrg_max_blocksize": ["2"],
    "strategy::dmrg_env_expand_mode": ["H1|H2|NSITE|FORWARD"],
    "strategy::initial_state": ["PRODUCT_STATE_NEEL"],
    "strategy::bond_increase_when": ["SATURATED"],
    "strategy::bond_increase_rate": ["2.0"],
    "strategy::trnc_decrease_when": ["STUCK"],
    "strategy::trnc_decrease_rate": ["0.25"],
    "precision::eigs_iter_min": ["1000"],
    "precision::eigs_iter_max": ["50000"],
    "precision::eigs_iter_gain": ["5.0"],
    "precision::eigs_iter_gain_policy": ["SAT_VAR"],
    "precision::svd_truncation_lim": ['1e-8'],
    "precision::svd_truncation_init": ['1e-8'],
    "precision::svd_switchsize_bdc": ['16'],
    "precision::variance_convergence_threshold": ['1e-13'],
    "precision::use_energy_shifted_mpo": ["false"],
    "precision::use_parity_shifted_mpo": ["false"],
    "precision::use_parity_shifted_mpo_squared": ["true"],
    "model::model_type": ['ising_majorana'],
    "model::model_size": ['20'],
    "model::ising_majorana::g": ['0.100'],
    "model::ising_majorana::delta": ['-4.00', '-3.00',
                                     # '-2.00', '+2.00',
                                     '+3.00', '+4.00'],
    "xdmrg::energy_spectrum_shift": ['0.0'],
    "xdmrg::iter_min": ['1'],
    "xdmrg::iter_max": ['500'],
    "xdmrg::bond_max": ['8192'],
    "xdmrg::bond_min": ['48'],
}

configs = get_config_product(config_ranges, config_paths)
for config in configs:
    # Write the config file
    config['filename'] = get_config_filename(config, config_ranges, config_paths)
    config_template = config_paths['config_template']
    write_config_file(config, config_template, config['filename'])

batch_setup = get_xdmrg_batch_setup('xdmrg4-exp-3site-forw')
write_batch_files(batch_setup=batch_setup, configs=configs, config_paths=config_paths)
update_batch_status(config_paths=config_paths)
move_directories(batch_setup=batch_setup, config_paths=config_paths)