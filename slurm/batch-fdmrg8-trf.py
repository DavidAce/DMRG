from utils.generators import get_config_product, write_config_file, write_batch_files, move_directories
from utils.fdmrg import get_output_filepath, get_config_filename, update_batch_status
from batches_fdmrg import get_fdmrg_batch_setup
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
    "model::model_type": ['ising_majorana'],
    "model::model_size": ['32','48','64'],
    "model::ising_majorana::g": ['0.500'],
    "model::ising_majorana::delta": ['-9.00','-6.00','-3.00',  '+3.00', '+4.00', '+5.00', '+6.00','+9.00'],
    "storage::output_filepath": [get_output_filepath],
    "storage::resume_policy": ['IF_UNSUCCESSFUL'],
    "storage::file_collision_policy": ['REVIVE'],
    "storage::temp_dir": [config_paths['temp_dir']],
    "storage::mps::state_emin::policy": ["NONE"],
    "storage::table::opdm::policy": ["FINISH"],
    "storage::table::opdm_spectrum::policy": ["FINISH"],
    "storage::dataset::subsystem_entanglement_entropies::bits_err": ["1e-8"],
    "storage::dataset::subsystem_entanglement_entropies::eig_size": ["8192"],
    "storage::dataset::subsystem_entanglement_entropies::bond_lim": ["4096"],
    "storage::dataset::subsystem_entanglement_entropies::trnc_lim": ["1e-8"],
    "precision::svd_truncation_lim": ['1e-9'],
    "precision::svd_truncation_init": ['1e-8'],
    "precision::variance_convergence_threshold": ['1e-12'],
    "strategy::initial_state": ["PRODUCT_STATE_NEEL"],
    "strategy::dmrg_blocksize_policy": ["MAXSTUCK"],
    "strategy::dmrg_min_blocksize": ["1"],
    "strategy::dmrg_max_blocksize": ["2"],
    "strategy::dmrg_max_prob_size": ["268435456"],
    "precision::use_parity_shifted_mpo": ["true"],
    "precision::use_parity_shifted_mpo_squared": ["true"],
    "fdmrg::ritz": ['SR'],
    "fdmrg::iter_min": ['6'],
    "fdmrg::iter_max": ['30'],
    "fdmrg::bond_max": ['1024'],
    "fdmrg::bond_min": ['16'],
    "fdmrg::print_freq": ['1'],
    "fdmrg::store_wavefn": ['false'],
}

configs = get_config_product(config_ranges, config_paths)
for config in configs:
    # Write the config file
    config['filename'] = get_config_filename(config, config_ranges, config_paths)
    config_template = config_paths['config_template']
    write_config_file(config, config_template, config['filename'])

batch_setup = get_fdmrg_batch_setup('fdmrg8-trf')
write_batch_files(batch_setup=batch_setup, configs=configs, config_paths=config_paths)
update_batch_status(config_paths=config_paths)
move_directories(batch_setup=batch_setup, config_paths=config_paths)