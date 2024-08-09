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
    "storage::output_filepath": [get_output_filepath],
    "storage::resume_policy": ['IF_UNSUCCESSFUL'],
    "storage::file_collision_policy": ['REVIVE'],
    "storage::temp_dir": [config_paths['temp_dir']],
    "storage::mps::state_emid::policy": ["ITER|FINISH|REPLACE"],
    "storage::table::opdm::policy": ["FINISH"],
    "storage::table::opdm_spectrum::policy": ["FINISH"],
    "storage::dataset::subsystem_entanglement_entropies::bits_err": ["1e-6"],
    "storage::dataset::subsystem_entanglement_entropies::eig_size": ["8192"],
    "storage::dataset::subsystem_entanglement_entropies::bond_lim": ["4096"],
    "storage::dataset::subsystem_entanglement_entropies::trnc_lim": ["1e-6"],
    "console::loglevel": ['2'],
    "precision::eigs_iter_min" :["5000"],
    "precision::eigs_iter_max" :["20000"],
    "precision::eigs_iter_gain" :["1.25"],
    "precision::eigs_iter_gain_policy" : ["SATURATED"],
    "precision::svd_truncation_min": ['1e-8'],
    "precision::svd_truncation_max": ['1e-3'],
    "precision::svd_switchsize_bdc": ['16'],
    "precision::variance_convergence_threshold": ['1e-12'],
    "strategy::iter_max_warmup": ['0'],
    "strategy::initial_state": ["PRODUCT_STATE_NEEL"],
    "strategy::dmrg_blocksize_policy": ["ICOM|SAT_VAR|FIN_BON|FIN_TRNC"],
    "strategy::dmrg_min_blocksize": ["1"],
    "strategy::dmrg_max_blocksize": ["16"],
    "strategy::bond_increase_when": ["NEVER"],
    "strategy::bond_increase_rate": ["8.0"],
    "strategy::trnc_decrease_when": ["STUCK"],
    "strategy::trnc_decrease_rate": ["0.5"],
    "precision::use_energy_shifted_mpo": ["false"],
    "precision::use_parity_shifted_mpo": ["true"],
    "precision::use_parity_shifted_mpo_squared": ["true"],
    "model::model_type": ['ising_majorana'],
    "model::model_size": ['20', '22', '24'],
    "model::ising_majorana::g": ['0.000', '0.005', '0.010','0.015', '0.020','0.025', '0.030', '0.040','0.060', '0.080', '0.100' ],
    # "model::ising_majorana::delta": ['-6.00', '-5.00', '-4.00', '-3.00', '-2.00', '-1.00', '-0.75', '-0.50' , '-0.25', '+0.00', '+0.25',  '+0.50' ,  '+0.75', '+1.00', '+2.00', '+3.00', '+4.00', '+5.00', '+6.00'],
    "model::ising_majorana::delta": ['-6.00', '-5.00', '-4.00', '-3.00', '-2.00', '+2.00', '+3.00', '+4.00', '+5.00', '+6.00'],
    "xdmrg::energy_spectrum_shift": ['1e-6'],
    "xdmrg::iter_min": ['1'],
    "xdmrg::iter_max": ['1000'],
    "xdmrg::bond_max": ['8192'],
    "xdmrg::bond_min": ['16'],
    "xdmrg::print_freq": ['1'],
    "xdmrg::try_directx2_when_stuck": ['false'],
    "xdmrg::try_shifting_when_degen": ['0'],
}

configs = get_config_product(config_ranges, config_paths)
for config in configs:
    # Write the config file
    config['filename'] = get_config_filename(config, config_ranges, config_paths)
    config_template = config_paths['config_template']
    write_config_file(config, config_template, config['filename'])

batch_setup = get_xdmrg_batch_setup('xdmrg3-letsgo')
write_batch_files(batch_setup=batch_setup, configs=configs, config_paths=config_paths)
update_batch_status(config_paths=config_paths)
move_directories(batch_setup=batch_setup, config_paths=config_paths)