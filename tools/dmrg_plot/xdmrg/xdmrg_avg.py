from glob import glob
from os import remove
from os import path
from statistics.write_statistics import *
from dmrg_plot.common.io.parse import parse
from dmrg_plot.common.io.h5ops import *

np.set_printoptions(linewidth=320)  # Adjust pycharm width


def xdmrg_avg(args):
    projdir = args.basedir
    batchdirs = []
    for batch in args.batches:
        batchdirs.extend(glob('{}/{}*'.format(projdir, batch)))

    for batchdir in batchdirs:
        analysisdir = batchdir + '/analysis'
        plotdir = analysisdir + '/plots'
        datadir = analysisdir + '/data'

        src = datadir + '/merged.h5'
        tgt = datadir + '/averaged.h5'
        if path.isfile(tgt):
            if args.clear:
                print("Removing file: {}".format(tgt))
                remove(tgt)
            else:
                continue
        h5_tgt = h5open(tgt, 'a')
        h5close(h5_tgt)

        data_props = {
            'algo_list': ['fDMRG', 'xDMRG'],
            'state_list': ['state_'],
            'point_list': ['finished'],
            'tables': {'measurements': ['iter',
                                        'energy',
                                        'energy_dens',
                                        'energy_variance',
                                        'entanglement_entropy_midchain',
                                        'bond_mid', 'bond_dimension_midchain',
                                        'spin_components',
                                        'truncation_error',
                                        'algorithm_time'],
                       # 'profiling' : ['t_tot', 't_sim'],
                       'status': ['iter',
                                  'bond_lim', 'bond_limit', 'bond_dimension_limit',
                                  'bond_max', 'bond_dimension_max', 'bond_dimension_maximum',
                                  'energy_dens',
                                  'wall_time'],
                       'mem_usage': ['rss', 'hwm'],
                       # 'bond_dimensions' : 'ALL',
                       'bond_dims': 'ALL',
                       'entanglement_entropies': 'ALL',
                       'renyi_2': 'ALL',
                       'renyi_3': 'ALL',
                       'renyi_4': 'ALL',
                       'truncation_errors': 'ALL',
                       'schmidt_midchain': 'ALL',
                       },
            'dsets': [
                'correlation_matrix_sx',
                'correlation_matrix_sy',
                'correlation_matrix_sz',
            ],

            'fes': {'measurements': ['iter',
                                     'energy',
                                     'energy_dens',
                                     'energy_variance',
                                     'entanglement_entropy_midchain',
                                     'bond_mid', 'bond_dimension_midchain',
                                     'spin_components',
                                     'truncation_error',
                                     'algorithm_time'],
                    # 'profiling' : ['t_tot', 't_sim'],
                    'status': ['iter',
                               'bond_lim', 'bond_limit', 'bond_dimension_limit',
                               'bond_max', 'bond_dimension_max', 'bond_dimension_maximum',
                               'energy_dens',
                               'wall_time'],
                    # 'mem_usage': ['rss', 'hwm'],
                    # 'bond_dimensions': 'ALL',
                    'bond_dims': 'ALL',
                    'entanglement_entropies': 'ALL',
                    'renyi_2': 'ALL',
                    'renyi_3': 'ALL',
                    'renyi_4': 'ALL',
                    'truncation_errors': 'ALL',
                    'schmidt_midchain': 'ALL',
                    },
            'bondpoint': {'measurements': ['iter',
                                           'energy',
                                           'energy_dens',
                                           'energy_variance',
                                           'entanglement_entropy_midchain',
                                           'bond_mid', 'bond_dimension_midchain',
                                           'spin_components',
                                           'truncation_error',
                                           'algorithm_time'],
                          # 'profiling' : ['t_tot', 't_sim'],
                          'status': ['iter',
                                     'bond_lim', 'bond_limit', 'bond_dimension_limit',
                                     'bond_max', 'bond_dimension_max', 'bond_dimension_maximum',
                                     'energy_dens',
                                     'wall_time'],
                          # 'mem_usage': ['rss', 'hwm'],
                          # 'bond_dimensions': 'ALL',
                          'bond_dims': 'ALL',
                          'entanglement_entropies': 'ALL',
                          'renyi_2': 'ALL',
                          'renyi_3': 'ALL',
                          'renyi_4': 'ALL',
                          'truncation_errors': 'ALL',
                          'schmidt_midchain': 'ALL',
                          },
        }

        write_statistics(src, tgt, data_props)


if __name__ == '__main__':
    args = parse('xDMRG')
    args.clear = True
    xdmrg_avg(args)
