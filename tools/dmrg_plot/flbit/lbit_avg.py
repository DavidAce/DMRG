from glob import glob
from os import remove
from os import path
from dmrg_plot.common.io.h5ops import *
from dmrg_plot.common.io.parse import parse
from statistics.write_statistics import *


def lbit_avg(args):
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
            'dsets': {  # For time independent data (or at the last time step)
                'schmidt_midchain': {},
                'decay': {'axis': 0, },
                # 'number_probabilities': {
                #     'copy': False,  # Copy the dataset as is, without averaging
                #     'hartley': True # Compute the hartley number entropy
                # },
            },

            'tables': {  # For data at the last time step
                'status': ['iter',
                           'chi_lim', 'bond_lim', 'bond_limit', 'bond_dimension_limit',
                           'bond_max', 'bond_dimension_max', 'bond_dimension_maximum',
                           'phys_time', 'algo_time', 'delta_t'],
                'mem_usage': 'ALL',
                'bond_dimensions': 'ALL',
                'truncation_errors': 'ALL',
            },

            'cronos':  # For data at each time step
                {
                    'measurements': ['iter',
                                     'entanglement_entropy_midchain',
                                     'number_entropy_midchain',
                                     'hartley_number_entropy_midchain',
                                     'bond_mid', 'bond_dimension_midchain',
                                     'truncation_error',
                                     'algorithm_time',
                                     'physical_time'],
                    'bond_dims': 'ALL',
                    'bond_dimensions': 'ALL',
                    'entanglement_entropies': 'ALL',
                    'number_entropies': 'ALL',
                    'truncation_errors': 'ALL',
                    'number_probabilities': 'ALL',
                    '__save_data__': ['entanglement_entropies', 'number_entropies']
                },

        }

        write_statistics(src, tgt, data_props)


if __name__ == '__main__':
    args = parse()
    args.clear = True
    lbit_avg(args)
