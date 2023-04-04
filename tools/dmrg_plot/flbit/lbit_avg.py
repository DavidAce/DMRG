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
                'model/hamiltonian': {'copy': True, },
                'model/model_size': {'copy': True, },
                'model/lbits/cls_avg_fit': {'axis': 0, },
                'model/lbits/cls_avg_rms': {'axis': 0, },
                'model/lbits/cls_avg_rsq': {'axis': 0, },
                'model/lbits/cls_typ_fit': {'axis': 0, },
                'model/lbits/cls_typ_rms': {'axis': 0, },
                'model/lbits/cls_typ_rsq': {'axis': 0, },
                'model/lbits/corravg': {'axis': 0, },
                'model/lbits/corrtyp': {'axis': 0, },
                'model/lbits/correrr': {'axis': 0, },
                'model/lbits/decay_avg': {'axis': 0, },
                'model/lbits/decay_err': {'axis': 0, },
                'model/lbits/corrmat': {
                    'copy': True,
                },
                'model/lbits/corroff': {
                    'copy': True,
                },
                'model/lbits/data': {
                    'copy': True,
                },
                'model/lbits/data_shifted': {
                    'copy': True,
                },
                # 'number_probabilities': {
                #     'copy': False,  # Copy the dataset as is, without averaging
                #     'hartley': True,  # Compute the hartley number entropy
                # },
            },

            'tables': {  # For data at the last time step
                'status': ['iter',
                           'bond_lim', 'bond_max',
                           'phys_time', 'algo_time', 'delta_t'],
                'mem_usage': 'ALL',
                'bond_dimensions': 'ALL',
                'truncation_errors': 'ALL',
            },

            'cronos':  # For data at each time step
                {
                    'measurements': ['iter',
                                     'entanglement_entropy',
                                     'number_entropy',
                                     'hartley_number_entropy',
                                     'bond_mid', 'bond_lim',
                                     'truncation_error',
                                     'algorithm_time',
                                     'physical_time'],
                    'bond_dims': 'ALL',
                    'bond_dimensions': 'ALL',
                    'entanglement_entropies': 'ALL',
                    'number_entropies': 'ALL',
                    'truncation_errors': 'ALL',
                    'number_probabilities': 'ALL',
                    # If a table with this name exists, save its midchain column to a new dataset with the same name
                    '__save_mid__': ['entanglement_entropies', 'number_entropies'],
                    # If a table has this column, save all columns to a new dataset with the same name
                    '__save_col__': ['entanglement_entropy', 'number_entropy', 'bond_mid', 'algorithm_time'],
                },

        }

        write_statistics(src, tgt, data_props)


if __name__ == '__main__':
    args = parse('fLBIT', ['lbit82'])
    args.clear = True
    lbit_avg(args)
