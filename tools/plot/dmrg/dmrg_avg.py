from glob import glob
from src.h5ops import *
from batches import get_batches
from src.write_statistics import *

def dmrg_avg(args):
    projdir = args.basedir
    batchdirs = []
    for batch in args.batches:
        batchdirs.extend(glob('{}/{}*'.format(projdir, batch)))

    for batchdir in batchdirs:
        print(f'Averaging: {batchdir}')
        analysisdir = batchdir + '/analysis'
        plotdir = analysisdir + '/plots'
        datadir = analysisdir + '/data'

        src = datadir + '/merged.h5'
        tgt = datadir + '/averaged.h5'
        if os.path.isfile(tgt):
            if args.clear:
                print("Removing file: {}".format(tgt))
                os.remove(tgt)
            else:
                continue
        h5_tgt = h5open(tgt, 'a')
        h5close(h5_tgt)

        data_props = {
            'dsets': {  # For time independent data (or at the last time step)
                'model/hamiltonian': {'copy': True, },
                'opdm' : {'copy': True, },
            },

            'tables': {  # For data at the last time step
                'status': ['iter', 'bond_lim', 'bond_max', 'trnc_lim','algo_time'],
                # 'mem_usage': 'ALL',
                'measurements': 'ALL',
                'memory': 'ALL',
                'bond_dimensions': 'ALL',
                'opdm_spectrum': 'ALL',
                'truncation_errors': 'ALL',
                # If the following tables exist, save it under <tablename>/data
                '__save_data__': ['mem_usage', 'measurements', 'opdm_spectrum'],
            },


        }

        write_statistics(src, tgt, data_props)


if __name__ == '__main__':
    batch = get_batches('xDMRG', ['xdmrg3-letsgo'], states=['state_emid'], basedir='/mnt/WDB-AN1500/mbl_transition')
    batch.clear = True
    dmrg_avg(batch)
