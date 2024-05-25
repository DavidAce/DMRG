import argparse


def get_batches(algo, batches=None, basedir='/mnt/WDB-AN1500/mbl_transition', states=None):
    if batches is not None:
        if not isinstance(batches, list):
            raise TypeError("batches must be a list")
    fdmrg_batches = ['fdmrg-see-test4',  # Subsystem entanglement entropy test for L=25
    ]
    xdmrg_batches = ['xdmrg1-fse',  # First test after spending years on the lbit project
                     'xdmrg2-tuned',  # Fixed various issues that arised in xdmrg1-fse.
                     'xdmrg3-letsgo',  # Fixed various issues that arised in xdmrg1-fse.
                     ]



    parser = argparse.ArgumentParser(description='dmrg-plot')
    parser.add_argument('--clear', action='store_true', help='Remake averaged.h5')
    parser.add_argument('--basedir', type=str, help='Main directory for all mbl batch data', default=basedir)
    parser.add_argument('--algos', type=list, help='List of algorithms to plot data for', choices=['xDMRG','fDMRG'], default=[algo])
    parser.add_argument('--states', type=list, help='List of states to plot data for', default=states if states is not None else ['state_emid', 'state_emin', 'state_emax'])
    parser.add_argument('--models', type=list, help='List of models to plot data for', default=['analysis', 'model'])
    parser.add_argument('--points', type=list, help='List of points to plot data for', default=['tables'])
    if algo == 'fDMRG':
        parser.add_argument('--batches', type=list, help='List of batches to plot data for',
                            default=[fdmrg_batches[-1]] if batches is None else batches)
    elif algo == 'xDMRG':
        parser.add_argument('--batches', type=list, help='List of batches to plot data for',
                            default=[xdmrg_batches[-1]] if batches is None else batches)

    args = parser.parse_args()
    return args
