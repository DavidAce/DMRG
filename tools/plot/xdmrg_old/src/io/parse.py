import argparse


def parse():
    parser = argparse.ArgumentParser(description='lbit-plot')
    parser.add_argument('--clear', action='store_true', help='Remake averaged.h5')
    parser.add_argument('--basedir', type=str, help='Main directory for all mbl batch data', default='/mnt/WDB-AN1500/mbl_transition')
    parser.add_argument('--algos', type=list, help='List of algorithms to plot data for', default=['fLBIT'])
    parser.add_argument('--states', type=list, help='List of states to plot data for', default=['state_real'])
    parser.add_argument('--models', type=list, help='List of models to plot data for', default=['analysis', 'model'])
    parser.add_argument('--points', type=list, help='List of points to plot data for', default=['tables'])
    parser.add_argument('--batches', type=list, help='List of batches to plot data for',
                        default=[
                            # 'lbit25',
                            # 'lbit26',
                            # 'lbit27',
                            # 'lbit28',
                            # 'lbit29', # swap-less: Tried to reduce the number of swap gates; L:24, w2:0.50,0.25, x:0.5, 1.0, f:0.15,0.25, u:2,4,6
                            # 'lbit30', # time-less: Skip time-evo when gates are ~1. L:24, w2:0.50,0.25, x:0.5, 1.0, f:0.15,0.25, u:2,4,6
                            # 'lbit31', # Try higher f to promote lbit movement. L:24, w2:0.25, x:0.5, f:0.30, u:5
                            # 'lbit32', # Check new numb.entr, time-evo and xi. L:24, w2:0.25, x:0.5, f:0.25, u:5
                            # 'lbit33', # Push to L32 with new numb.entr, time-evo and xi.  L:32, w2:0.25, x:0.5, f:0.25, u:5
                            # 'lbit34', # Compare u4 and u4: L:24, w2:0.25, x:0.5, f:0.25, u:4,5
                            # 'lbit35', # unwave: Check if more precision can remove waviness (NO):  L:24, w2:0.25, x:0.5, f:0.25, u:4
                            # 'lbit36', # dearest-neighbor: Update exp decay to exp(-(r-1)/xi. L:24, w2:0.25, x:0.5, f:0.25, u:4
                            # 'lbit37', # unwave-bump-x: Check if f ~ 0.15 can remove waviness (NO). L:24, w2:0.25, x:0.5, f:0.15, u:4
                            # 'lbit38', # unwave-bump-x: Check if x ~ 1.00 can remove waviness (YES). L:24, w2:0.25, x:1.0, f:0.25, u:4
                            # 'lbit39', # serious-trial-L28: Check loglog window. Predict ~5 decades of S_N growth (YES!):  L:12-28, w2:0.25, x:0.8, f:0.25, u:4
                            # 'lbit40', # Mini test to check progress on new storage form
                            # 'lbit41', # Mini test to check progress on new storage form
                            # 'lbit42', # Big run, 5000 realizations L = 8..24
                            # 'lbit43', # Small test w-dependence, L = 16
                            # 'lbit44',  # Small test of normally distributed random circuit L = [4...16]
                            # 'lbit45',  # Small test of squared distributed random circuit L = [4...16]
                            # 'lbit46',  # Small test of squared distributed random circuit L = [8...16]
                            # 'lbit47',  # Small test of big bias on the J2 and J3 parameters L = [8...16]
                            # 'lbit48',  # Small test of choked unitaries on the J2 and J3 parameters L = [8...16]
                            'lbit49',  # Small test of choked unitaries on the J2 and J3 parameters L = [8...16]
                        ])

    args = parser.parse_args()
    return args
