import argparse
import flbit.lbit_plot

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='dmrg-plot')
    parser.add_argument('--algo', type=list, help='List of algorithms to plot data for', choices=['fLBIT', 'xDMRG'], required=True)
    args, unknownargs = parser.parse_known_args()

    # if args.algo == 'fLBIT':
    #     lbit_plot
