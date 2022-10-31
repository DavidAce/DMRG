from src.plotting.tools import *
from src.io.h5ops import *
import numpy
import matplotlib.pyplot as plt
from src.plotting.filter import *
from tqdm import tqdm
import warnings


def multiplot_schmidt(src, plotdir='', algo_filter='', state_filter=''):
    print('Plotting: Schmidt values')
    h5_src = h5open(src, 'r')
    path_L = h5py_unique_finder(h5_src, filter='L_', dep=1)
    path_l = h5py_unique_finder(h5_src, filter='l_', dep=2)
    path_d = h5py_unique_finder(h5_src, filter='d_', dep=3)
    # One figure per unique l and unique d
    for l in path_l:
        for d in path_d:
            rows, cols = get_optimal_subplot_num(len(path_L))
            fig, axes = plt.subplots(nrows=rows, ncols=cols, figsize=(4.5 * cols, 4.5 * rows), sharey='all')
            fig.tight_layout(pad=5, w_pad=1.0, h_pad=1.0)
            fig.subplots_adjust(wspace=0.2, hspace=0.2)
            # chain_length = h5_src[L][path_l[0]][path_d[0]]['schmidt'].attrs['chain_length']
            used_ax = 0
            delt = 0
            lamb = 0
            for ax, L in zip(np.ravel(axes), path_L):
                current_palette = itertools.cycle(sns.color_palette())
                basenode = h5_src[L][l][d]
                chain_length = basenode.attrs['model_size']
                delt = basenode.attrs['delta']
                lamb = basenode.attrs['lambda']
                for algokey, algopath, algonode in h5py_group_iterator(g=basenode, filter=algo_filter, dep=1):
                    for statekey, statepath, statenode in h5py_group_iterator(g=algonode, filter=state_filter, dep=1):
                        for datakey, datapath, datanode in h5py_node_finder(g=statenode, filter='schmidt_midchain',
                                                                            dep=8):
                            data = np.array(datanode['data'])
                            num = datanode['num'][()]
                            nonzeros = np.count_nonzero(data, axis=1)
                            nonzeros[nonzeros == 0] = 1
                            data[data == 0] = np.nan
                            # I expect to see RuntimeWarnings in this block
                            with warnings.catch_warnings():
                                warnings.simplefilter("ignore", category=RuntimeWarning)
                                ydata = np.nanmean(data, axis=1)
                                edata = np.nanstd(data, axis=1) / np.sqrt(nonzeros)
                            xdata = range(len(ydata))
                            color = next(current_palette)
                            nicename = re.sub(r'[\W_]', ' ', str(statekey))
                            nicename = nicename + ' (' + str(num) + ')'
                            ax.fill_between(x=xdata, y1=ydata - edata, y2=ydata + edata, color=color, alpha=0.2,
                                            label=None)
                            ax.plot(xdata, ydata, color=color, marker='.', markersize=1,
                                    label=nicename)

                ax.set_yscale("log", nonpositive='clip')
                ax.set_ylabel('$\langle \lambda_n \\rangle$ ')
                ax.set_xlabel('$n$')
                ax.set_title('$L = ' + str(chain_length) + '$')
                ax.legend()
                used_ax = used_ax + 1
            fig.suptitle('Average Schmidt values @ $\Delta = ' + str(delt) + '\quad \lambda = ' + str(lamb) + '$')
            for ax in np.ravel(axes)[len(path_L):]:
                fig.delaxes(ax)
            if plotdir != '':
                plt.savefig(plotdir + '/schmidt_' + l + '_' + d + '.pdf', format='pdf')
                plt.savefig(plotdir + '/schmidt_' + l + '_' + d + '.png', format='png')

    h5close(h5_src)
