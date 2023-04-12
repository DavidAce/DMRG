from .tools import *
from src.io.h5ops import *
import matplotlib.pyplot as plt
from .filter import *


def multiplot_mem_distribution(src, plotdir='', algo_filter='', state_filter=''):
    print('Plotting: Peak memory distribution for: ', algo_filter, state_filter)
    h5_src = h5open(src, 'r')
    path_L = h5py_unique_finder(h5_src, filter='L_', dep=1)
    path_l = h5py_unique_finder(h5_src, filter='l_', dep=2)
    path_d = h5py_unique_finder(h5_src, filter='d_', dep=3)
    # One figure per unique l, unique d
    for l in path_l:
        for d in path_d:
            # In each figure we want one subplot per unique_L
            rows, cols = get_optimal_subplot_num(len(path_L))
            fig, axes = plt.subplots(nrows=rows, ncols=cols, figsize=(7 * cols, 7 * rows))
            fig.tight_layout(pad=5, w_pad=1.0, h_pad=1.0)
            fig.subplots_adjust(wspace=0.3, hspace=0.3)
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
                        for win_idx, win in enumerate(variance_window_limits):
                            idx = get_v_filtered_index_list(statenode, win)
                            for datakey, datapath, datanode in h5py_node_finder(g=statenode,
                                                                                filter='hwm',
                                                                                dep=8):
                                if not idx:
                                    data = datanode['data']
                                else:
                                    data = datanode['data'][idx]
                                num = datanode['num'][()]
                                avg = datanode['avg'][()]
                                nicename = re.sub(r'[\W_]', ' ', str(algokey + " " + statekey))
                                nicename = nicename + ' (' + str(num) + ')'
                                color = next(current_palette)
                                ax.hist(x=np.array(data), color=color, bins=30, linewidth=1,
                                        histtype='step', density=True, label=nicename)
                                ax.axvline(avg, color=color, linestyle='dashed', linewidth=1, label='Avg ')
                ax.set_xlabel('RAM [MiB]')
                ax.set_ylabel('Histogram')
                ax.set_title('$L = ' + str(chain_length) + '$')

                used_ax = used_ax + 1
                ax.legend()
            fig.suptitle('Distribution of Peak Memory (HWM) @ $\Delta = $' + str(delt) + '$\lambda = $' + str(lamb))
            for ax in np.ravel(axes)[used_ax:]:
                fig.delaxes(ax)
            if plotdir != '':
                plt.savefig(plotdir + '/Mem_distribution_' + l + '_' + d + '.pdf', format='pdf')
                plt.savefig(plotdir + '/Mem_distribution_' + l + '_' + d + '.png', format='png')
    h5close(h5_src)
