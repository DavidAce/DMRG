from .tools import *
from src.io.h5ops import *
import matplotlib.pyplot as plt
from .filter import *


def multiplot_chi_distribution(src, plotdir='', algo_inc='', state_inc=''):
    print('Plotting: Bond dimension distribution for: ', algo_inc, state_inc)
    h5_src = h5open(src, 'r')
    path_L = h5py_unique_finder(h5_src, keypattern='L_', dep=1)
    path_l = h5py_unique_finder(h5_src, keypattern='l_', dep=2)
    path_d = h5py_unique_finder(h5_src, keypattern='d_', dep=3)
    # One figure per unique l, unique d
    for l in path_l:
        for d in path_d:
            # In each figure we want one subplot per unique_L
            rows, cols = get_optimal_subplot_num(len(path_L))
            fig, axes = plt.subplots(nrows=rows, ncols=cols, figsize=(4.5 * cols, 4.5 * rows))
            fig.tight_layout(pad=5, w_pad=1.0, h_pad=1.0)
            fig.subplots_adjust(wspace=0.2, hspace=0.2)
            used_ax = 0
            delt = 0
            lamb = 0
            chi_max = 0
            for ax, L in zip(np.ravel(axes), path_L):
                current_palette = itertools.cycle(sns.color_palette())
                key_num = 0
                basenode = h5_src[L][l][d]
                chain_length = basenode.attrs['model_size']
                delt = basenode.attrs['delta']
                lamb = basenode.attrs['lambda']
                midx = int(chain_length / 2)
                mlbl = 'L_{}'.format(midx)
                for algokey, algopath, algonode in h5py_group_iterator(node=basenode, keypattern=algo_inc, dep=1):
                    for statekey, statepath, statenode in h5py_group_iterator(node=algonode, keypattern=state_inc, dep=1):
                        # chi_result = h5py_node_finder(node=statenode,keypattern=['chi_lim','status'],num=1,dep=8)
                        # chi_max  = chi_result[0][2]["avg"][()]
                        for win_idx, win in enumerate(variance_window_limits):
                            idx = get_v_filtered_index_list(statenode, win)
                            for datakey, datapath, datanode in h5py_node_finder(node=statenode,
                                                                                keypattern=['bond_dimensions'],
                                                                                dep=8):
                                if "checkpoint" in datapath:
                                    continue
                                if 'states' in statekey:
                                    data = datanode['data'][midx, :]
                                    avg = datanode['avg'][midx]
                                    num = datanode['num'][()]
                                else:
                                    data = get_data(datanode['data'], keys=mlbl, dtype='f8')[idx]
                                    avg = np.mean(data)
                                    num = np.size(data)
                                    chi_max = np.max([chi_max, np.max(data)])
                                nicename = re.sub(r'[\W_]', ' ', str(algokey + " " + statekey))
                                nicename = nicename + ' (' + str(num) + ')'
                                color = next(current_palette)
                                ax.hist(x=np.array(data), color=color, bins=60, linewidth=1,
                                        histtype='step', density=True, label=nicename)
                                ax.axvline(avg, color=color, linestyle='dashed', linewidth=1, label='Avg ')
                ax.set_xlabel('$\chi$')
                ax.set_ylabel('Histogram')
                ax.set_title('$L = ' + str(chain_length) + '$')
                ax.set_xlim(left=0, right=chi_max)

                used_ax = used_ax + 1
                ax.legend()
            fig.suptitle('Distribution of mid-chain bond dimension @ $\Delta = $' + str(delt) + '$\lambda = $' + str(lamb))
            for ax in np.ravel(axes)[used_ax:]:
                fig.delaxes(ax)
            if plotdir != '':
                plt.savefig(plotdir + '/Chi_distribution_' + l + '_' + d + '.pdf', format='pdf')
                plt.savefig(plotdir + '/Chi_distribution_' + l + '_' + d + '.png', format='png')
    h5close(h5_src)
