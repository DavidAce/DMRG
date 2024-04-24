from .tools import *
from src.io.h5ops import *
import matplotlib.pyplot as plt
from .filter import *


def multiplot_S_distribution(src, plotdir='', algo_inc='', state_inc='', bins=200):
    print('Plotting: Entanglement entropy distribution for: ', algo_inc, state_inc)
    h5_src = h5open(src, 'r')
    path_L = h5py_unique_finder(h5_src, keypattern='L_', dep=1)
    path_l = h5py_unique_finder(h5_src, keypattern='l_', dep=2)
    path_d = h5py_unique_finder(h5_src, keypattern='d_', dep=4)

    # One figure per unique_l, unique_J and unique_h
    for l in path_l:
        for d in path_d:
            # In each figure we want one subplot per unique_L
            rows, cols = get_optimal_subplot_num(len(path_L))
            fig, axes = plt.subplots(nrows=rows, ncols=cols, figsize=(4.5 * cols, 4.5 * rows), sharex='all', sharey='all')
            fig.tight_layout(pad=5, w_pad=1.0, h_pad=1.0)
            fig.subplots_adjust(wspace=0.2, hspace=0.2)
            used_ax = 0
            delt = 0
            lamb = 0
            for ax, L in zip(np.ravel(axes), path_L):
                if not L in h5_src or not l in h5_src[L] or not d in h5_src[L][l]:
                    continue
                basenode = h5_src[L][l][d]
                ed_palette = itertools.cycle(sns.color_palette("Set2"))
                current_palette = itertools.cycle(sns.color_palette("colorblind", 5))
                chain_length = basenode.attrs['model_size']
                delt = basenode.attrs['delta']
                lamb = basenode.attrs['lambda']
                midx = int(chain_length / 2)
                mlbl = 'L_{}'.format(midx)
                max_S = 0
                for algokey, algopath, algonode in h5py_group_iterator(node=basenode, keypattern=algo_inc, dep=1):
                    for statekey, statepath, statenode in h5py_group_iterator(node=algonode, keypattern=state_inc, dep=1):
                        for win_idx, win in enumerate(variance_window_limits):
                            idx = get_v_filtered_index_list(statenode, win)
                            for datakey, datapath, datanode in h5py_node_finder(node=statenode, keypattern='entanglement_entropies', dep=3):
                                if "checkpoint" in datapath:
                                    continue
                                if not idx or 'states' in statekey:
                                    data = datanode['data'][midx, :]
                                    num = datanode['num'][()]
                                    avg = datanode['avg'][midx]
                                else:
                                    data = get_data(datanode['data'], keys=mlbl, dtype='f8')[idx]
                                    avg = np.mean(data)
                                    num = np.size(data)
                                if np.any(np.isnan(data)):
                                    raise ValueError("Data contains nan's")
                                datarange = [np.min(data), np.max(data)]
                                hist, edges = np.histogram(data, bins=bins, range=datarange, density=False)
                                bincentres = [(edges[i] + edges[i + 1]) / 2. for i in range(len(edges) - 1)]
                                widths = np.diff(edges)
                                norm = np.dot(hist, widths)
                                if "states" in statekey:
                                    color = next(ed_palette)
                                    nicename = "ED e=[" + statenode.attrs["efmt"] + "]"
                                    lwidth = 2.4
                                    lalpha = 0.8
                                else:
                                    color = next(current_palette)
                                    nicename = re.sub(r'[\W_]', ' ', str(algokey + " " + statekey))
                                    lwidth = 1.4
                                    lalpha = 0.9
                                    max_S = np.max([max_S, np.max(data)])
                                nicename = nicename + ' (' + str(num) + ')'
                                ax.step(bincentres, hist / norm, where='mid', label=nicename, linewidth=lwidth, alpha=lalpha, color=color)
                                ax.axvline(avg, linestyle='dashed', linewidth=lwidth, alpha=lalpha, color=color)

                    if max_S > 0:
                        ax.set_xlim(0, max_S)
                    ax.set_xlim(0, 2)
                    ax.set_xlabel('$S_E$')
                    ax.set_ylabel('$P(S_E)$')
                    ax.set_title('$L = ' + str(chain_length) + '$')
                    # ax.set_xlim(1e-21,100)
                used_ax = used_ax + 1
                ax.legend()
            fig.suptitle(
                'Distribution of mid-chain entanglement entropy @ $\Delta = $' + str(delt) + '$\lambda = $' + str(lamb))
            for ax in np.ravel(axes)[used_ax:]:
                fig.delaxes(ax)
            if plotdir != '':
                Jh = re.sub('/', '_', str(d))
                plt.savefig(plotdir + '/S_distribution_' + l + '_' + d + '.pdf', format='pdf')
                plt.savefig(plotdir + '/S_distribution_' + l + '_' + d + '.png', format='png')
    h5close(h5_src)
