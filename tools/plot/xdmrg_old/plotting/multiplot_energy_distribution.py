from .tools import *
from src.io.h5ops import *
import matplotlib.pyplot as plt
from .filter import *


def multiplot_energy_distribution(src, plotdir='', algo_inc='', state_inc=''):
    print('Plotting: energy density distribution for: ', algo_inc, state_inc)
    h5_src = h5open(src, 'r')
    path_L = h5py_unique_finder(h5_src, keypattern='L_', dep=1)
    path_l = h5py_unique_finder(h5_src, keypattern='l_', dep=2)
    path_d = h5py_unique_finder(h5_src, keypattern='d_', dep=3)
    timecounter = 0
    simcounter = 0
    # One figure per unique_l, unique_J and unique_h
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
                ed_palette = itertools.cycle(sns.color_palette("Set2"))
                current_palette = itertools.cycle(sns.color_palette())
                key_num = 0
                basenode = h5_src[L][l][d]
                chain_length = basenode.attrs['model_size']
                delt = basenode.attrs['delta']
                lamb = basenode.attrs['lambda']
                for algokey, algopath, algonode in h5py_group_iterator(node=basenode, keypattern=algo_inc, dep=1):
                    for statekey, statepath, statenode in h5py_group_iterator(node=algonode, keypattern=state_inc, dep=1):
                        idx = get_v_filtered_index_list(statenode, variance_window_limits[0])
                        for datakey, datapath, datanode in h5py_node_finder(node=statenode, keypattern=['status', 'energy_dens'], dep=8):
                            if "checkpoint" in datapath:
                                continue
                            if 'states' in statekey and 'energy_dens' in datakey:
                                data = np.array(datanode['data'])
                            else:
                                data = np.array(datanode['data']['energy_dens'])[idx]

                            if np.any(np.isnan(data)):
                                raise ValueError("Data contains nan's")
                            key_num = key_num + 1
                            hist, edges = np.histogram(data, bins=100, density=True)
                            num = np.size(data)
                            avg = np.mean(data)
                            bincentres = [(edges[i] + edges[i + 1]) / 2. for i in range(len(edges) - 1)]
                            widths = np.diff(edges)
                            if "states" in statekey:
                                color = next(ed_palette)
                                # nicename = "ED e=[" + statenode.attrs["efmt"] + "]"
                                nicename = "ED"
                                lwidth = 3.5
                                lalpha = 1.0
                            else:
                                color = next(current_palette)
                                # nicename = re.sub(r'[\W_]', ' ', str(algokey + " " + statekey))
                                nicename = re.sub(r'[\W_]', ' ', str(statekey))
                                lwidth = 0.9
                                lalpha = 0.7
                            nicename = nicename + ' (' + str(num) + ')'
                            ax.step(bincentres, hist, where='mid', label=nicename,
                                    color=color, alpha=lalpha, linewidth=lwidth)
                            ax.axvline(avg, color=color, linestyle='dashed', linewidth=lwidth)
                            timecounter = timecounter + np.sum(data)
                            simcounter = simcounter + num
                ax.set_xlabel('Energy density')
                ax.set_ylabel('Histogram')
                ax.set_title('$L = ' + str(chain_length) + '$')
                ax.set_xlim(left=0, right=1)
                # ax.set_yscale('log')
                # ax.legend()
                ax.legend(framealpha=0.2, fontsize='small', labelspacing=0.25, ncol=2, loc='lower center')
                used_ax = used_ax + 1
            fig.suptitle('Energy density distribution @ $\Delta = ' + str(delt) + '\quad \lambda = ' + str(lamb) + '$')
            for ax in np.ravel(axes)[used_ax:]:
                fig.delaxes(ax)
            if plotdir != '':
                plt.savefig(plotdir + '/energy_dens_distribution_' + l + '_' + d + '.pdf', format='pdf')
                plt.savefig(plotdir + '/energy_dens_distribution_' + l + '_' + d + '.png', format='png')
    h5close(h5_src)
