from .tools import *
from src.io.h5ops import *
import matplotlib.pyplot as plt
from .filter import *


def multiplot_iter_distribution(src, plotdir='', algo_filter='', state_filter=''):
    print('Plotting: Iteration distribution for: ', algo_filter, state_filter)
    h5_src = h5open(src, 'r')
    path_L = h5py_unique_finder(h5_src, filter='L_', dep=1)
    path_l = h5py_unique_finder(h5_src, filter='l_', dep=2)
    path_d = h5py_unique_finder(h5_src, filter='d_', dep=3)
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
                for algokey, algopath, algonode in h5py_group_iterator(g=basenode, filter=algo_filter, dep=1):
                    for statekey, statepath, statenode in h5py_group_iterator(g=algonode, filter=state_filter, dep=1):
                        # chi_result = h5py_node_finder(g=statenode,filter="chi_lim",num=1,dep=8)
                        # chi_max  = chi_result[0][2]["avg"][()]
                        for win_idx, win in enumerate(variance_window_limits):
                            idx = get_v_filtered_index_list(statenode, win)
                            for datakey, datapath, datanode in h5py_node_finder(g=statenode,
                                                                                filter='iter',
                                                                                dep=8):
                                # mid_idx = int(chain_length/2)
                                if not idx:
                                    data = datanode['data']
                                else:
                                    data = datanode['data'][idx]
                                num = datanode['num'][()]
                                avg = datanode['avg'][()]
                                nicename = re.sub(r'[\W_]', ' ', str(algokey + " " + statekey))
                                nicename = nicename + ' (' + str(num) + ')'
                                color = next(current_palette)
                                ax.hist(x=np.array(data), color=color, bins=60, linewidth=1,
                                        histtype='step', density=False, label=nicename)
                                ax.axvline(avg, color=color, linestyle='dashed', linewidth=1, label='Avg ')
                ax.set_xlabel('Iteration')
                ax.set_ylabel('Histogram')
                ax.set_title('$L = ' + str(chain_length) + '$')
                if (chi_max > 0):
                    ax.set_xlim(left=0, right=chi_max)

                # h5keys = [item for item in h5_src[L][l][d].keys() if any(s in item for s in key_list)]
                # key_sorted = sorted(h5keys, key=natural_keys)
                # for key in key_sorted:
                #     try:
                #         node = h5_src[L][l][d][key]['bond_dimensions']
                #     except:
                #         continue
                #     key_num = key_num+1
                #     length = node.attrs['chain_length']
                #     delt   = node.attrs['delta']
                #     lamb   = node.attrs['lambda']
                #     chi_max= node.attrs['chi_max']
                #     middle = int(length/2)
                #     for win_idx, win in enumerate(variance_window_limits):
                #         idx = get_v_filtered_index_list(h5_src[L][l][d][key], win)
                #         data = np.array(node['data'][idx, middle])
                #         if np.any(np.isnan(data)):
                #             raise ValueError("Data contains nan's")
                #         hist, edges = np.histogram(data, bins=60, density=False)
                #         num_elems = len(data)
                #         avg = np.mean(data, axis=0)
                #         nicename = re.sub(r'[\W_]', ' ', str(key)) + ' (' + str(num_elems) +') '
                #         # nicename = nicename + ' (' + str(num_elems) +') ' + variance_window_names[win_idx][0]
                #         color = next(current_palette)
                #         ax.hist(x=np.array(data), color=color, bins=np.array(edges), linewidth=1, histtype='step',density=False, label=nicename)
                #         ax.axvline(avg, color=color, linestyle='dashed', linewidth=1, label='Avg '+str(key))
                #     ax.set_xlabel('$\chi$')
                #     ax.set_ylabel('Histogram')
                #     ax.set_title('$L = ' + str(length) + '$')
                #     if(chi_max > 0):
                #         ax.set_xlim(left=0,right=chi_max)
                used_ax = used_ax + 1
                ax.legend()
            fig.suptitle('Distribution of Iterations  @ $\Delta = $' + str(delt) + '$\lambda = $' + str(lamb))
            for ax in np.ravel(axes)[used_ax:]:
                fig.delaxes(ax)
            if plotdir != '':
                plt.savefig(plotdir + '/Iteration_distribution_' + l + '_' + d + '.pdf', format='pdf')
                plt.savefig(plotdir + '/Iteration_distribution_' + l + '_' + d + '.png', format='png')
    h5close(h5_src)
