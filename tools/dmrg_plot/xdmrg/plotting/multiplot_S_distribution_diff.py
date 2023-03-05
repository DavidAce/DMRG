from src.plotting.tools import *
from src.io.h5ops import *
import matplotlib.pyplot as plt
from src.plotting.filter import *


def multiplot_S_distribution_diff(src, plotdir='', algo_filter='', state_filter='', key_comp='', normalized=True, bins=200):
    print('Plotting: Difference in Entanglement entropy distribution for: ', key_comp, ' vs ', algo_filter, state_filter)
    h5_src = h5open(src, 'r')
    path_L = h5py_unique_finder(h5_src, keypattern='L_', dep=1)
    path_l = h5py_unique_finder(h5_src, keypattern='l_', dep=2)
    path_d = h5py_unique_finder(h5_src, keypattern='d_', dep=3)

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
                current_palette = itertools.cycle(sns.color_palette())
                basenode = h5_src[L][l][d]
                chain_length = basenode.attrs['model_size']
                delt = basenode.attrs['delta']
                lamb = basenode.attrs['lambda']
                # h5keys = [item for item in h5_src[L][l][d].keys() if any(s in item for s in key_list)]

                if isinstance(key_comp, list) and len(key_comp) > 1:
                    raise ValueError("This function only takes one comparator key, got: " + str(key_comp) + "\nTry being more specific")
                cmkeys = [item for item in basenode.keys() if any(s in item for s in [key_comp])]
                if len(cmkeys) == 0:
                    print("Could not match the given pattern: " + str(key_comp) + " to any of the keys: " + str(basenode.keys()))
                    continue
                if len(cmkeys) > 1:
                    print("The pattern provided matched several keys: " + str(cmkeys) + "\nTry being more specific")
                    continue
                compresult = h5py_node_finder(node=basenode[cmkeys[0]], keypattern="entanglement_entropy", num=1, dep=8)
                if not compresult:
                    continue
                compkey, comppath, compnode = compresult[0]
                for algokey, algopath, algonode in h5py_group_iterator(node=basenode, keypattern=algo_filter, dep=1):
                    for statekey, statepath, statenode in h5py_group_iterator(node=algonode, keypattern=state_filter, dep=1):
                        for win_idx, win in enumerate(variance_window_limits):
                            idx = get_v_filtered_index_list(statenode, win)
                            compidx = get_v_filtered_index_list(basenode[cmkeys[0]], win)
                            for datakey, datapath, datanode in h5py_node_finder(node=statenode, keypattern='entanglement_entropy', dep=8):
                                if "checkpoint" in datapath:
                                    continue
                                if not idx:
                                    data = datanode['data']
                                else:
                                    data = datanode['data'][idx]
                                color = next(current_palette)
                                if np.any(np.isnan(data)):
                                    raise ValueError("Data contains nan's")
                                datarange = [np.min(data), np.max(data)]
                                hist, edges = np.histogram(data, range=datarange, bins=bins, density=False)
                                num_elems = len(data)
                                avg = np.mean(data)
                                nicename = re.sub(r'[\W_]', ' ', str(algokey + " " + statekey))
                                ax.axvline(avg, color=color, linestyle='dashed', linewidth=1)

                                # bincentres = [(edges[i] + edges[i + 1]) / 2. for i in range(len(edges) - 1)]
                                widths = np.diff(edges)
                                norm = np.dot(hist, widths)
                                hist1 = hist / norm
                                # Now get the the data from the comparator key to subtract
                                if not compidx:
                                    comp = compnode['data']
                                else:
                                    comp = compnode['data'][compidx]
                                if np.any(np.isnan(comp)):
                                    raise ValueError("Data contains nan's")
                                hist, edges = np.histogram(comp, range=datarange, bins=bins, density=False)
                                num_elems = len(comp)
                                avg = np.mean(comp)
                                # nicename2 = re.sub(r'[\W_]', ' ', str(cmkeys[0]))
                                nicename2 = str(cmkeys[0][:2])
                                ax.axvline(avg, color='grey', linestyle='dashed', linewidth=1)
                                bincentres = [(edges[i] + edges[i + 1]) / 2. for i in range(len(edges) - 1)]
                                widths = np.diff(edges)
                                norm = np.dot(hist, widths)
                                hist2 = hist / norm
                                hist2_divisor = hist2
                                hist2_divisor[hist2 == 0] = 1
                                hist_diff = (hist2 - hist1)
                                if normalized:
                                    hist_diff = hist_diff / hist2_divisor
                                ax.plot(bincentres, hist_diff, label='[' + nicename2 + '] - [' + nicename + ']', linewidth=0.5, color=color, marker='o',
                                        markersize=5)
                                # ax.step(bincentres, hist_diff, where='mid',
                                #         label='[' + nicename2 + '] - [' + nicename + ']', linewidth=1, color=color)
                                # ax.step(bincentres, np.cumsum(hist2-hist1),where='mid',label='cumsum (['+ nicename2 +'] - ['+ nicename + '])', linewidth=4,color=color,alpha=0.5)
                    if normalized:
                        ax.set_ylim([-0.9, 0.9])
                        ax.set_ylabel('$\\frac{P(S_\mathrm{ED}) - P(S_\mathrm{xDMRG})}{P(S_\mathrm{ED})} $')

                    else:
                        ax.set_ylim([-0.25, 0.25])
                        ax.set_ylabel('$P(S_\mathrm{ED}) - P(S_\mathrm{xDMRG})$')
                    ax.set_xlabel('$S_E$')
                    ax.set_title('$L = ' + str(chain_length) + '$')
                    # ax.set_xlim(1e-21,100)
                used_ax = used_ax + 1
                ax.legend()
            if normalized:
                fig.suptitle("Difference of distributions\n"
                             "No linearFit means no bias\n"
                             "$\Delta = " + str(delt) + '\lambda = ' + str(lamb) + "$"
                             )
            else:
                fig.suptitle("Difference of distributions (not normalized)\n"
                             "No linearFit means no bias\n"
                             "$\Delta = " + str(delt) + '\lambda = ' + str(lamb) + "$"
                             )
            for ax in np.ravel(axes)[used_ax:]:
                fig.delaxes(ax)
            if plotdir != '':
                if normalized:
                    plt.savefig(plotdir + '/S_distribution_diff_' + l + '_' + d + '.pdf', format='pdf')
                    plt.savefig(plotdir + '/S_distribution_diff_' + l + '_' + d + '.png', format='png')
                else:
                    plt.savefig(plotdir + '/S_distribution_diff_not_normalized_' + l + '_' + d + '.pdf', format='pdf')
                    plt.savefig(plotdir + '/S_distribution_diff_not_normalized_' + l + '_' + d + '.png', format='png')

            #
            #     key_sorted = sorted(basenode.keys(), key=natural_keys)
            #     for key in key_sorted:
            #         key_num = key_num+1
            #         node         = h5_src[L][l][d][key]['entanglement_entropies']
            #         # print("Num Lambdas: " , len(node['data'][0, :]))
            #         # print("Mid Lambda : " , int(len(node['data'][0, :])/2))
            #         # print("Middle : " , middle)
            #         for win_idx, win in enumerate(variance_window_limits):
            #             # First get the data off one of the given keys
            #             if np.shape(node['data'])[1] == chain_length + 1:
            #                 width = np.shape(node['data'])[1]
            #                 middle = int(width / 2)
            #                 idx = get_v_filtered_index_list(h5_src[L][l][d][key], win)
            #                 data = node['data'][idx, middle]
            #             elif np.shape(node['data'])[0] == chain_length - 1:
            #                 width = np.shape(node['data'])[0]
            #                 middle = int(width / 2)
            #                 data = node['data'][middle, :]
            #             else:
            #                 raise TypeError("Data shape couldnt be matched: [" + str(np.shape(node['data'])) + "]")
            #             if np.any(np.isnan(data)):
            #                 raise ValueError("Data contains nan's")
            #             color = variance_colors[win_idx]
            #             data_range = [np.min(data), np.max(data)]
            #             hist, edges = np.histogram(data,range=data_range, bins=bins, density=False)
            #             num_elems = len(data)
            #             avg = np.mean(data)
            #             nicename = re.sub(r'[\W_]', ' ', str(key))
            #             ax.axvline(avg, color=color, linestyle='dashed', linewidth=1, label='Avg ' + nicename + ' (' +str(num_elems)+ ')')
            #
            #             bincentres = [(edges[i] + edges[i + 1]) / 2. for i in range(len(edges) - 1)]
            #             widths = np.diff(edges)
            #             norm = np.dot(hist, widths)
            #             hist1 = hist/norm
            #
            #             # Now get the the data from the comparator key to subtract
            #             node_comp = h5_src[L][l][d][cmkeys[0]]['entanglement_entropies']
            #             if np.shape(node_comp['data'])[1] == chain_length + 1:
            #                 width = np.shape(node_comp['data'])[1]
            #                 middle = int(width / 2)
            #                 idx = get_v_filtered_index_list(h5_src[L][l][d][key], win)
            #                 data_comp = node_comp['data'][idx, middle]
            #             elif np.shape(node_comp['data'])[0] == chain_length - 1:
            #                 width = np.shape(node_comp['data'])[0]
            #                 middle = int(width / 2)
            #                 data_comp = node_comp['data'][middle, :]
            #             else:
            #                 raise TypeError("Data shape couldnt be matched: [" + str(np.shape(node_comp['data'])) + "]")
            #             if np.any(np.isnan(data_comp)):
            #                 raise ValueError("Data contains nan's")
            #             hist, edges = np.histogram(data_comp,range=data_range, bins=bins, density=False)
            #             num_elems = len(data_comp)
            #             avg = np.mean(data_comp)
            #             # nicename2 = re.sub(r'[\W_]', ' ', str(cmkeys[0]))
            #             nicename2 = str(cmkeys[0][:2])
            #             ax.axvline(avg, color='grey', linestyle='dashed', linewidth=1, label='Avg ' + nicename2 + ' (' +str(num_elems)+ ')')
            #             bincentres = [(edges[i] + edges[i + 1]) / 2. for i in range(len(edges) - 1)]
            #             widths = np.diff(edges)
            #             norm = np.dot(hist, widths)
            #             hist2 = hist/norm
            #             hist_diff = (hist2-hist1)/hist2
            #             ax.step(bincentres, hist_diff, where='mid', label='['+ nicename2 +'] - ['+ nicename + ']', linewidth=1,color=color)
            #             # ax.step(bincentres, np.cumsum(hist2-hist1),where='mid',label='cumsum (['+ nicename2 +'] - ['+ nicename + '])', linewidth=4,color=color,alpha=0.5)
            #         ax.set_ylim([-0.7,0.7])
            #         ax.set_xlabel('$S_E$')
            #         ax.set_ylabel('$\\frac{P(S_\mathrm{ED}) - P(S_\mathrm{xDMRG})}{P(S_\mathrm{ED})} $')
            #         ax.set_title('$L = ' + str(chain_length) + '$')
            #         # ax.set_xlim(1e-21,100)
            #     used_ax = used_ax + 1
            #     ax.legend()
            # fig.suptitle('Difference of distributions \nMid-chain entanglement entropy @ $\Delta = $' + str(delt) + '$\lambda = $' + str(lamb))
            # for ax in np.ravel(axes)[used_ax:]:
            #     fig.delaxes(ax)
            # if plotdir != '':
            #     Jh = re.sub('/', '_', str(d))
            #     plt.savefig(plotdir + '/S_distribution_diff_' + l + '_' + str(Jh) + '.pdf', format='pdf')
            #     plt.savefig(plotdir + '/S_distribution_diff_' + l + '_' + str(Jh) + '.png', format='png')
    h5close(h5_src)
