from src.plotting.tools import *
from src.io.h5ops import *
from src.plotting.filter import *
import numpy as np
import matplotlib.pyplot as plt


def multiplot_S_vs_ewin_scatter(src, plotdir='', key_list=''):
    print('Plotting: Cumulative avg of midchain Entanglement entropy for: ', key_list)
    h5_src = h5open(src, 'r')
    path_L = h5py_unique_finder(h5_src, filter='L_', dep=1)
    path_l = h5py_unique_finder(h5_src, filter='l_', dep=2)
    path_J = h5py_unique_finder(h5_src, filter='J_', dep=3)
    path_h = h5py_unique_finder(h5_src, filter='h_', dep=4)
    path_d = []

    S_ED = {'L_16': 0.5965248564095034, 'L_20': 0.6180112356038007, 'L_24': 0.6299770044181504, 'L_28': 0.6436141776418879, 'L_32': 0.6519452957939131,
            'L_36': 0.6581362571591568}
    S_ER = {'L_16': 0.0013482694624003474, 'L_20': 0.0013576325487211977, 'L_24': 0.0013648335485899527, 'L_28': 0.0013881518030119846,
            'L_32': 0.0013917464514767426, 'L_36': 0.0014060116906587084}

    for J in path_J:
        for h in path_h:
            path_d.append(J + '/' + h)

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
                max_S = 0
                h5keys = [item for item in h5_src[L][l][d].keys() if any(s in item for s in key_list)]
                key_sorted = sorted(h5keys, key=natural_keys)
                for key in key_sorted:
                    s_avg = []
                    s_ste = []
                    e_mid = []
                    for ewin in ewin_array:
                        node = h5_src[L][l][d][key]['entanglement_entropies']
                        length = node.attrs['chain_length']
                        middle = int(length / 2)

                        idx = get_e_filtered_index_list(h5_src[L][l][d][key], ewin)
                        if len(idx) == 0:
                            continue
                        idx = idx[idx < len(node['data'])]  # Filter out indices larger than current data
                        data = np.array(node['data'][idx, middle])
                        if np.any(np.isnan(data)):
                            raise ValueError("Data contains nan's")
                        s_avg.append(np.mean(data, axis=0))
                        s_ste.append(np.std(data, axis=0) / np.sqrt(len(idx)))
                        e_mid.append(np.mean(ewin))
                    if np.size(s_avg) == 0:
                        continue
                    nicename = re.sub(r'[\W_]', ' ', str(key))
                    # ax.scatter(e_mid,s_avg, s=3, label=nicename)
                    ax.errorbar(e_mid, s_avg, yerr=s_ste, linestyle="None", fmt='.', capsize=0.1, elinewidth=0.3, markeredgewidth=0.8, label=nicename)
                    node_e = h5_src[L][l][d][key]['entanglement_entropies']
                    length = node_e.attrs['chain_length']
                    delt = node_e.attrs['delta']
                    lamb = node_e.attrs['lambda']
                    ax.set_xlabel('Energy center')
                    ax.set_ylabel('$\langle S_E \\rangle$')
                    # ax.set_xscale('log')
                    ax.set_title('$L = ' + str(length) + '$')
                    if L in S_ED.keys():
                        ax.fill_between([0, 1], S_ED[L] - S_ER[L], S_ED[L] + S_ER[L], color='grey', alpha=0.5, linewidth=1.4, label='ED')
                    max_S = np.max([max_S, np.max(s_avg)])
                    if not np.isnan(max_S):
                        ax.set_ylim(0, max_S * 1.2)
                        ax.set_xlim(ewin_array[0][0], ewin_array[-1][1])
                used_ax = used_ax + 1
                # for win_idx, win in enumerate(variance_window_limits):
                #     for lim_idx, lim in enumerate(win):
                #         ax.axvline(lim, color=variance_colors[win_idx], linestyle='dashed', linewidth=1,
                #                    label=variance_window_names[win_idx][lim_idx])
                ax.legend()

            fig.suptitle('Average of mid-chain entanglement entropy\nin energy windows of width $w = $' + str(winsize) + '\n@ $\Delta = $' + str(
                delt) + '$\lambda = $' + str(lamb))
            for ax in np.ravel(axes)[used_ax:]:
                fig.delaxes(ax)
            if plotdir != '':
                Jh = re.sub('/', '_', str(d))
                plt.savefig(plotdir + '/S_vs_Ewin_scatter_' + l + '_' + str(Jh) + '.pdf', format='pdf')
                plt.savefig(plotdir + '/S_vs_Ewin_scatter_' + l + '_' + str(Jh) + '.png', format='png')
    h5close(h5_src)
