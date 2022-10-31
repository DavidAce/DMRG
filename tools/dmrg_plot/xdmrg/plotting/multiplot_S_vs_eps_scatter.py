from src.plotting.tools import *
from src.io.h5ops import *
from src.plotting.filter import *
import numpy as np
import matplotlib.pyplot as plt


def multiplot_S_vs_eps_scatter(src, plotdir='', key_list=''):
    print('Plotting: Cumulative avg of midchain Entanglement entropy')
    h5_src = h5open(src, 'r')
    path_L = h5py_unique_finder(h5_src, filter='L_', dep=1)
    path_l = h5py_unique_finder(h5_src, filter='l_', dep=2)
    path_J = h5py_unique_finder(h5_src, filter='J_', dep=3)
    path_h = h5py_unique_finder(h5_src, filter='h_', dep=4)
    path_d = []
    for J in path_J:
        for h in path_h:
            path_d.append(J + '/' + h)
    S_ED = {'L_16': 0.5965248564095034, 'L_20': 0.6180112356038007, 'L_24': 0.6299770044181504, 'L_28': 0.6436141776418879, 'L_32': 0.6519452957939131,
            'L_36': 0.6581362571591568}
    S_ER = {'L_16': 0.0013482694624003474, 'L_20': 0.0013576325487211977, 'L_24': 0.0013648335485899527, 'L_28': 0.0013881518030119846,
            'L_32': 0.0013917464514767426, 'L_36': 0.0014060116906587084}

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
                    s_avg1 = []
                    s_ste1 = []
                    s_avg2 = []
                    s_ste2 = []
                    s_avg3 = []
                    s_ste3 = []
                    for eps in eps_array:
                        node = h5_src[L][l][d][key]['entanglement_entropies']
                        chain_length = node.attrs['chain_length']

                        if np.shape(node['data'])[1] == chain_length + 1:
                            width = np.shape(node['data'])[1]
                            middle = int(width / 2)
                            idx = get_v_filtered_index_list(h5_src[L][l][d][key], [0, eps])
                            if (len(idx) > 0):
                                data = np.array(node['data'][idx, middle])
                                if np.any(np.isnan(data)):
                                    raise ValueError("Data1 contains nan's")
                                s_avg1.append(np.mean(data, axis=0))
                                s_ste1.append(np.std(data, axis=0) / np.sqrt(len(idx)))
                            else:
                                s_avg1.append(0)
                                s_ste1.append(0)
                            idx = get_v_filtered_index_list(h5_src[L][l][d][key], [eps, 1])
                            if (len(idx) > 0):
                                data = np.array(node['data'][idx, middle])
                                if np.any(np.isnan(data)):
                                    raise ValueError("Data2 contains nan's")
                                s_avg2.append(np.mean(data, axis=0))
                                s_ste2.append(np.std(data, axis=0) / np.sqrt(len(idx)))
                            else:
                                s_avg2.append(0)
                                s_ste2.append(0)
                            eps_min = eps / 10
                            eps_max = eps * 10
                            idx = get_v_filtered_index_list(h5_src[L][l][d][key], [eps_min, eps_max])
                            if (len(idx) > 0):
                                data = np.array(node['data'][idx, middle])
                                if np.any(np.isnan(data)):
                                    raise ValueError("Data3 contains nan's")
                                s_avg3.append(np.mean(data, axis=0))
                                s_ste3.append(np.std(data, axis=0) / np.sqrt(len(idx)))
                            else:
                                s_avg3.append(0)
                                s_ste3.append(0)

                        elif np.shape(node['data'])[0] == chain_length - 1:
                            raise TypeError("This doesn't make sense for ed data")
                        else:
                            raise TypeError("Data shape couldnt be matched: [" + str(np.shape(node['data'])) + "]")
                        if np.any(np.isnan(data)):
                            raise ValueError("Data contains nan's")

                        # middle = int(chain_length / 2)
                        # idx     = get_v_filtered_index_list(h5_src[L][l][d][key], [0 , eps])
                        # if (len(idx) > 0):
                        #     data    = np.array(node['data'][idx, middle])
                        #     if np.any(np.isnan(data)):
                        #         raise ValueError("Data1 contains nan's")
                        #     s_avg1.append(np.mean(data, axis=0))
                        #     s_ste1.append(np.std(data,axis=0) / np.sqrt(len(idx)))
                        # else:
                        #     s_avg1.append(0)
                        #     s_ste1.append(0)

                        # idx = get_v_filtered_index_list(h5_src[L][l][d][key], [eps, 1])
                        # if (len(idx) > 0):
                        #     data    = np.array(node['data'][idx, middle])
                        #     if np.any(np.isnan(data)):
                        #         raise ValueError("Data2 contains nan's")
                        #     s_avg2.append(np.mean(data, axis=0))
                        #     s_ste2.append(np.std(data,axis=0) / np.sqrt(len(idx)))
                        # else:
                        #     s_avg2.append(0)
                        #     s_ste2.append(0)
                        # eps_min = eps/10
                        # eps_max = eps*10
                        # idx = get_v_filtered_index_list(h5_src[L][l][d][key], [eps_min, eps_max])
                        # if (len(idx) > 0):
                        #     data    = np.array(node['data'][idx, middle])
                        #     if np.any(np.isnan(data)):
                        #         raise ValueError("Data3 contains nan's")
                        #     s_avg3.append(np.mean(data, axis=0))
                        #     s_ste3.append(np.std(data,axis=0) / np.sqrt(len(idx)))
                        # else:
                        #     s_avg3.append(0)
                        #     s_ste3.append(0)

                    nicename = re.sub(r'[\W_]', ' ', str(key))
                    ax.errorbar(eps_array, s_avg1, yerr=s_ste1, linestyle="None", fmt='.', capsize=0.1, elinewidth=0.3, markeredgewidth=0.8,
                                label=nicename + " ]$-\infty,\epsilon$]")
                    ax.errorbar(eps_array, s_avg2, yerr=s_ste2, linestyle="None", fmt='.', capsize=0.1, elinewidth=0.3, markeredgewidth=0.8,
                                label=nicename + " [$\epsilon, 1$]")
                    ax.errorbar(eps_array, s_avg3, yerr=s_ste3, linestyle="None", fmt='.', capsize=0.1, elinewidth=0.3, markeredgewidth=0.8,
                                label=nicename + " [$\epsilon/10, \epsilon*10$]")
                    # ax.scatter(eps_array,s_avg1, s=3, label=nicename + " ]$-\infty,\epsilon$]")
                    # ax.scatter(eps_array,s_avg2, s=3, label=nicename + " [$epsilon, 1$]")
                    node_e = h5_src[L][l][d][key]['entanglement_entropies']
                    chain_length = node_e.attrs['chain_length']
                    delt = node_e.attrs['delta']
                    lamb = node_e.attrs['lambda']
                    ax.set_xlabel('$\epsilon$')
                    ax.set_ylabel('$S_E$')
                    ax.set_xscale('log')
                    ax.set_title('$L = ' + str(chain_length) + '$')
                    max_S = np.max([max_S, np.max(s_avg1), np.max(s_avg2)])
                    if not np.isnan(max_S):
                        ax.set_ylim(0, max_S * 1.2)
                        ax.set_xlim(1e-18, 10)
                    if L in S_ED.keys():
                        xlim = ax.get_xlim()
                        ax.fill_between(xlim, S_ED[L] - S_ER[L], S_ED[L] + S_ER[L], color='grey', alpha=0.5, linewidth=1.4, label='ED')

                used_ax = used_ax + 1
                # for win_idx, win in enumerate(variance_window_limits):
                #     for lim_idx, lim in enumerate(win):
                #         ax.axvline(lim, color=variance_colors[win_idx], linestyle='dashed', linewidth=1,
                #                    label=variance_window_names[win_idx][lim_idx])
                ax.legend()

            fig.suptitle('Average mid-chain entanglement entropy after filtering variance,  @ $\Delta = $' + str(delt) + '$\lambda = $' + str(lamb))
            for ax in np.ravel(axes)[used_ax:]:
                fig.delaxes(ax)
            if plotdir != '':
                Jh = re.sub('/', '_', str(d))
                plt.savefig(plotdir + '/S_vs_Eps_scatter_' + l + '_' + str(Jh) + '.pdf', format='pdf')
                plt.savefig(plotdir + '/S_vs_Eps_scatter_' + l + '_' + str(Jh) + '.png', format='png')
    h5close(h5_src)
