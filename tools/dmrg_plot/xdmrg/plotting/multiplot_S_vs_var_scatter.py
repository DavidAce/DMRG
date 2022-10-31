from src.plotting.tools import *
from src.io.h5ops import *
from src.plotting.filter import *
import numpy as np
import matplotlib.pyplot as plt


def multiplot_S_vs_var_scatter(src, key_list, plotdir=''):
    print('Plotting: Histogram of midchain Entanglement entropy')
    h5_src = h5open(src, 'r')
    path_L = h5py_unique_finder(h5_src, filter='L_', dep=1)
    path_l = h5py_unique_finder(h5_src, filter='l_', dep=2)
    path_J = h5py_unique_finder(h5_src, filter='J_', dep=3)
    path_h = h5py_unique_finder(h5_src, filter='h_', dep=4)
    path_d = []
    for J in path_J:
        for h in path_h:
            path_d.append(J + '/' + h)

    # One figure per unique_l, unique_J and unique_h
    for l in path_l:
        for d in path_d:
            # In each figure we want one subplot per unique_L
            rows, cols = get_optimal_subplot_num(len(path_L))
            fig, axes = plt.subplots(nrows=rows, ncols=cols, figsize=(3.5 * cols, 3.5 * rows))
            fig.tight_layout(pad=5, w_pad=1.0, h_pad=1.0)
            fig.subplots_adjust(wspace=0.3, hspace=0.3)
            used_ax = 0
            delt = 0
            lamb = 0
            for ax, L in zip(np.ravel(axes), path_L):
                key_num = 0
                key_sorted = sorted(h5_src[L][l][d].keys(), key=natural_keys)
                key_sorted = [item for item in key_sorted if key_list in key_sorted]
                for key in key_sorted:
                    key_num = key_num + 1
                    node_e = h5_src[L][l][d][key]['entanglement_entropies']
                    node_v = h5_src[L][l][d][key]['variance']
                    num = node_e['num'][()]
                    length = node_e.attrs['chain_length']
                    delt = node_e.attrs['delta']
                    lamb = node_e.attrs['lambda']

                    assert (num == node_v['num'][()])
                    assert (length == node_v.attrs['chain_length'])
                    assert (delt == node_v.attrs['delta'])
                    assert (lamb == node_v.attrs['lambda'])
                    var_data = node_v['data']
                    ent_data = node_e['data'][:, int(length / 2)]
                    nicename = re.sub(r'[\W_]', ' ', str(key))
                    nicename = nicename + ' (' + str(num) + ')'
                    ax.scatter(var_data, ent_data, s=3, label=nicename)
                    ax.set_xlabel('$\log \sigma(E)^2$')
                    ax.set_ylabel('$S_E$')
                    ax.set_xscale('log')
                    ax.set_title('$L = ' + str(length) + '$')
                used_ax = used_ax + 1
                for win_idx, win in enumerate(variance_window_limits):
                    for lim_idx, lim in enumerate(win):
                        ax.axvline(lim, color=variance_colors[win_idx], linestyle='dashed', linewidth=1,
                                   label=variance_window_names[win_idx][lim_idx])
                ax.legend()

            fig.suptitle('Distribution of mid-chain entanglement entropy @ $\Delta = $' + str(delt) + '$\lambda = $' + str(lamb))
            for ax in np.ravel(axes)[used_ax:]:
                fig.delaxes(ax)
            if plotdir != '':
                Jh = re.sub('/', '_', str(d))
                plt.savefig(plotdir + '/S_vs_Var_scatter_' + l + '_' + str(Jh) + '.pdf', format='pdf')
    h5close(h5_src)
