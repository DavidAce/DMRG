from .tools import *
from src.io.h5ops import *
import matplotlib.pyplot as plt
from .filter import *


def multiplot_S_vs_Delta(src, plotdir='', key_list=''):
    print('Plotting: Entanglement entropy vs Delta for: ', key_list)
    h5_src = h5open(src, 'r')
    path_L = h5py_unique_finder(h5_src, filter='L_', dep=1)
    path_l = h5py_unique_finder(h5_src, filter='l_', dep=2)
    path_J = h5py_unique_finder(h5_src, filter='J_', dep=3)
    path_h = h5py_unique_finder(h5_src, filter='h_', dep=4)
    path_d = []
    for J in path_J:
        for h in path_h:
            path_d.append(J + '/' + h)
    for l in path_l:
        # In each figure we want one subplot per unique_L
        rows, cols = get_optimal_subplot_num(len(path_L))
        fig, axes = plt.subplots(nrows=rows, ncols=cols, figsize=(7 * cols, 7 * rows))
        fig.tight_layout(pad=5, w_pad=1.0, h_pad=1.0)
        fig.subplots_adjust(wspace=0.3, hspace=0.3)
        used_ax = 0
        delt = 0
        lamb = 0
        for ax, L in zip(np.ravel(axes), path_L):

            key_num = 0
            chain_length = 0
            max_S = 0
            stats = {}
            for d in path_d:
                h5keys = [item for item in h5_src[L][l][d].keys() if any(s in item for s in key_list)]
                key_sorted = sorted(h5keys, key=natural_keys)
                for key in key_sorted:
                    stats[key] = {'x': [], 'y': [], 'e': [], 'n': []}

            for d in path_d:
                h5keys = [item for item in h5_src[L][l][d].keys() if any(s in item for s in key_list)]
                key_sorted = sorted(h5keys, key=natural_keys)
                for key in key_sorted:
                    try:
                        node = h5_src[L][l][d][key]['entanglement_entropies']
                    except:
                        # print('[ simulation_time ] missing in key: ', key)
                        continue
                    key_num = key_num + 1
                    chain_length = node.attrs['chain_length']
                    lamb = node.attrs['lambda']
                    delta = node.attrs['delta']
                    if np.shape(node['data'])[1] == chain_length + 1:
                        width = np.shape(node['data'])[1]
                        middle = int(width / 2)
                        idx = get_v_filtered_index_list(h5_src[L][l][d][key], variance_window_limits[0])
                        data = node['data'][idx, middle]
                    elif np.shape(node['data'])[0] == chain_length - 1:
                        width = np.shape(node['data'])[0]
                        middle = int(width / 2)
                        data = node['data'][middle, :]
                    else:
                        raise TypeError("Data shape couldnt be matched: [" + str(np.shape(node['data'])) + "]")
                    if np.any(np.isnan(data)):
                        raise ValueError("Data contains nan's")

                    # middle = int(chain_length / 2)
                    # idx = get_v_filtered_index_list(h5_src[L][l][d][key], variance_window_limits[0])
                    # data = node['data'][idx, middle]
                    # if np.any(np.isnan(data)):
                    #     raise ValueError("Data contains nan's")
                    num_elems = len(data)
                    avg = np.nanmean(data, axis=0)
                    ste = np.nanstd(data, axis=0) / np.sqrt(num_elems)
                    stats[key]['x'].append(chain_length)
                    stats[key]['y'].append(avg)
                    stats[key]['e'].append(ste)
                    stats[key]['n'].append(num_elems)
            key_sorted = sorted(stats.keys(), key=natural_keys)
            for key in key_sorted:
                val = stats[key]
                nicekey = re.sub(r'[\W_]', ' ', str(key))
                ax.errorbar(x=val['x'], y=val['y'], yerr=val['e'], label=nicekey, capsize=2, elinewidth=0, markeredgewidth=0.8)
                max_S = np.max([max_S, np.max(val['y'])])
                # for i, txt in enumerate(val['n']):
                #     ax.annotate(txt, (val['x'][i], val['y'][i]))
            ax.legend()
            ax.set_xlabel('$\Delta$')
            ax.set_ylabel('$\langle S_E(L/2)  \\rangle$')
            ax.set_title('$L=' + str(chain_length) + '$')
            if not np.isnan(max_S):
                ax.set_ylim(bottom=0, top=max_S * 1.2)
            used_ax = used_ax + 1
        fig.suptitle('Mid-chain entanglement entropy vs $\Delta$ @ $\lambda = ' + str(lamb) + '$')
        for ax in np.ravel(axes)[used_ax:]:
            fig.delaxes(ax)
        if plotdir != '':
            plt.savefig(plotdir + '/S_vs_Delta_' + l + '.pdf', format='pdf')
            plt.savefig(plotdir + '/S_vs_Delta_' + l + '.png', format='png')

    h5close(h5_src)
