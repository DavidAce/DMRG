from .tools import *
from src.io.h5ops import *
import matplotlib.pyplot as plt
from .filter import *
from tqdm import tqdm


def multiplot_Sq_vs_Length(src, key_list, plotdir=''):
    print('Plotting: Entanglement entropy vs Length')
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
        rows, cols = get_optimal_subplot_num(len(path_d))
        fig, axes = plt.subplots(nrows=rows, ncols=cols, figsize=(7 * cols, 7 * rows))
        fig.tight_layout(pad=5, w_pad=1.0, h_pad=1.0)
        fig.suptitle('Mid-chain entanglement entropy vs chain length')
        fig.subplots_adjust(wspace=0.3, hspace=0.3)
        lamb = None
        delta = None
        for id, (ax, d) in enumerate(zip(np.ravel(axes), path_d)):
            stats = {}
            for L in path_L:
                h5keys = [item for item in h5_src[L][l][d].keys() if any(s in item for s in key_list)]
                key_sorted = sorted(h5keys, key=natural_keys)
                for key in key_sorted:
                    stats[key] = {'x': [], 'y': [], 'e': [], 'n': []}
            for L in path_L:
                h5keys = [item for item in h5_src[L][l][d].keys() if any(s in item for s in key_list)]
                key_sorted = sorted(h5keys, key=natural_keys)
                for key in key_sorted:
                    node = h5_src[L][l][d][key]['entanglement_entropies']
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
                    num_elems = len(data)
                    avg = np.nanmean(data, axis=0)
                    std = np.nanstd(data, axis=0)
                    ste = np.nanstd(data, axis=0) / np.sqrt(num_elems)
                    med = np.nanmedian(data, axis=0)
                    stats[key]['x'].append(chain_length)
                    stats[key]['y'].append(avg)
                    stats[key]['e'].append(ste)
                    stats[key]['n'].append(num_elems)

            key_sorted = sorted(stats.keys(), key=natural_keys)
            for key in key_sorted:
                val = stats[key]
                nicekey = re.sub(r'[\W_]', ' ', str(key))
                ax.errorbar(x=val['x'], y=val['y'], yerr=val['e'], label=nicekey, capsize=2, elinewidth=0.3, markeredgewidth=0.8)
                for i, txt in enumerate(val['n']):
                    ax.annotate(txt, (val['x'][i], val['y'][i]))
            ax.set_ylabel('$\langle S (L/2) \\rangle$')
            ax.set_xlabel('Length')
            ax.set_title('$\lambda = ' + str(lamb) + ' \qquad \Delta = ' + str(delta) + '$')
            ax.legend()

        for ax in np.ravel(axes)[len(path_d):]:
            fig.delaxes(ax)

        if plotdir != '':
            plt.savefig(plotdir + '/S_vs_Length_' + l + '.pdf', format='pdf')
            plt.savefig(plotdir + '/S_vs_Length_' + l + '.png', format='png')
    h5close(h5_src)
