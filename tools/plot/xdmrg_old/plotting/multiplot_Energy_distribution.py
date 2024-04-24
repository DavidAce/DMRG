from .tools import *
from src.io.h5ops import *
import matplotlib.pyplot as plt
from .filter import *


def multiplot_Energy_distribution(src, plotdir='', key_list='', type='typical'):
    print('Plotting: Energy distribution')
    h5_src = h5open(src, 'r')
    path_L = h5py_unique_finder(h5_src, filter='L_', dep=1)
    path_l = h5py_unique_finder(h5_src, filter='l_', dep=2)
    path_J = h5py_unique_finder(h5_src, filter='J_', dep=3)
    path_h = h5py_unique_finder(h5_src, filter='h_', dep=4)
    path_d = []
    for J in path_J:
        for h in path_h:
            path_d.append(J + '/' + h)
    timecounter = 0
    simcounter = 0
    # One figure per unique_l, unique_J and unique_h
    try:
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
                    key_num = 0
                    h5keys = [item for item in h5_src[L][l][d].keys() if any(s in item for s in key_list)]
                    key_sorted = sorted(h5keys, key=natural_keys)
                    for key in key_sorted:
                        try:
                            node = h5_src[L][l][d][key]['energy']
                        except Exception as er:
                            continue
                            # raise LookupError('[ simulation_time ] missing in key: ' + str(key) + '\nReason: ' + str(er))
                        key_num = key_num + 1
                        chain_length = node.attrs['chain_length']
                        delt = node.attrs['delta']
                        lamb = node.attrs['lambda']
                        if 'ed' in key:
                            data = np.array(node['data']).T / chain_length
                        else:
                            idx = get_v_filtered_index_list(h5_src[L][l][d][key], variance_window_limits[0])
                            data = np.array(node['data'])[idx]
                        if np.any(np.isnan(data)):
                            raise ValueError("Data contains nan's")

                        # idx = get_v_filtered_index_list(h5_src[L][l][d][key], variance_window_limits[0])
                        # data   = np.array(node['data'])
                        # middle = int((len(data[0]) - 1) / 2)
                        hist, edges = np.histogram(data, bins=100, density=True)
                        num = np.shape(data)[0]
                        avg = np.mean(data)
                        nicename = re.sub(r'[\W_]', ' ', str(key))
                        nicename = nicename + ' (' + str(num) + ')'
                        bincentres = [(edges[i] + edges[i + 1]) / 2. for i in range(len(edges) - 1)]
                        widths = np.diff(edges)
                        ax.bar(bincentres, hist, width=widths, align='center', label=nicename, alpha=0.2 + 0.7 / key_num, linewidth=0.2)
                        ax.axvline(avg, color='grey', linestyle='dashed', linewidth=1)
                        timecounter = timecounter + np.sum(data)
                        simcounter = simcounter + num

                    ax.set_xlabel('Energy')
                    ax.set_ylabel('Histogram')
                    ax.set_title('$L = ' + str(chain_length) + '$')
                    # ax.set_xlim(left=0, right=1)
                    # ax.set_yscale('log')
                    ax.legend()
                    used_ax = used_ax + 1
                fig.suptitle('Energy density distribution @ $\Delta = ' + str(delt) + '\quad \lambda = ' + str(lamb) + '$')
                for ax in np.ravel(axes)[used_ax:]:
                    fig.delaxes(ax)
                if plotdir != '':
                    Jh = re.sub('/', '_', str(d))
                    plt.savefig(plotdir + '/energy_distribution_' + l + '_' + str(Jh) + '.pdf', format='pdf')
                    plt.savefig(plotdir + '/energy_distribution_' + l + '_' + str(Jh) + '.png', format='png')
    except Exception as er:
        print("Could not plot Energy. \nReason: " + str(er))
    h5close(h5_src)
