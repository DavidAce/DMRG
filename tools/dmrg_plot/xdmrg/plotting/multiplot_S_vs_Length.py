from src.plotting.tools import *
from src.io.h5ops import *
import matplotlib.pyplot as plt
from src.plotting.filter import *
from tqdm import tqdm


def multiplot_S_vs_Length(src, plotdir='', algo_inc='', state_inc='', type='average'):
    print('Plotting: Entanglement entropy vs Length for: ', algo_inc, state_inc)
    h5_src_others = []
    if isinstance(src, list):
        h5_src = h5open(src[0], 'r')
        for i in range(1, len(src)):
            h5_src_others[i] = h5open(src[i], 'r')
    else:
        h5_src = h5open(src, 'r')
    path_L = h5py_unique_finder(h5_src, keypattern='L_', dep=1)
    path_l = h5py_unique_finder(h5_src, keypattern='l_', dep=2)
    path_d = h5py_unique_finder(h5_src, keypattern='d_', dep=3)

    rows, cols = get_optimal_subplot_num(len(path_d))
    fig, axes = plt.subplots(nrows=rows, ncols=cols, figsize=(7 * cols, 7 * rows), sharey='all',
                             sharex='all')
    fig.tight_layout(pad=5, w_pad=1.0, h_pad=1.0)
    fig.suptitle('Mid-chain entanglement entropy vs chain length')
    fig.subplots_adjust(wspace=0.2, hspace=0.2)
    axes_used = []
    max_S = 0
    min_S = 10
    for id, (ax, d) in enumerate(zip(np.ravel(axes), path_d)):
        ed_palette = itertools.cycle(sns.color_palette("Set2"))
        current_palette = itertools.cycle(sns.color_palette())
        lamb = None
        delta = None
        plotdict = {}
        plot_id = 0
        for L_id, (L_key, L_node) in enumerate(h5_src.items()):
            for l_id, (l_key, l_node) in enumerate(L_node.items()):
                if not d in l_node:
                    continue
                basenode = l_node[d]
                delta = basenode.attrs['delta']
                for algo_id, (algokey, algopath, algonode) in enumerate(
                        h5py_group_iterator(node=basenode, keypattern=algo_inc, dep=1)):
                    for state_id, (statekey, statepath, statenode) in enumerate(
                            h5py_group_iterator(node=algonode, keypattern=state_inc, dep=1)):
                        state_id = state_id + 1
                        for v_win_idx, v_win in enumerate(variance_window_limits):
                            for e_win_idx, e_win in enumerate(energy_window_limits):
                                idx = get_v_e_filtered_index_list(statenode, v_win, e_win)
                                for datakey, datapath, datanode in h5py_node_finder(node=statenode,
                                                                                    keypattern='entanglement_entropy_midchain',
                                                                                    dep=8):
                                    if "checkpoint" in datapath:
                                        continue
                                    plotkey = "{}/{}/{}/{}/{}".format(l_key, algokey, statekey, v_win_idx, e_win_idx)
                                    ndata = datanode['num'][()]
                                    if not plotkey in plotdict:
                                        plot_id = plot_id + 1
                                        plotdict[plotkey] = {'x': [], 'y': [], 'e': [], 'n': []}
                                        plotdict[plotkey]['lambda'] = basenode.attrs['lambda']
                                        plotdict[plotkey]['delta'] = basenode.attrs['delta']
                                        plotdict[plotkey]['algo_id'] = algo_id
                                        plotdict[plotkey]['state_id'] = state_id
                                        plotdict[plotkey]['plot_id'] = plot_id

                                        if "states" in statekey:
                                            plotdict[plotkey]['color'] = next(ed_palette)
                                            plotdict[plotkey]['nicename'] = "ED e=[" + statenode.attrs[
                                                "efmt"] + "] $\lambda =" + str(basenode.attrs['lambda']) + "$"
                                            plotdict[plotkey]['lwidth'] = 3.0
                                            plotdict[plotkey]['lalpha'] = 1.0
                                            plotdict[plotkey]['mstyle'] = '^'
                                            plotdict[plotkey]['lstyle'] = 'solid'
                                        else:
                                            plotdict[plotkey]['color'] = next(current_palette)
                                            plotdict[plotkey]['nicename'] = re.sub(r'[\W_]', ' ',
                                                                                   str(algokey + " " + statekey)) \
                                                                            + " $\lambda =" + str(
                                                basenode.attrs['lambda']) + "$"
                                            plotdict[plotkey]['lwidth'] = 1.4
                                            plotdict[plotkey]['lalpha'] = 1.0
                                            plotdict[plotkey]['mstyle'] = '.'
                                            plotdict[plotkey]['lstyle'] = 'dotted'

                                    plotdict[plotkey]['x'].append(basenode.attrs['model_size'])
                                    plotdict[plotkey]['y'].append(datanode['avg'][()])
                                    plotdict[plotkey]['e'].append(datanode['ste'][()])
                                    plotdict[plotkey]['n'].append(ndata)

        for key, val in plotdict.items():
            ax.errorbar(x=val['x'], y=val['y'], yerr=val['e'], label=val['nicename'], capsize=2, elinewidth=0.3,
                        markeredgewidth=0.8,
                        color=val['color'], marker=val['mstyle'], linestyle=val['lstyle'], alpha=val['lalpha'],
                        linewidth=val['lwidth'])
            max_S = np.max([max_S, np.max(val['y'])])
            min_S = np.min([min_S, np.min(val['y'])])
            axes_used.append(id) if not id in axes_used else axes_used
        min_S = 0.54
        max_S = 0.65
        for key, val in plotdict.items():
            for i, txt in enumerate(val['n']):
                ytext = min_S - 0.06 * (max_S - min_S) * val['plot_id']
                ax.annotate(txt, (val['x'][i], val['y'][i]), textcoords='data', xytext=[val['x'][i] - 0.2, ytext],
                            color=val['color'], fontsize='x-small')

        ax.set_ylabel('$\langle S (L/2) \\rangle$')
        ax.set_xlabel('Length')
        ax.set_title('$\Delta = ' + str(delta) + '$')

    prettify_plot(fig, axes, rows=rows, cols=cols, axes_used=axes_used, ymin=min_S * 0.9, ymax=max_S * 1.05)
    #
    # for ax in np.ravel(axes)[len(path_d):]:
    #     fig.delaxes(ax)
    # for ax in np.ravel(axes):
    #     ax.set_ylim(ymin=min_S * 0.7, ymax=max_S * 1.1)

    if plotdir != '':
        plt.savefig(plotdir + '/S_vs_Length.pdf', format='pdf')
        plt.savefig(plotdir + '/S_vs_Length.png', format='png')

    h5close(h5_src)


def multiplot_S_vs_Length_old(src, key_list, plotdir='', algo_inc='', state_inc='', type='average'):
    print('Plotting: Entanglement entropy vs Length for: ', algo_inc, state_inc)
    h5_src = h5open(src, 'r')
    path_L = h5py_unique_finder(h5_src, keypattern='L_', dep=1)
    path_l = h5py_unique_finder(h5_src, keypattern='l_', dep=2)
    path_J = h5py_unique_finder(h5_src, keypattern='J_', dep=3)
    path_h = h5py_unique_finder(h5_src, keypattern='h_', dep=4)
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
                        idx = get_v_inced_index_list(h5_src[L][l][d][key], variance_window_limits[0])
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
                ax.errorbar(x=val['x'], y=val['y'], yerr=val['e'], label=nicekey, capsize=2, elinewidth=0.3,
                            markeredgewidth=0.8)
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
