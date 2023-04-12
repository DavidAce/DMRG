from .tools import *
import matplotlib.pyplot as plt
# from src.database.database import *
from src.general.filter import *
import itertools


def multiplot_var_distribution(h5, db, meta, vwin):
    print('Plotting: Energy Variance distribution')
    algoinc = meta.get('inc').get('algo')
    stateinc = meta.get('inc').get('state')
    pointinc = meta.get('inc').get('point')

    # One figure per L
    # One subplot per delta
    # One line per lambda per state
    for key0 in db['keys']['L']:
        figrows, figcols = get_optimal_subplot_num(1 + len(db['keys']['d']))  # Add one for the legend
        fig, axes = plt.subplots(nrows=figrows, ncols=figcols, figsize=(4 * figcols, 4 * figrows), sharey='all', sharex='all')
        fig.tight_layout(pad=5, w_pad=1.0, h_pad=1.0)
        fig.subplots_adjust(wspace=0.2, hspace=0.2)
        axes_used = []
        vmin = None
        vmax = None
        hmax = None
        for deltaidx, (key1, ax) in enumerate(zip(db['keys']['d'], np.ravel(axes))):
            ed_palette = itertools.cycle(sns.color_palette("Set2"))
            numcolors = len(db['keys']['state']) * len(db['keys']['l'])
            current_palette = itertools.cycle(sns.color_palette("colorblind", numcolors))
            lstyles = itertools.cycle(['-.', '-', '--', ':', ])
            mstyles = itertools.cycle(('.', ',', '+', 'o', '*'))
            delt = None
            size = None
            for algokey in iterate(db['keys']['algo'], algoinc):
                for statekey in iterate(db['keys']['state'], stateinc):
                    for pointkey in iterate(db['keys']['point'], pointinc):
                        findlist = [key0, key1, algokey, statekey, pointkey, meta['groupname']]
                        datapaths = [val['path']['data'] for key, val in db['dsets'].items() if
                                     all(findlist[-1] in key and '/' + k + '/' in key for k in findlist[:-1])]

                        if not datapaths:
                            continue

                        lstyle = next(lstyles)
                        mstyle = next(mstyles)
                        # Now we have a set of setkeys with fixed L,d,algo and state, varying l.
                        for datapath in datapaths:
                            dset = db['dsets'][datapath]
                            midx = dset['vals']['midx']
                            lamb = dset['vals']['l']
                            delt = dset['vals']['d']
                            size = dset['vals']['L']
                            b = dset['vals']['b']
                            x = dset['vals']['x']
                            ndata = dset['num']
                            style = dset['style']
                            if "states" in statekey:
                                continue
                            data = h5[datapath]['data'][meta['colname']][()]
                            # data = dset['datanode']['data'][meta['colname']][()]
                            bins = np.logspace(start=-18, stop=0, num=64, endpoint=True)
                            hist, edges = np.histogram(data, bins=bins, density=False)
                            nicename = ''
                            if 'algo' in meta['legendcols']:
                                nicename += ' {}'.format(algokey)
                            if 'state' in meta['legendcols']:
                                nicename += ' {}'.format(statekey.replace('_', ''))
                            if 'b' in meta['legendcols']:
                                nicename += ' {}'.format(dset['tex']['eqs']['b'])
                            if 'x' in meta['legendcols']:
                                nicename += ' {}'.format(dset['tex']['eqs']['x'])
                            if 'time' in meta['legendcols'] and 'time' in dset['tex']['eqs']:
                                nicename += ' $\\bar t{}${}'.format('{:}', dset['tex']['eqs']['time'])
                            nicename += ' {}'.format(dset['tex']['eqs']['l'])
                            # nicename = nicename + ' $\Delta={:.3f}'.format(delt)
                            if 'num' in meta['legendcols']:
                                nicename = nicename + ' $n{}{}$'.format('{:}', ndata)

                            color = next(current_palette)
                            line = ax.hist(x=np.array(data), bins=np.array(edges), linewidth=1, histtype='step', density=False,
                                           label=nicename, color=color)
                            vmax = np.max([vmax, np.max(data)]) if vmax else np.max(data)
                            vmin = np.min([vmin, np.min(data)]) if vmin else np.min(data)
                            hmax = np.max([hmax, np.max(hist)]) if hmax else np.max(hist)

            if vwin:
                ax.axvline(vwin[0], color='grey', linestyle='dashed', linewidth=1.5)
                ax.axvline(vwin[1], color='grey', linestyle='dashed', linewidth=1.5,
                           label='Var$(H) \in [{:.2e}, {:.2e}]$'.format(vwin[0], vwin[1]))
            ax.set_xlim(left=np.min([np.max([vmin, 1e-18]), 1e-16]), right=np.max([vmax, 1e0]))
            ax.set_title('$\Delta = {:.3f}$'.format(delt))
            # ax.set_ylim(bottom=-0.01*hmax)
            # ax.set_yscale('log')
            # for win_idx, win in enumerate(variance_window_limits):
            #     for lim_idx, lim in enumerate(win):
            #         ax.axvline(lim, color='grey', linestyle='dashed', linewidth=1.5,
            #                    label=variance_window_names[win_idx][lim_idx])

            ax.set_xlabel('Var$(H)$')
            ax.set_ylabel('Histogram')
            ax.set_xscale('log')
            ax.legend(loc='upper right', framealpha=0.7, fontsize='x-small', labelspacing=0.25)
            fig.suptitle('$L = {}$'.format(size))
            fig.suptitle('Distribution of energy variance @ $L = {} $'.format(size))
            # $\Delta = $' + str(delt) + '$\lambda = $' + str(lamb))
            if not deltaidx in axes_used:
                axes_used.append(deltaidx)

        remove_empty_subplots(fig=fig, axes=axes, axes_used=axes_used)
        figname = "{}/{}_var_distribution_{}".format(meta['plotdir'], meta['plotprefix'], key0)
        plt.savefig('{}.pdf'.format(figname), format='pdf')
        plt.savefig('{}.png'.format(figname), format='png')
        plt.savefig('{}.svg'.format(figname), format='svg')

    return


def multiplot_var_distribution_old(src, plotdir='', algo_inc='', state_inc='', bins=200):
    print('Plotting: Energy Variance distribution for: ', algo_inc, state_inc)
    h5_src = h5open(src, 'r')
    path_L = h5py_unique_finder(h5_src, filter='L_', dep=1)
    path_l = h5py_unique_finder(h5_src, filter='l_', dep=2)
    path_d = h5py_unique_finder(h5_src, filter='d_', dep=3)
    # One figure per path_l, path_d and path_d
    for l in path_l:
        for d in path_d:
            # In each figure we want one subplot per path_l
            rows, cols = get_optimal_subplot_num(len(path_L))
            fig, axes = plt.subplots(nrows=rows, ncols=cols, figsize=(7 * cols, 7 * rows))
            fig.tight_layout(pad=5, w_pad=1.0, h_pad=1.0)
            fig.subplots_adjust(wspace=0.3, hspace=0.3)
            used_ax = 0
            delt = 0
            lamb = 0
            for ax, L in zip(np.ravel(axes), path_L):
                current_palette = itertools.cycle(sns.color_palette())
                basenode = h5_src[L][l][d]
                chain_length = basenode.attrs['model_size']
                delt = basenode.attrs['d']
                lamb = basenode.attrs['l']
                # h5keys = [item for item in h5_src[L][l][d].keys() if any(s in item for s in key_list)]
                # key_sorted = sorted(h5keys, key=natural_keys)
                for algokey, algopath, algonode in h5py_group_iterator(g=basenode, filter=algo_inc, dep=1):
                    for statekey, statepath, statenode in h5py_group_iterator(g=algonode, filter=state_inc, dep=1):
                        for win_idx, win in enumerate(variance_window_limits):
                            idx = get_v_filtered_index_list(statenode, win)
                            for datakey, datapath, datanode in h5py_node_finder(g=statenode, filter='energy_variance', dep=8):
                                if not idx:
                                    data = datanode['data']
                                else:
                                    data = datanode['data'][idx]
                                num = datanode['num'][()]
                                bins = np.logspace(start=-18, stop=0, num=64, endpoint=True)
                                hist, edges = np.histogram(data, bins=bins, density=False)
                                nicename = re.sub(r'[\W_]', ' ', str(algokey + " " + statekey))
                                nicename = nicename + ' (' + str(num) + ')'
                                color = next(current_palette)
                                ax.hist(x=np.array(data), bins=np.array(edges), linewidth=1, histtype='step', density=False,
                                        label=nicename, color=color)
                                ax.set_xlabel('$\sigma(H)^2$')
                                ax.set_ylabel('Histogram')
                                ax.set_xscale('log')
                                ax.set_title('$L = ' + str(chain_length) + '$')
                used_ax = used_ax + 1
                # color = next(current_palette)
                for win_idx, win in enumerate(variance_window_limits):
                    color = next(current_palette)
                    for lim_idx, lim in enumerate(win):
                        # ax.axvline(lim, color=variance_colors[win_idx], linestyle='dashed', linewidth=1.5, label='Avg ' + nicename)
                        ax.axvline(lim, color='grey', linestyle='dashed', linewidth=1.5,
                                   label=variance_window_names[win_idx][lim_idx])
                ax.legend(loc='upper right')
            fig.suptitle('Distribution of energy variance @ $\Delta = $' + str(delt) + '$\lambda = $' + str(lamb))
            for ax in np.ravel(axes)[used_ax:]:
                fig.delaxes(ax)
            if plotdir != '':
                plt.savefig(plotdir + '/Var_distribution_' + l + '_' + d + '.pdf', format='pdf')
                plt.savefig(plotdir + '/Var_distribution_' + l + '_' + d + '.png', format='png')

    h5close(h5_src)
