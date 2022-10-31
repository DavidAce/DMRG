from src.plotting.tools import *
import matplotlib.pyplot as plt
from src.database.database import *
from src.general.filter import *
import itertools
import datetime


def multiplot_time_distribution(h5, db, meta):
    print('Plotting: Simulation time distribution')
    algoinc = meta.get('inc').get('algo')
    stateinc = meta.get('inc').get('state')
    pointinc = meta.get('inc').get('point')
    timecounter = 0
    simcounter = 0
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
                            num = dset['num']
                            style = dset['style']

                            if "states" in statekey:
                                continue
                            print(datapath)
                            data = h5[datapath]['data'][meta['colname']][()]
                            # avg = h5[datapath]['avg'][meta['colname']][()]
                            hist, edges = np.histogram(data / 60, bins=50, density=False)
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
                                nicename += ' {}'.format(dset['tex']['eqs']['time'])
                            nicename += ' {}'.format(dset['tex']['eqs']['l'])
                            # nicename = nicename + ' $\Delta={:.3f}'.format(delt)
                            if 'num' in meta['legendcols']:
                                nicename = nicename + ' $n{}{}$'.format('{:}', num)

                            bincentres = [(edges[i] + edges[i + 1]) / 2. for i in range(len(edges) - 1)]
                            widths = np.diff(edges)
                            color = next(current_palette)
                            lwidth = 1.25
                            lalpha = 1.0
                            hmax = np.max(hist)
                            ax.step(bincentres, hist, where='mid', label=nicename,
                                    color=color, alpha=lalpha, linewidth=lwidth)
                            # ax.axvline(avg / 60, color=color, linestyle='dashed', linewidth=lwidth)
                            timecounter = timecounter + np.sum(data)
                            simcounter = simcounter + num

            # ax.set_title('$\Delta = {:.3f}$'.format(delt))
            # ax.set_ylim(bottom=-0.01*hmax)
            # ax.set_yscale('log')
            # for win_idx, win in enumerate(variance_window_limits):
            #     for lim_idx, lim in enumerate(win):
            #         ax.axvline(lim, color='grey', linestyle='dashed', linewidth=1.5,
            #                    label=variance_window_names[win_idx][lim_idx])

            ax.set_xlabel('Time [min]')
            ax.set_ylabel('Histogram')
            ax.set_yscale('log')
            ax.set_xscale('log')
            if hmax:
                ax.set_ylim(ymin=0.1, ymax=1.1 * hmax)
            ax.legend(loc='upper right', framealpha=0.7, fontsize='x-small', labelspacing=0.25)
            ax.set_title('$\Delta = {:.4f}$'.format(delt))
            fig.suptitle('$L = {}$'.format(size))
            fig.suptitle('Distribution of Simulation Time @ $L = {} $'.format(size))
            # $\Delta = $' + str(delt) + '$\lambda = $' + str(lamb))
            if not deltaidx in axes_used:
                axes_used.append(deltaidx)

        remove_empty_subplots(fig=fig, axes=axes, axes_used=axes_used)
        figname = "{}/{}_distribution_{}".format(meta['plotdir'], meta['plotprefix'], key0)
        plt.savefig('{}.pdf'.format(figname), format='pdf')
        plt.savefig('{}.png'.format(figname), format='png')
        plt.savefig('{}.svg'.format(figname), format='svg')
    time_tot = str(datetime.timedelta(seconds=timecounter))
    time_sim = str(datetime.timedelta(seconds=timecounter / simcounter))
    print("Total sim time      : ", time_tot)
    print("Total num sims      : ", simcounter)
    print("Average time per sim: ", time_sim)

    return


def multiplot_time_distribution_old(h5, db, meta):
    print('Plotting: Time distribution')
    timecounter = 0
    simcounter = 0
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
                current_palette = itertools.cycle(sns.color_palette())
                key_num = 0
                basenode = h5_src[L][l][d]
                chain_length = basenode.attrs['model_size']
                delt = basenode.attrs['delta']
                lamb = basenode.attrs['lambda']
                for algokey, algopath, algonode in h5py_group_iterator(node=basenode, keypattern=algo_inc, dep=1):
                    for statekey, statepath, statenode in h5py_group_iterator(node=algonode, keypattern=state_inc, dep=1):
                        for win_idx, win in enumerate(variance_window_limits):
                            idx = get_v_filtered_index_list(statenode, win)
                            for datakey, datapath, datanode in h5py_node_finder(node=statenode,
                                                                                keypattern=['measurements', 'algorithm_time'],
                                                                                dep=8):
                                if "checkpoint" in datapath:
                                    continue
                                if 'states' in statekey and 'algorithm_time' in datakey:
                                    data = np.asarray(datanode['data'])
                                else:
                                    data = np.array(datanode['data']['algorithm_time'])[idx]
                                if np.any(np.isnan(data)):
                                    raise ValueError("Data contains nan's")
                                key_num = key_num + 1
                                num = np.size(data)
                                avg = np.mean(data)
                                hist, edges = np.histogram(data / 60, bins=20, density=False)
                                nicename = re.sub(r'[\W_]', ' ', str(algokey + " " + statekey))
                                nicename = nicename + ' (' + str(num) + ')'
                                bincentres = [(edges[i] + edges[i + 1]) / 2. for i in range(len(edges) - 1)]
                                widths = np.diff(edges)
                                color = next(current_palette)
                                lwidth = 1.25
                                lalpha = 1.0
                                ax.step(bincentres, hist, where='mid', label=nicename,
                                        color=color, alpha=lalpha, linewidth=lwidth)
                                ax.axvline(avg / 60, color=color, linestyle='dashed', linewidth=lwidth)
                                timecounter = timecounter + np.sum(data)
                                simcounter = simcounter + num
                ax.set_xlabel('Time [minutes]')
                ax.set_ylabel('Histogram')
                ax.set_title('$L = ' + str(chain_length) + '$')
                ax.set_xlim(left=0)
                ax.set_yscale('log')
                ax.legend()
                used_ax = used_ax + 1
            fig.suptitle(
                'Distribution Simulation time @ $\Delta = ' + str(delt) + '\quad \lambda = ' + str(lamb) + '$')
            for ax in np.ravel(axes)[used_ax:]:
                fig.delaxes(ax)

            if plotdir != '':
                plt.savefig(plotdir + '/Time_distribution_' + l + '_' + d + '.pdf', format='pdf')
                plt.savefig(plotdir + '/Time_distribution_' + l + '_' + d + '.png', format='png')

    time_tot = str(datetime.timedelta(seconds=timecounter))
    time_sim = str(datetime.timedelta(seconds=timecounter / simcounter))
    print("Total sim time      : ", time_tot)
    print("Total num sims      : ", simcounter)
    print("Average time per sim: ", time_sim)

    h5close(h5_src)
