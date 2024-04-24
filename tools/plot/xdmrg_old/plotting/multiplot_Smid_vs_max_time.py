from .tools import *
import matplotlib.pyplot as plt
# from src.database.database import *
from src.general.filter import *
import itertools


def multiplot_Smid_vs_tmax_v1(h5_src, db=None, plotdir='', algo_inc=None, state_inc=None):
    print('Plotting: S(L/2) vs max time for : ', algo_inc, state_inc)
    if not db:
        db = load_database(h5_src, 'entanglement_entropy_midchain', algo_inc, state_inc)
    # One figure per L
    # One subplot per delta
    # One line per lambda per state
    for Lkey in db['keys']['L']:
        figrows, figcols = get_optimal_subplot_num(1 + len(db['keys']['d']))  # Add one for the legend
        fig, axes = plt.subplots(nrows=figrows, ncols=figcols, figsize=(5 * figcols, 5 * figrows), sharey='all',
                                 sharex='all')
        fig.tight_layout(pad=5, w_pad=1.0, h_pad=1.0)
        fig.subplots_adjust(wspace=0.2, hspace=0.2)
        axes_used = []
        max_S = None
        min_S = None
        for didx, (dkey, ax) in enumerate(zip(db['keys']['d'], np.ravel(axes))):
            ed_palette = itertools.cycle(sns.color_palette("Set2"))
            numcolors = len(db['keys']['state']) * len(db['keys']['l'])
            current_palette = itertools.cycle(sns.color_palette("colorblind", numcolors))
            lstyles = itertools.cycle(['-.', '-', '--', ':', ])
            mstyles = itertools.cycle(('.', ',', '+', 'o', '*'))
            delt = None
            size = None
            nlabel = {'line': [], 'text': []}
            for algoidx, algokey in enumerate(db['keys']['algo']):
                for stateidx, statekey in enumerate(db['keys']['state']):
                    if not contains(algokey, algo_inc) or not contains(statekey, state_inc):
                        continue
                    dsetkeys = [x for x in db['dsets'] if
                                Lkey in x and dkey in x and algokey in x and statekey in x]
                    if not dsetkeys:
                        continue
                    lstyle = next(lstyles)
                    mstyle = next(mstyles)
                    # Now we have a set of setkeys with fixed L,d,algo and state, varying l.
                    for idx, dsetkey in enumerate(dsetkeys):
                        meta = db['dsets'][dsetkey]
                        midx = meta['midx']
                        lamb = meta['l']
                        delt = meta['d']
                        size = meta['L']
                        ndata = meta['num']
                        style = meta['style']
                        if "states" in statekey:
                            ydata = meta['datanode']['avg'][midx]
                            lgnd = "ED $\lambda = {:.3f}$ ({})".format(lamb, ndata)
                            color = next(ed_palette)
                            ax.axhline(ydata, color=color, linestyle='dashed', linewidth=1, label=lgnd)
                        else:
                            xdata, tdata, ydata, edata = [], [], [], []
                            color = next(current_palette)
                            lgnd = "{} $\lambda = {:.3f}$".format(re.sub(r'[\W_]', ' ', str(algokey + " " + statekey)), lamb)
                            totnum = meta['datanode']['num'][()]
                            oldnum = 0
                            t_inf = get_time_maximum(meta['pointnode'])
                            for t_max in np.linspace(0.1, t_inf, 10):
                                idx = get_time_filtered_index_list(meta['pointnode'], [0, t_max])
                                newnum = len(idx[0])
                                if newnum == oldnum:
                                    continue
                                if (len(idx[0]) == 0):
                                    tdata.append(t_max)
                                    ydata.append(0)
                                    edata.append(0)
                                    continue
                                oldnum = newnum
                                data = meta['datanode']['data'][midx][idx[0]]
                                axis = 0
                                num = np.shape(data)[axis]
                                tdata.append(t_max)
                                ydata.append(np.nanmean(data, axis=axis))
                                xdata.append(newnum / totnum * 100)
                                edata.append(np.nanstd(data, axis=axis) / np.sqrt(num))

                            if np.any(np.isnan(ydata)):
                                raise ValueError("Data contains nan's")
                            if np.any(np.isnan(ydata)):
                                raise ValueError("Standard error contains nan's")
                            tdata = np.asarray(tdata)
                            ydata = np.asarray(ydata)
                            edata = np.asarray(edata)
                            half = int(len(tdata) / 2)
                            qrtr = int(len(tdata) / 4)
                            max_S = np.max([max_S, np.max(ydata[half:])]) if max_S else np.max(ydata[half:])
                            min_S = np.min([min_S, np.min(ydata[qrtr:])]) if min_S else np.min(ydata[qrtr:])
                            # Note the comma below! It's important!
                            line, = ax.plot(tdata / 60, ydata, label=lgnd, color=color,
                                            # capsize=1, elinewidth=0.3, markeredgewidth=0.8,
                                            linestyle=lstyle, marker=mstyle, markersize=5.0,  # style['mstyle'],
                                            linewidth=style['lwidth'], alpha=style['lalpha'])
                            ax.fill_between(tdata / 60, ydata - edata, ydata + edata, color=color, alpha=0.1, antialiased=True)
                            nlabel['line'].append(line)
                            nlabel['text'].append('{}'.format(ndata))
            ax.set_title('$\Delta = {:.4f}$'.format(delt))
            ax.set_xlabel('Simulation time [minutes]')
            ax.set_ylabel('$S_E(L/2)$')
            ax.legend(nlabel['line'], nlabel['text'], title='Realizations',
                      loc='lower right', framealpha=0.2, fontsize='small', labelspacing=0.25, ncol=2)
            fig.suptitle('$L = {}$'.format(size))
            axes_used.append(didx) if not didx in axes_used else axes_used

        prettify_plot(fig, axes, rows=figrows, cols=figcols, axes_used=axes_used, ymin=0.8 * min_S, ymax=1.2 * max_S)
        if plotdir != '':
            plt.savefig(plotdir + '/Smid_vs_tmax_' + Lkey + '.pdf', format='pdf')
            plt.savefig(plotdir + '/Smid_vs_tmax_' + Lkey + '.png', format='png')


def multiplot_Smid_vs_tmax_v2(h5_src, db=None, plotdir='', algo_inc=None, state_inc=None):
    print('Plotting: S(L/2) vs max time for : ', algo_inc, state_inc)
    if not db:
        db = load_database(h5_src, 'entanglement_entropy_midchain', algo_inc, state_inc)
    # One figure per L
    # One subplot per delta
    # One line per lambda per state
    for Lkey in db['keys']['L']:
        figrows, figcols = get_optimal_subplot_num(1 + len(db['keys']['d']))  # Add one for the legend
        fig, axes = plt.subplots(nrows=figrows, ncols=figcols, figsize=(5 * figcols, 5 * figrows), sharey='all',
                                 sharex='all')
        fig.tight_layout(pad=5, w_pad=1.0, h_pad=1.0)
        fig.subplots_adjust(wspace=0.2, hspace=0.2)
        axes_used = []
        max_S = None
        min_S = None
        for didx, (dkey, ax) in enumerate(zip(db['keys']['d'], np.ravel(axes))):
            ed_palette = itertools.cycle(sns.color_palette("Set2"))
            numcolors = len(db['keys']['state']) * len(db['keys']['l'])
            current_palette = itertools.cycle(sns.color_palette("colorblind", numcolors))
            lstyles = itertools.cycle(['-.', '-', '--', ':', ])
            mstyles = itertools.cycle(('.', ',', '+', 'o', '*'))
            delt = None
            size = None
            nlabel = {'line': [], 'text': []}
            for algoidx, algokey in enumerate(db['keys']['algo']):
                for stateidx, statekey in enumerate(db['keys']['state']):
                    if not contains(algokey, algo_inc) or not contains(statekey, state_inc):
                        continue
                    dsetkeys = [x for x in db['dsets'] if
                                Lkey in x and dkey in x and algokey in x and statekey in x]
                    if not dsetkeys:
                        continue
                    lstyle = next(lstyles)
                    mstyle = next(mstyles)
                    # Now we have a set of setkeys with fixed L,d,algo and state, varying l.
                    for idx, dsetkey in enumerate(dsetkeys):
                        meta = db['dsets'][dsetkey]
                        midx = meta['midx']
                        lamb = meta['l']
                        delt = meta['d']
                        size = meta['L']
                        ndata = meta['num']
                        style = meta['style']
                        if "states" in statekey:
                            ydata = meta['datanode']['avg'][midx]
                            lgnd = "ED $\lambda = {:.3f}$ ({})".format(lamb, ndata)
                            color = next(ed_palette)
                            ax.axhline(ydata, color=color, linestyle='dashed', linewidth=1, label=lgnd)
                        else:
                            xdata, tdata, ydata, edata = [], [], [], []
                            color = next(current_palette)
                            lgnd = "{} $\lambda = {:.3f}$".format(re.sub(r'[\W_]', ' ', str(algokey + " " + statekey)), lamb)
                            totnum = meta['datanode']['num'][()]
                            oldnum = 0
                            t_inf = get_time_maximum(meta['pointnode'])
                            for t_max in np.linspace(0.1, t_inf, 10):
                                idx = get_time_filtered_index_list(meta['pointnode'], [0, t_max])
                                newnum = len(idx[0])
                                if newnum == oldnum:
                                    continue
                                if (len(idx[0]) == 0):
                                    tdata.append(t_max)
                                    ydata.append(0)
                                    edata.append(0)
                                    continue
                                oldnum = newnum
                                data = meta['mmntnode']['data']['entanglement_entropy_midchain'][idx[0]]
                                axis = 0
                                num = np.shape(data)[axis]
                                tdata.append(t_max)
                                ydata.append(np.nanmean(data, axis=axis))
                                xdata.append(newnum / totnum * 100)
                                edata.append(np.nanstd(data, axis=axis) / np.sqrt(num))

                            if np.any(np.isnan(ydata)):
                                raise ValueError("Data contains nan's")
                            if np.any(np.isnan(ydata)):
                                raise ValueError("Standard error contains nan's")
                            tdata = np.asarray(tdata)
                            ydata = np.asarray(ydata)
                            edata = np.asarray(edata)
                            half = int(len(tdata) / 2)
                            qrtr = int(len(tdata) / 4)
                            max_S = np.max([max_S, np.max(ydata[half:])]) if max_S else np.max(ydata[half:])
                            min_S = np.min([min_S, np.min(ydata[qrtr:])]) if min_S else np.min(ydata[qrtr:])
                            # Note the comma below! It's important!
                            line, = ax.plot(tdata / 60, ydata, label=lgnd, color=color,
                                            # capsize=1, elinewidth=0.3, markeredgewidth=0.8,
                                            linestyle=lstyle, marker=mstyle, markersize=5.0,  # style['mstyle'],
                                            linewidth=style['lwidth'], alpha=style['lalpha'])
                            ax.fill_between(tdata / 60, ydata - edata, ydata + edata, color=color, alpha=0.1, antialiased=True)
                            nlabel['line'].append(line)
                            nlabel['text'].append('{}'.format(ndata))
            ax.set_title('$\Delta = {:.4f}$'.format(delt))
            ax.set_xlabel('Simulation time [minutes]')
            ax.set_ylabel('$S_E(L/2)$')
            ax.legend(nlabel['line'], nlabel['text'], title='Realizations',
                      loc='lower right', framealpha=0.2, fontsize='small', labelspacing=0.25, ncol=2)
            fig.suptitle('$L = {}$'.format(size))
            axes_used.append(didx) if not didx in axes_used else axes_used

        prettify_plot(fig, axes, rows=figrows, cols=figcols, axes_used=axes_used, ymin=0.8 * min_S, ymax=1.2 * max_S)
        if plotdir != '':
            plt.savefig(plotdir + '/Smid_vs_tmax_' + Lkey + '.pdf', format='pdf')
            plt.savefig(plotdir + '/Smid_vs_tmax_' + Lkey + '.png', format='png')
