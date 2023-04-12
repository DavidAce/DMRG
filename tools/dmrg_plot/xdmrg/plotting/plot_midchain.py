from matplotlib.ticker import MaxNLocator
# from src.database.database import *
from .tools import *
from src.general.filter import *
from itertools import product
from copy import deepcopy
import numpy as np
import seaborn as sns
import matplotlib.patheffects as pe
from pathlib import Path
import tikzplotlib


def plot_midchain(h5, db, meta, sub1, l1, x1):
    if len(sub1) != 1:
        raise AssertionError("sub must have length 1")
    if len(l1) != 1:
        raise AssertionError("l must have length 1")
    if len(x1) != 1:
        raise AssertionError("x must have length 1")

    algoinc = meta.get('inc').get('algo')
    stateinc = meta.get('inc').get('state')
    pointinc = meta.get('inc').get('point')

    if 'mplstyle' in meta:
        plt.style.use(meta['mplstyle'])
    if 'plotdir' in meta and 'mplstyle' in meta:
        if Path(meta['plotdir']).stem != Path(meta['mplstyle']).stem:
            meta['plotdir'] = Path(meta['plotdir'], Path(meta['mplstyle']).stem)
            Path(meta['plotdir']).mkdir(parents=True, exist_ok=True)
            print("Setting plotdir: ", meta['plotdir'])
    if 'mplstyle' in meta and 'slack' in meta['mplstyle']:
        palette_name = "husl"
        path_effects = [pe.SimpleLineShadow(offset=(0.6, -0.6), alpha=0.35), pe.Normal()]
    else:
        palette_name = "husl"
        path_effects = None

    prb_style = 'prb' in meta['mplstyle'] if 'mplstyle' in meta else False
    print(l1)

    legend_col_keys = l1.copy()
    if legendcols := meta['legendcols']:
        for col in legendcols:
            if not col.split(':')[0] in [l.split(':')[0] for l in sub1 + x1 + l1]:
                print('appending col: {}'.format(col))
                legend_col_keys.append(col)

    numplots = len(get_keys(db, sub1[0]))

    f = get_fig_meta(numplots, meta=meta)
    for idx, ((key0), ax) in enumerate(zip(get_keys(db, sub1[0]), f['ax'])):
        lstyles = get_linestyles()
        for algokey in iterate(db['keys']['algo'], algoinc):
            for statekey in iterate(db['keys']['state'], stateinc):
                lstyle = next(lstyles)
                mstyles = get_markerstyles()
                dbval = None
                for pointkey in iterate(db['keys']['point'], pointinc):
                    ed_palette = itertools.cycle(sns.color_palette("Set2"))
                    xd_palette = itertools.cycle(sns.color_palette(palette_name, len(get_keys(db, l1[0]))))
                    mstyle = next(mstyles)
                    for key1 in get_keys(db, l1[0]):
                        xdata, ydata, edata, ndata, ldata = [], [], [], [], []
                        for key2 in get_keys(db, x1[0]):
                            findlist = [key0, key1, key2, algokey, statekey, pointkey, meta['ygrp']]
                            datapaths = [value['path']['data'] for key, value in db['dsets'].items() if
                                         all(findlist[-1] in key and '/' + k + '/' in key for k in findlist[:-1])]
                            # print("found", len(datapaths), "datapaths -- findlist: ", findlist)
                            if not datapaths:
                                continue
                            print("found", len(datapaths), "data paths -- findlist: ", findlist)
                            for datapath in datapaths:
                                print('-- {}'.format(datapath))
                            x, y, e, n, l = [], [], [], [], []
                            for datapath in datapaths:
                                datanode = h5[datapath]
                                dbval = db['dsets'][datanode.name]
                                midx = dbval['vals']['midx']
                                y.append(get_data(datanode['avg'], 'L_', 'f8')[midx])
                                e.append(get_data(datanode['ste'], 'L_', 'f8')[midx])
                                n.extend(get_data(datanode['num'], 'L_', 'f8')[()])
                                x.append(get_vals(dbval, x1[0]))
                                l.append(get_legend_row(db=db, datanode=datanode, legend_col_keys=legend_col_keys))
                            if not x or not y or not e or not n:
                                continue
                            sort = np.argsort(x)
                            xdata.append(np.asarray(x)[sort])
                            ydata.append(np.asarray(y)[sort])
                            edata.append(np.asarray(e)[sort])
                            ndata.append(np.asarray(n)[sort])
                            ldata.append(np.asarray(l)[sort])

                        if not xdata or not ydata or not edata or not ndata:
                            continue
                        if "states" in statekey:
                            color = next(ed_palette)
                            mstyle = 'o'
                        else:
                            color = next(xd_palette)

                        for x, y, e, n, l in zip(np.asarray(xdata).T, np.asarray(ydata).T, np.asarray(edata).T, np.asarray(ndata).T, np.asarray(ldata).T):
                            line, _, _ = ax.errorbar(x=x, y=y, yerr=e, color=color, marker=mstyle, linestyle=lstyle,
                                                     path_effects=path_effects)
                            # A plot line needs a single legend row, so we need to reduce l down to just 1 row.
                            # The simplest is to take the last one, but in general we could average over them.
                            for icol, (col, key) in enumerate(zip(l[-1], legend_col_keys)):
                                key, fmt = key.split(':') if ':' in key else [key, '']
                                f['legends'][idx][icol]['handle'].append(line)
                                f['legends'][idx][icol]['label'].append(col)
                                f['legends'][idx][icol]['title'] = db['tex'][key]

                        f['axes_used'].append(idx)
                        ax.set_xlabel(get_tex(db, x1[0]))
                        if ylabel := meta.get('ylabel'):
                            rangle2 = '\\rangle\\rangle'
                            if rangle2 in ylabel:
                                ylabel = ylabel.replace(rangle2, '(L/2)' + rangle2)
                                ax.set_ylabel(ylabel)
                            else:
                                ax.set_ylabel(ylabel + '$(L/2)$')
                        if xlim := meta.get('xlim', {}).get(x1[0]):
                            ax.set_xlim(xmin=xlim[0], xmax=xlim[1], auto=True)
                        if ylim := meta.get('ylim'):
                            ax.set_ylim(ymin=ylim[0], ymax=ylim[1], auto=True)
                        if dbval:
                            ax.set_title(get_title(dbval, sub1), x=0.5, horizontalalignment='left')
    if titlename := meta.get('titlename'):
        f['fig'].suptitle('{} vs {}'.format(titlename, db['tex'][x1[0]]))

    prettify_plot5(fmeta=f)
    figname = "{}/{}_midchain({})_sub({})_line({})".format(meta['plotdir'], meta['plotprefix'], x1[0], sub1[0], l1[0])
    f['fig'].savefig('{}.pdf'.format(figname), format='pdf')
    f['fig'].savefig('{}.png'.format(figname), format='png')
    f['fig'].savefig('{}.svg'.format(figname), format='svg')
    f['fig'].savefig('{}.pgf'.format(figname), format='pgf')
    tikzplotlib.save('{}.tex'.format(figname))


def multiplot_Smid_sub1_l1_x1(h5_avg, db, meta, plotdir='', g3=['algo', 'state', 'point'], algo_inc='', state_inc='',
                              sub1=['L'], l1=['l'], x1=['d'], vwin=None):
    print('Plotting: S vs Site for: ', algo_inc, state_inc)
    if len(sub1) != 1:
        raise AssertionError("sub2 must have length 1")
    if len(l1) != 1:
        raise AssertionError("l1 must have length 1")
    if len(x1) != 1:
        raise AssertionError("x1 must have length 1")
    # One figure per L
    # One subplot per delta
    # One line per lambda per state
    # One figure per L
    # One subplot per delta
    # One line per lambda per state
    numsub = len(db['keys'][sub1[0]])
    figrows, figcols = get_optimal_subplot_num(numsub)
    fig, axes = plt.subplots(nrows=figrows, ncols=figcols, figsize=(5 * figcols, 5 * figrows), sharey='all')
    fig.tight_layout(pad=5, w_pad=1.0, h_pad=1.0)
    fig.subplots_adjust(wspace=0.2, hspace=0.2)
    axes_used = []
    for subidx, (key0, ax) in enumerate(zip(db['keys'][sub1[0]], np.ravel(axes))):
        ed_palette = itertools.cycle(sns.color_palette("Set2"))
        numcolors = len(db['keys'][g3[0]]) * len(db['keys'][g3[1]]) * len(db['keys'][g3[2]])
        current_palette = itertools.cycle(sns.color_palette("colorblind", numcolors))
        lstyles = itertools.cycle(['-.', '-', '--', ':', ])
        mstyles = itertools.cycle(('.', ',', '+', 'o', '*'))
        delt = None
        size = None
        for lidx, key1 in enumerate(db['keys'][l1[0]]):
            for gidx, (key2, key3, key4) in enumerate(product(db['keys'][g3[0]], db['keys'][g3[1]], db['keys'][g3[2]])):
                xdata, ydata, edata, ndata = [], [], [], []
                for key5 in db['keys'][x1[0]]:
                    findlist = [key0, key1, key2, key3, key4, key5, meta['dsetname']]
                    datanode = [value['datanode'] for key, value in db['dsets'].items() if
                                all(k in key for k in findlist)]
                    if len(datanode) != 1:
                        continue
                    if not any(k in datanode[0].name for k in algo_inc):
                        continue
                    if not any(k in datanode[0].name for k in state_inc):
                        continue

                        # raise LookupError("Found incorrect number of datanodes")
                    datanode = datanode[0]
                    dbval = db['dsets'][datanode.name]
                    midx = dbval['midx']
                    lamb = dbval['l']
                    delt = dbval['d']
                    size = dbval['L']
                    ydata.append(get_data(datanode['avg'], 'L_', 'f8')[midx])
                    edata.append(get_data(datanode['ste'], 'L_', 'f8')[midx])
                    ndata.append(get_data(datanode['num'], 'L_', 'f8')[()])
                    xdata.append(dbval[x1[0]])
                    print("ndata: {}".format(ndata))
                if not xdata or not ydata:
                    continue
                sort = np.argsort(xdata)
                ydata = np.array(ydata)[sort]
                edata = np.array(edata)[sort]
                xdata = np.array(xdata)[sort]
                print("ydata after:", np.shape(ydata))
                print("edata after:", np.shape(edata))
                print("xdata after:", np.shape(xdata))
                print("ndata after:", np.shape(ndata))
                ndata = np.min(ndata)

                statekey = dbval['keys']['state']
                algokey = dbval['keys']['algo']

                style = dbval['style']

                lstyle = next(lstyles)
                mstyle = next(mstyles)

                if vwin:
                    ydata, xdata, edata, ndata = get_v_filtered_edata(ydata, xdata, edata, ndata, dbval, vwin)

                if "states" in statekey:
                    color = next(ed_palette)
                    mstyle = None
                    lstyle = 'solid'
                    lgnd = "ED ({})".format(ndata)
                    # nicename = "ED e=[" + statenode.attrs["efmt"] + "]"
                else:
                    color = next(current_palette)
                    lgnd = "{} {} ({})".format(re.sub(r'[\W_]', ' ', "{} {}".format(algokey, statekey)),
                                               db['tex'][key1], ndata)
                print("Plotting xdata {}, legend {}".format(np.shape(xdata), np.shape(lgnd)))
                ax.errorbar(x=xdata, y=ydata, yerr=edata, label=lgnd, capsize=2,
                            color=color, elinewidth=0.3, markeredgewidth=0.8,
                            marker=mstyle, markersize=4.0, linestyle=lstyle,
                            linewidth=style['lwidth'], alpha=style['lalpha'])

        ax.set_xlabel(db['tex'][x1[0]])
        ax.set_ylabel(meta['ylabel'] + "(L/2)")
        ax.set_title('{}'.format(db['tex'][key0]))
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        ax.legend(loc='lower center', fontsize='small', labelspacing=0.25)
        fig.suptitle('Entanglement entropy vs {} @ {}'.format(db['tex'][x1[0]], db['tex'][key0]))

        # $\Delta = $' + str(delt) + '$\lambda = $' + str(lamb))
        axes_used.append(subidx) if not subidx in axes_used else axes_used

        remove_empty_subplots(fig=fig, axes=axes, axes_used=axes_used)
        if plotdir != '':
            plt.savefig('{}/Smid_vs_Delta_{}.pdf'.format(plotdir, key0), format='pdf')
            plt.savefig('{}/Smid_vs_Delta_{}.png'.format(plotdir, key0), format='png')

    return

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
            max_S = 0
            for ax, L in zip(np.ravel(axes), path_L):
                ed_palette = itertools.cycle(sns.color_palette("Set2"))
                current_palette = itertools.cycle(sns.color_palette())
                basenode = h5_src[L][l][d]
                size = basenode.attrs['L']
                delt = basenode.attrs['d']
                lamb = basenode.attrs['l']
                for algokey, algopath, algonode in h5py_group_iterator(g=basenode, filter=algo_inc, dep=1):
                    for statekey, statepath, statenode in h5py_group_iterator(g=algonode, filter=state_inc,
                                                                              dep=1):
                        for v_win_idx, v_win in enumerate(variance_window_limits):
                            for e_win_idx, e_win in enumerate(energy_window_limits):
                                idx = get_v_e_filtered_index_list(statenode, v_win, e_win)
                                for datakey, datapath, datanode in h5py_node_finder(g=statenode,
                                                                                    filter='entanglement_entropies',
                                                                                    dep=8):
                                    ndata = datanode['num'][()]
                                    ydata = np.array(datanode['avg'])
                                    edata = np.array(datanode['ste'])
                                    xdata = range(len(ydata))

                                    if np.any(np.isnan(ydata)):
                                        raise ValueError("Data contains nan's")
                                    if np.any(np.isnan(ydata)):
                                        raise ValueError("Standard error contains nan's")

                                    if "states" in statekey:
                                        color = next(ed_palette)
                                        nicename = "ED e=[" + statenode.attrs["efmt"] + "]"
                                        lwidth = 3.0
                                        lalpha = 0.9
                                        mstyle = None
                                        lstyle = 'solid'
                                        # Test taking a subset of states
                                        # edata = bootstrap_sterr(data=datanode['data'],chunksize=2000)

                                        # idx = random.sample(range(0,ndata-1),2000)
                                        # idx.sort()
                                        # ydata = np.nanmean(np.array(datanode['data'][:,idx]),axis=1)


                                    else:
                                        color = next(current_palette)
                                        nicename = re.sub(r'[\W_]', ' ', str(algokey + " " + statekey))
                                        lwidth = 1.4
                                        lalpha = 1.0
                                        mstyle = '.'
                                        lstyle = 'dotted'

                                    nicename = nicename + ' (' + str(ndata) + ')'
                                    ax.errorbar(x=xdata, y=ydata, yerr=edata, label=nicename, capsize=2,
                                                color=color,
                                                elinewidth=0.3, markeredgewidth=0.8, marker=mstyle,
                                                linestyle=lstyle,
                                                linewidth=lwidth, alpha=lalpha)
                                    max_S = np.max([max_S, np.max(ydata)])

                                    # Select a random subset of the data
                                    if "states" in statekey:
                                        continue
                                    data = np.array(datanode['data'])
                                    ndata = int(data.shape[1] * 0.5)
                                    for rep in range(0):
                                        idx = np.random.choice(data.shape[1], ndata, replace=False)
                                        subdata = data[:, idx]
                                        std = np.nanstd(subdata, axis=1)
                                        ydata = np.nanmean(subdata, axis=1)
                                        edata = std / np.sqrt(ndata)
                                        nicename = re.sub(r'[\W_]', ' ',
                                                          str(algokey + " " + statekey)) + ' (random ' + str(
                                            ndata) + ')'
                                        ax.errorbar(x=xdata, y=ydata, yerr=edata, label=nicename, capsize=2,
                                                    elinewidth=0.3, markeredgewidth=0.8, linewidth=0.3)
                                        max_S = np.max([max_S, np.max(ydata)])

                ax.set_xlabel('Site $l$')
                ax.set_ylabel('$S_E(l)$')
                ax.set_title('$L = ' + str(size) + '$')
                used_ax = used_ax + 1
                ax.legend()

            for ax in np.ravel(axes)[used_ax:]:
                fig.delaxes(ax)
            for ax in np.ravel(axes):
                ax.set_ylim(ymin=0, ymax=max_S * 1.2)
            fig.suptitle('Entanglement entropy vs site @ $\Delta = ' + str(delt) + '\quad \lambda = ' + str(
                lamb) + '$')
            if plotdir != '':
                plt.savefig(plotdir + '/S_vs_Site_' + l + '_' + d + '.pdf', format='pdf')
                plt.savefig(plotdir + '/S_vs_Site_' + l + '_' + d + '.png', format='png')

    h5close(h5_src)
