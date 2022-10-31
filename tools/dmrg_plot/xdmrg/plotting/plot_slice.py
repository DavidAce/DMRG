from matplotlib.ticker import MaxNLocator
from src.database.database import *
from src.plotting.tools import *
from src.general.filter import *
from itertools import product
from copy import deepcopy
import numpy as np
import seaborn as sns
import matplotlib.patheffects as pe
from pathlib import Path
import tikzplotlib


def plot_slice(h5, db, meta, sub1, l1, x1):
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
                        ykeys = [key0, key1, get_keys(db, x1[0]), algokey, statekey, pointkey, meta['ygrp']]
                        ygrps = get_matching_path(db, ykeys)

                        if not ygrps:
                            raise KeyError(
                                "Could not find ygrp [{}] in: {}".format(meta['ygrp'], ykeys))
                        ydata = get_matching_data(h5, ygrps, 'avg', meta['ycol'], 'f8')
                        edata = get_matching_data(h5, ygrps, 'ste', meta['ycol'], 'f8')
                        ndata = get_matching_data(h5, ygrps, 'num', meta['ycol'], 'f8')
                        xdata = get_matching_vals(h5, ygrps, db=db, keys=x1[0])
                        ldata = get_legend_rows(db=db, h5=h5, grps=ygrps, legend_col_keys=legend_col_keys)

                        sort = np.argsort(xdata)
                        xdata = np.asarray(xdata)[sort]
                        ydata = np.asarray(ydata)[sort]
                        edata = np.asarray(edata)[sort]
                        ndata = np.asarray(ndata)[sort]
                        ldata = np.asarray(ldata)[sort]

                        if any(len(d) == 0 for d in [xdata, ydata, edata, ndata, ldata]):
                            continue

                        if ydiv := meta.get('ydiv'):
                            ddata = get_matching_vals(h5, ygrps, db=db, keys=ydiv)
                            ddata = np.asarray(ddata)[sort]
                            ydata /= ddata

                        if "states" in statekey:
                            color = next(ed_palette)
                            mstyle = 'o'
                        else:
                            color = next(xd_palette)

                        line, _, _ = ax.errorbar(x=xdata, y=ydata, yerr=edata, color=color, marker=mstyle, linestyle=lstyle,
                                                 path_effects=path_effects)
                        for icol, (col, key) in enumerate(zip(ldata[-1], legend_col_keys)):  # Just take the last legend row
                            key, fmt = key.split(':') if ':' in key else [key, '']
                            f['legends'][idx][icol]['handle'].append(line)
                            f['legends'][idx][icol]['label'].append(col)
                            f['legends'][idx][icol]['title'] = db['tex'][key]

                        if not idx in f['axes_used']:
                            f['axes_used'].append(idx)

                        dbval = db['dsets'][ygrps[0]]
                        ax.set_title(get_title(dbval, sub1), x=0.5, horizontalalignment='left')
                        if ylabel := meta.get('ylabel'):
                            ax.set_ylabel(ylabel)
                        if xlabel := meta.get('xlabel'):
                            ax.set_ylabel(xlabel)
                        else:
                            ax.set_xlabel(get_tex(db, x1[0]))
                        if xlim := meta.get('xlim', {}).get(x1[0]):
                            ax.set_xlim(xmin=xlim[0], xmax=xlim[1], auto=True)
                        if ylim := meta.get('ylim'):
                            ax.set_ylim(ymin=ylim[0], ymax=ylim[1], auto=True)

    if titlename := meta.get('titlename'):
        f['fig'].suptitle('{} vs {}'.format(titlename, db['tex'][x1[0]]))

    prettify_plot5(fmeta=f)
    figname = "{}/{}_slice({})_sub({})_line({})".format(meta['plotdir'], meta['plotprefix'], x1[0], sub1[0], l1[0])
    f['fig'].savefig('{}.pdf'.format(figname), format='pdf')
    f['fig'].savefig('{}.png'.format(figname), format='png')
    f['fig'].savefig('{}.svg'.format(figname), format='svg')
    f['fig'].savefig('{}.pgf'.format(figname), format='pgf')
    tikzplotlib.save('{}.tex'.format(figname))
