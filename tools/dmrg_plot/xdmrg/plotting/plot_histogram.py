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


def plot_histogram(h5, db, meta, sub1, l1):
    if len(sub1) != 1:
        raise AssertionError("sub must have length 1")
    if len(l1) != 1:
        raise AssertionError("l must have length 1")

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

    legend_col_keys = l1.copy()
    if legendcols := meta['legendcols']:
        for col in legendcols:
            if not col.split(':')[0] in [l.split(':')[0] for l in sub1 + l1]:
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
                        findlist = [key0, key1, algokey, statekey, pointkey]
                        ygrps = [value['path']['data'] for key, value in db['dsets'].items() if
                                 all(meta['ygrp'] in key and '/{}/'.format(k) in key for k in findlist)]

                        if not ygrps:
                            raise KeyError("Could not find ygrps: \n findlist: {}".format(findlist))
                        for ygrp in ygrps:
                            if not 'data' in h5[ygrp]:
                                raise KeyError('Could not find dataset [data] in {}'.format(ygrp))
                            dbval = db['dsets'][ygrp]
                            data = get_data(h5[ygrp]['data'], meta['ycol'], 'f8')
                            ldata = get_legend_row(db=db, datanode=h5[ygrp], legend_col_keys=legend_col_keys)

                            if "states" in statekey:
                                color = next(ed_palette)
                                mstyle = 'o'
                            else:
                                color = next(xd_palette)

                            hist, edges = np.histogram(data, bins=meta['bins'], density=False)
                            line = ax.step(x=edges[:-1], y=hist, where='mid', color=color, linestyle=lstyle,
                                           path_effects=path_effects)
                            for icol, (col, key) in enumerate(zip(ldata, legend_col_keys)):
                                key, fmt = key.split(':') if ':' in key else [key, '']
                                f['legends'][idx][icol]['handle'].append(line[0])
                                f['legends'][idx][icol]['label'].append(col)
                                f['legends'][idx][icol]['title'] = db['tex'][key]

                            f['axes_used'].append(idx)
                            if ylabel := meta.get('ylabel'):
                                ax.set_ylabel(ylabel)
                            if xlabel := meta.get('xlabel'):
                                ax.set_xlabel(xlabel)
                            if xlim := meta.get('xlim', {}):
                                ax.set_xlim(xmin=xlim[0], xmax=xlim[1], auto=True)
                            if ylim := meta.get('ylim'):
                                ax.set_ylim(ymin=ylim[0], ymax=ylim[1], auto=True)
                            if dbval:
                                ax.set_title(get_title(dbval, sub1), x=0.5, horizontalalignment='left')

    if titlename := meta.get('titlename'):
        f['fig'].suptitle(titlename)

    prettify_plot5(fmeta=f)
    figname = "{}/{}_histogram({})_sub({})_line({})".format(meta['plotdir'], meta['plotprefix'], meta['ycol'], sub1[0], l1[0])
    f['fig'].savefig('{}.pdf'.format(figname), format='pdf')
    f['fig'].savefig('{}.png'.format(figname), format='png')
    f['fig'].savefig('{}.svg'.format(figname), format='svg')
    f['fig'].savefig('{}.pgf'.format(figname), format='pgf')
    tikzplotlib.save('{}.tex'.format(figname))
