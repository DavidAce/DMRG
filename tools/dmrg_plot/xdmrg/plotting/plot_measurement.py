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


def plot_measurement(h5, db, meta, sub1, l1, x1):
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

    legend_col_keys = l1.copy()
    if legendcols := meta['legendcols']:
        for col in legendcols:
            if not col in [l.split(':')[0] for l in sub1 + x1 + l1]:
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
                            findlist = [key0, key1, key2, algokey, statekey, pointkey, meta['load']]
                            datapaths = [value['path']['data'] for key, value in db['dsets'].items() if
                                         all(findlist[-1] in key and '/' + k + '/' in key for k in findlist[:-1])]
                            # print("found", len(datapaths), "datapaths -- findlist: ", findlist)
                            if not datapaths:
                                continue
                            for datapath in datapaths:
                                print('-- {}'.format(datapath))
                            x, y, e, n, l = [], [], [], [], []
                            for datapath in datapaths:
                                datanode = h5[datapath]
                                dbval = db['dsets'][datanode.name]
                                stats = meta['statname'] if 'statname' in meta else 'avg'
                                y.extend(get_data(datanode[stats], meta['colname'], 'f8')[()])
                                e.extend(get_data(datanode['ste'], meta['colname'], 'f8')[()])
                                n.extend(get_data(datanode['num'], meta['colname'], 'f8')[()])
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
                        for x, y, e, n, l in zip(zip(*xdata), zip(*ydata), zip(*edata), zip(*ndata), zip(*ldata)):
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
                        if xlabel := db.get('texdef').get(x1[0]) or db.get('tex').get(x1[0]):
                            ax.set_xlabel(xlabel)
                        if ylabel := meta.get('ylabel'):
                            ax.set_ylabel(ylabel)
                        if dbval:
                            ax.set_title(get_title(dbval, sub1), x=0.5, horizontalalignment='left')
    if titlename := meta.get('titlename'):
        f['fig'].suptitle('{} vs {}'.format(titlename, get_tex(db, x1[0])))

    prettify_plot5(fmeta=f)
    figname = "{}/{}_measurement({})_sub({})_line({})".format(meta['plotdir'], meta['plotprefix'], x1[0], sub1[0], l1[0])
    f['fig'].savefig(figname + '.pdf', format='pdf')
    f['fig'].savefig(figname + '.png', format='png')
    f['fig'].savefig(figname + '.svg', format='svg')
    tikzplotlib.save('{}.tex'.format(figname))
