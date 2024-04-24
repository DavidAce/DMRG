from matplotlib.ticker import MaxNLocator
# from src.database.database import *
from .tools import *
from src.general.filter import *
from itertools import product
from copy import deepcopy
from .style import *
import numpy as np
import seaborn as sns
import matplotlib.patheffects as pe
from pathlib import Path
import tikzplotlib
from scipy.optimize import curve_fit


def S_vs_l_effective_central_charge(x, ceff, c1):
    L = len(x) - 1
    L2 = int(L / 2)
    S = [0] * len(x)
    S[1:L2] = ceff / 3.0 * np.log(x[1:L2]) + c1
    S[L2:-1] = ceff / 3.0 * np.log(x[L2:0:-1]) + c1
    # print("len S: ", len(S))
    # print(x)
    # print(S)
    # if len(x) != len(S):
    #     raise AssertionError("x {} != S {}".format(len(x), len(S)))
    return S


def S_vs_l_effective_central_charge_folded(x, ceff, c1):
    S = [0] * len(x)
    S[1:] = ceff / 3.0 * np.log(x[1:]) + c1
    return S


def S_vs_l_effective_central_charge_v2(x, ceff, c1):
    L = len(x) - 1
    L2 = int(L / 2)
    S = [0] * len(x)
    S[1:L2] = ceff / 3.0 * np.log((2 * L / np.pi) * np.sin((np.pi / L) * np.asarray(x[1:L2]))) + c1
    S[L2:-1] = ceff / 3.0 * np.log((2 * L / np.pi) * np.sin((np.pi / L) * np.asarray(x[L2:0:-1]))) + c1
    return S


def S_vs_l_effective_central_charge_v2_folded(x, ceff, c1):
    L = 2 * (len(x) - 1)
    S = [0] * len(x)
    S[1:] = ceff / 3.0 * np.log((2 * L / np.pi) * np.sin((np.pi / L) * np.asarray(x[1:]))) + c1
    return S


def plot_table(h5, db, meta, sub1, l1):
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
        for algokey in db['keys']['algo']:
            for statekey in db['keys']['state']:
                lstyle = next(lstyles)
                mstyles = get_markerstyles()
                dbval = None
                for pointkey in db['keys']['point']:
                    ed_palette = itertools.cycle(sns.color_palette("Set2"))
                    xd_palette = itertools.cycle(sns.color_palette(palette_name, len(get_keys(db, l1[0]))))
                    mstyle = next(mstyles)
                    for key1 in get_keys(db, l1[0]):
                        findlist = [key0, key1, algokey, statekey, pointkey]
                        ygrps = [value['path']['data'] for key, value in db['dsets'].items() if
                                 all(meta['ygrp'] in key and '/{}/'.format(k) in key for k in findlist)]
                        if not ygrps:
                            raise KeyError("Could not find ygrp [{}] in findlist: {}".format(meta['ygrp'], findlist))

                        # We can plot full chain quantities, such as S vs L.
                        # To do that we set meta['xgrp'] = None
                        xgrps = [None]
                        if meta.get('xgrp') is not None:
                            xgrps = [value['path']['data'] for key, value in db['dsets'].items() if
                                     all(meta['xgrp'] in key and '/{}/'.format(k) in key for k in findlist)]
                            # if not xgrps:
                        #     raise KeyError("Could not find xgrps: \n findlist: {}".format(findlist))
                        for ygrp, xgrp in itertools.product(ygrps, xgrps):
                            dbval = db['dsets'][ygrp]
                            ldata = get_legend_row(db=db, datanode=h5[ygrp], legend_col_keys=legend_col_keys)
                            ydata = get_data(h5[ygrp]['avg'], meta['ycol'], 'f8')
                            edata = get_data(h5[ygrp]['ste'], meta['ycol'], 'f8')
                            xdata = list(range(len(ydata)))
                            if xgrp is not None:
                                xdata = get_data(h5[xgrp]['avg'], meta['xcol'], 'f8')
                                sort = np.argsort(xdata)
                                xdata = xdata[sort]
                                ydata = ydata[sort]
                                edata = edata[sort]

                            if xfold := meta.get('xfold'):
                                L = int(len(xdata) - 1)
                                L2 = int(L / 2)
                                xdata = xdata[0:L2 + 1]
                                ydata = 0.5 * np.asarray(ydata[0:L2 + 1] + ydata[L:L2 - 1:-1])
                                edata = 0.5 * np.asarray(edata[0:L2 + 1] + edata[L:L2 - 1:-1])

                            if "states" in statekey:
                                color = next(ed_palette)
                                mstyle = 'o'
                            else:
                                color = next(xd_palette)
                            line, _, _ = ax.errorbar(x=xdata, y=ydata, yerr=edata, color=color, marker=mstyle, linestyle=lstyle,
                                                     path_effects=path_effects)
                            for icol, (col, key) in enumerate(zip(ldata, legend_col_keys)):
                                key, fmt = key.split(':') if ':' in key else [key, '']
                                f['legends'][idx][icol]['handle'].append(line)
                                f['legends'][idx][icol]['label'].append(col)
                                f['legends'][idx][icol]['title'] = db['tex'][key]

                            if yfit := meta.get('yfit'):
                                if yfit == 'S_vs_l_effective_central_charge':
                                    bounds = ([0, 0], [10, 10])
                                    fit = S_vs_l_effective_central_charge_v2
                                    if xfold := meta.get('xfold'):
                                        fit = S_vs_l_effective_central_charge_v2_folded
                                    popt, pcov = curve_fit(f=fit,
                                                           xdata=xdata,
                                                           ydata=ydata,
                                                           bounds=bounds)
                                    print(key0, key1, popt)
                                    ax.plot(xdata, fit(xdata, *popt), marker=None,
                                            linewidth=0.8,
                                            linestyle='--', label=None, color=color,
                                            path_effects=path_effects)

                                    icol = len(legend_col_keys)
                                    f['legends'][idx][icol + 0]['handle'].append(line)
                                    f['legends'][idx][icol + 0]['label'].append('{:.2f}'.format(popt[0]))
                                    f['legends'][idx][icol + 0]['title'] = '$c_\mathrm{eff}$'
                                    f['legends'][idx][icol + 1]['handle'].append(line)
                                    f['legends'][idx][icol + 1]['label'].append('{:.2f}'.format(popt[1]))
                                    f['legends'][idx][icol + 1]['title'] = '$c_1$'

                            if not idx in f['axes_used']:
                                f['axes_used'].append(idx)
                if dbval is not None:
                    ax.set_title(get_title(dbval, sub1), x=0.5, horizontalalignment='left')
                if ylabel := meta.get('ylabel'):
                    ax.set_ylabel(ylabel)
                if xlabel := meta.get('xlabel'):
                    ax.set_xlabel(xlabel)
                if xlim := meta.get('xlim', {}):
                    ax.set_xlim(xmin=xlim[0], xmax=xlim[1], auto=True)
                if ylim := meta.get('ylim'):
                    ax.set_ylim(ymin=ylim[0], ymax=ylim[1], auto=True)

    if titlename := meta.get('titlename'):
        f['fig'].suptitle(titlename)

    prettify_plot5(fmeta=f)
    suffix = "" if not meta.get('xfold') else '_folded'
    figname = "{}/{}_table({}_vs_{})_sub({})_line({}){}".format(meta['plotdir'], meta['plotprefix'], meta['ytag'], meta['xtag'], sub1[0], l1[0], suffix)
    f['fig'].savefig('{}.pdf'.format(figname), format='pdf')
    f['fig'].savefig('{}.png'.format(figname), format='png')
    f['fig'].savefig('{}.svg'.format(figname), format='svg')
    f['fig'].savefig('{}.pgf'.format(figname), format='pgf')
    tikzplotlib.save('{}.tex'.format(figname))
