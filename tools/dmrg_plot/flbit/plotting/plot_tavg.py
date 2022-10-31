from itertools import product
from pathlib import Path

import matplotlib.patheffects as pe
import seaborn as sns
from scipy.optimize import curve_fit

from .tools import *


# def floglog(x, a, b,c):
#     return a + b*np.log(np.log(x + c))
def floglog_v2(x, a, b, c):
    with np.errstate(invalid='ignore'):
        return a + b * np.log(np.log(x - c))


def flinear(x, a, b):
    with np.errstate(invalid='ignore'):
        return a + b * x


def fpower(x, a, b):
    with np.errstate(invalid='ignore'):
        return a * x ** b


def plot_tavg_fig3_sub2_line1(db, meta, fig3, sub2, l1, x1, algo_filter=None, state_filter=None, point_filter=None, figs=None, palette_name=None):
    if len(fig3) != 3:
        raise AssertionError("fig must have length 3")
    if len(sub2) != 2:
        raise AssertionError("sub must have length 2")
    if len(l1) != 1:
        raise AssertionError("l1 must have length 1")
    if len(x1) != 1:
        raise AssertionError("x1 must have length 1")
    if 'mplstyle' in meta:
        plt.style.use(meta['mplstyle'])
    if 'plotdir' in meta and 'mplstyle' in meta:
        if Path(meta['plotdir']).stem != Path(meta['mplstyle']).stem:
            meta['plotdir'] = Path(meta['plotdir'], Path(meta['mplstyle']).stem)
            Path(meta['plotdir']).mkdir(parents=True, exist_ok=True)
    if 'mplstyle' in meta and 'slack' in meta['mplstyle']:
        # palette_name = "Spectral"
        if not palette_name:
            palette_name = "colorblind"
        # path_effects = [pe.SimpleLineShadow(offset=(0.5, -0.5), alpha=0.3), pe.Normal()]
        path_effects = None
    else:
        if not palette_name:
            palette_name = "colorblind"
        path_effects = None

    prb_style = 'prb' in meta['mplstyle'] if 'mplstyle' in meta else False

    # legend_col_keys = list(itertools.chain(l1, [col for col in meta['legendcols'] if 'legendcols' in meta]))
    legend_col_keys = []
    if legendcols := meta['legendcols']:
        for col in legendcols:
            if not col in [l.split(':')[0] for l in sub2]:
                legend_col_keys.append(col)

    for key0 in get_keys(db, fig3[0]):
        for key1 in get_keys(db, fig3[1]):
            for key2 in get_keys(db, fig3[2]):

                keyprod = list(product(*get_keys(db, sub2)))
                numplots = len(keyprod)
                if figs is None:
                    figs = []
                figs.append(get_fig_meta(numplots, meta=meta))
                f = figs[-1]
                # print(f['ax'])
                # print(keyprod)

                for idx, ((key3, key4), ax) in enumerate(zip(keyprod, f['ax'])):
                    popt = None
                    pcov = None
                    dbval = None
                    for algokey, statekey, cronokey in product(db['keys']['algo'], db['keys']['state'], db['keys']['crono']):
                        palette = sns.color_palette(palette=palette_name, n_colors=len(db['keys'][l1[0]]))
                        for key5, color in zip(get_keys(db, l1[0]), palette):
                            x1len = len(get_keys(db, x1[0]))
                            ydata = np.empty(shape=(x1len,), dtype=np.float)
                            xdata = np.empty(shape=(x1len,), dtype=np.float)
                            edata = np.empty(shape=(x1len,), dtype=np.float)
                            ndata = np.empty(shape=(x1len,), dtype=np.int)
                            for xidx, key6 in enumerate(get_keys(db, x1[0])):
                                findlist = [key0, key1, key2, key3, key4, key5, key6, algokey, statekey, cronokey, meta['groupname']]
                                datanode = [value['node']['data'] for key, value in db['dsets'].items() if
                                            all(k in key for k in findlist)]
                                if len(datanode) != 1:
                                    print("ERROR: found", len(datanode), "datanodes: ", datanode, " | findlist: ", findlist)
                                    continue
                                    # raise LookupError("Found incorrect number of datanodes")
                                datanode = datanode[0]
                                mmntnode = datanode.parent['measurements']
                                dbval = db['dsets'][datanode.name]
                                # Get this datapoint
                                t = mmntnode['avg']['physical_time'][()]
                                s = mmntnode['avg']['entanglement_entropy_midchain'][()]
                                y = datanode['data'][()]
                                n = datanode['avg']['num'][()][0]
                                idx_sat = find_saturation_idx3(t, s, dbval)

                                if meta.get('normpage'):
                                    y /= page_entropy(dbval['L'])
                                if normalize := meta.get('normalize'):
                                    y /= normalize
                                # Calculate the infinite time average (1/T) integral_0^T y(t) dt, where [0,T] is
                                # a range where y has saturated (i.e. after the transient)
                                ytavg = np.trapz(y=y[idx_sat:], x=t[idx_sat:], axis=0) / (t[-1] - t[idx_sat])
                                # Calculate the infinite time average (1/T) integral_0^T y(t) dt
                                ydata[xidx] = np.mean(ytavg)
                                xdata[xidx] = get_vals(dbval, x1[0])
                                edata[xidx] = np.std(ytavg) / np.sqrt(len(ytavg))
                                ndata[xidx] = n

                            line = ax.errorbar(x=xdata, y=ydata, yerr=edata, color=color, path_effects=path_effects)
                            # ax.fill_between(x=xdata, y1=ydata - edata, y2=ydata + edata, alpha=0.10, color=color)
                            # line, = ax.plot(xdata, ydata, marker=None, color=color, path_effects=path_effects)

                            legendrow = get_legend_row(db=db, datanode=datanode, legend_col_keys=legend_col_keys)
                            for icol, (col, key) in enumerate(zip(legendrow, legend_col_keys)):
                                key, fmt = key.split(':') if ':' in key else [key, '']
                                f['legends'][idx][icol]['handle'].append(line)
                                f['legends'][idx][icol]['title'] = db['tex'][key]
                                f['legends'][idx][icol]['label'].append(col)

                            if not idx in f['axes_used']:
                                f['axes_used'].append(idx)

                    if dbval:
                        ax.set_title(get_title(dbval, sub2),
                                     horizontalalignment='left', x=0.05,
                                     fontstretch="ultra-condensed",
                                     # bbox=dict(boxstyle='square,pad=0.15', facecolor='white', alpha=0.6)
                                     )

                if not prb_style and dbval:
                    f['fig'].suptitle('{} distribution\n{}'.format(meta['titlename'], get_title(dbval, fig3)))

                # prettify_plot4(fmeta=f, lgnd_meta=axes_legends)
                if not f['filename']:
                    suffix = ''
                    suffix = suffix + '_normpage' if 'normpage' in meta and meta['normpage'] else suffix
                    suffix = suffix + '_loglog' if 'timeloglevel' in meta and meta['timeloglevel'] >= 2 else suffix
                    f['filename'] = "{}/{}_tavg_fig({}_{}_{})_sub({}_{}){}".format(meta['plotdir'], meta['plotprefix'],
                                                                                   str(key0), str(key1), str(key2), sub2[0], sub2[1],
                                                                                   suffix)

    return figs
