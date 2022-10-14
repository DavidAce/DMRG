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


def plot_dist_fig4_sub3(db, meta, fig4, sub3, algo_filter=None, state_filter=None, point_filter=None, f=None, palette_name=None):
    if len(fig4) != 4:
        raise AssertionError("fig must have length 4")
    if len(sub3) != 3:
        raise AssertionError("sub must have length 3")
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
            if not col in [l.split(':')[0] for l in sub3]:
                legend_col_keys.append(col)

    for key0 in get_keys(db, fig4[0]):
        for key1 in get_keys(db, fig4[1]):
            for key2 in get_keys(db, fig4[2]):
                for key3 in get_keys(db, fig4[3]):

                    keyprod = list(product(*get_keys(db, sub3)))
                    numplots = len(keyprod)
                    if not f:
                        f = get_fig_meta(numplots, meta=meta)
                    for idx, ((key4, key5, key6), ax) in enumerate(zip(keyprod, f['ax'])):
                        popt = None
                        pcov = None
                        dbval = None
                        for algokey, statekey, cronokey in product(db['keys']['algo'], db['keys']['state'], db['keys']['crono']):
                            # palette = plt.rcParams['axes.prop_cycle'].by_key()['color']
                            palette = sns.color_palette(palette=palette_name, n_colors=len(meta['tidx']))
                            findlist = [key0, key1, key2, key3, key4, key5, key6, algokey, statekey, cronokey,
                                        meta['groupname']]
                            datanode = [value['node']['data'] for key, value in db['dsets'].items() if
                                        all(k in key for k in findlist)]
                            if len(datanode) != 1:
                                print("ERROR: found", len(datanode), "datanodes: ", datanode, " | findlist: ", findlist)
                                continue
                                # raise LookupError("Found incorrect number of datanodes")
                            datanode = datanode[0]
                            mmntnode = datanode.parent['measurements']
                            ndata = datanode['avg']['num'][()]
                            if np.min(ndata) < 10:
                                continue
                            for i, (tidx, color) in enumerate(zip(meta['tidx'], palette)):
                                dbval = db['dsets'][datanode.name]
                                data = datanode['data'][tidx, :]
                                dbval['vals']['t'] = mmntnode['avg']['physical_time'][tidx]
                                dbval['vals']['bavg'] = mmntnode['avg']['bond_mid'][tidx]
                                dbval['vals']['num'] = np.shape(datanode['data'])[1]

                                if meta.get('normpage') == True:
                                    p = page_entropy(dbval['L'])
                                    data /= page_entropy(dbval['L'])
                                if 'normalize' in meta:
                                    data /= meta['normalize']

                                hist, edges = np.histogram(data, bins=meta['bins'], density=True)
                                bincentres = [(edges[j] + edges[j + 1]) / 2. for j in range(len(edges) - 1)]
                                line, = ax.step(x=bincentres, y=hist, where='mid', label=None,
                                                color=color, path_effects=path_effects)

                                legendrow = get_legend_row(db=db, datanode=datanode, legend_col_keys=legend_col_keys)
                                for icol, (col, key) in enumerate(zip(legendrow, legend_col_keys)):
                                    key, fmt = key.split(':') if ':' in key else [key, '']
                                    print(icol, col, key, fmt)
                                    f['legends'][idx][icol]['handle'].append(line)
                                    f['legends'][idx][icol]['title'] = db['tex'][key]
                                    f['legends'][idx][icol]['label'].append(col)

                                    # if 'num' in key:
                                    #
                                    # else:
                                    #     f['legends'][idx][icol]['title'] = db['tex'][key]

                                if not idx in f['axes_used']:
                                    f['axes_used'].append(idx)
                        if dbval:
                            ax.set_title(get_title(dbval, sub3), x=0, horizontalalignment='left')
                        if xlabel := meta.get('xlabel'):
                            ax.set_xlabel(xlabel)

                        if ymin := meta.get('ymin'):
                            f['ymin'] = ymin
                            ax.set_ylim(ymin=ymin)
                        if ymax := meta.get('ymax'):
                            f['ymax'] = ymax
                            ax.set_ylim(ymax=ymax)
                        if xmin := meta.get('xmin'):
                            f['xmin'] = xmin
                            ax.set_xlim(xmin=xmin)
                        if xmax := meta.get('xmax'):
                            f['xmax'] = xmax
                            ax.set_xlim(xmax=xmax)

                if not prb_style and dbval:
                    f['fig'].suptitle('{} distribution\n{}'.format(meta['titlename'], get_title(dbval, fig4)))

                # prettify_plot4(fmeta=f, lgnd_meta=axes_legends)
                if not f['filename']:
                    suffix = ''
                    suffix = suffix + '_normpage' if 'normpage' in meta and meta['normpage'] else suffix
                    suffix = suffix + '_loglog' if 'timeloglevel' in meta and meta['timeloglevel'] >= 2 else suffix
                    f['filename'] = "{}/{}_dist_fig({}_{}_{}_{})_sub({}_{}_{}){}".format(meta['plotdir'], meta['plotprefix'],
                                                                                         str(key0), str(key1), str(key2), str(key3), sub3[0], sub3[1], sub3[2],
                                                                                         suffix)

    return f
