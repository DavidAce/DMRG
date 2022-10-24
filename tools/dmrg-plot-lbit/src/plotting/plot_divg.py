from itertools import product
from pathlib import Path
import matplotlib.transforms as transforms
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


def plot_divg_fig3_sub3_line1(db, meta, fig3, sub3, l1, algo_filter=None, state_filter=None, point_filter=None, figs=None, palette_name=None):
    if len(fig3) != 3:
        raise AssertionError("fig must have length 3")
    if len(sub3) != 3:
        raise AssertionError("sub must have length 3")
    if len(l1) != 1:
        raise AssertionError("l must have length 1")
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

    for key0 in get_keys(db, fig3[0]):
        for key1 in get_keys(db, fig3[1]):
            for key2 in get_keys(db, fig3[2]):

                keyprod = list(product(*get_keys(db, sub3)))
                numplots = len(keyprod)
                if figs is None:
                    figs = []
                figs.append(get_fig_meta(numplots, meta=meta))
                f = figs[-1]
                for idx, ((key3, key4, key5), ax) in enumerate(zip(keyprod, f['ax'])):
                    popt = None
                    pcov = None
                    dbval = None
                    for algokey, statekey, cronokey in product(db['keys']['algo'], db['keys']['state'], db['keys']['crono']):
                        palette = sns.color_palette(palette=palette_name, n_colors=len(db['keys'][l1[0]]))
                        for key6, color in zip(get_keys(db, l1[0]), palette):
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
                            ndata = datanode['avg']['num'][()]
                            if np.min(ndata) < 10:
                                continue
                            t = mmntnode['avg']['physical_time']
                            y = datanode['data']
                            n = datanode['avg']['num'][()][0]
                            if meta.get('normpage'):
                                y /= page_entropy(dbval['L'])
                            if normalize := meta.get('normalize'):
                                y /= normalize
                            # Calculate the infinite time average (1/T) integral_0^T y(t) dt
                            ytavg = np.trapz(y=y, x=t, axis=0) / t[-1]

                            hist, edges = np.histogram(ytavg, bins=meta['bins'], density=True)
                            bincentres = [(edges[j] + edges[j + 1]) / 2. for j in range(len(edges) - 1)]
                            line, = ax.step(x=bincentres, y=hist, where='mid', label=None,
                                            color=color, path_effects=path_effects)
                            ax.axvline(x=np.log(2), color='grey')
                            ax.axvline(x=np.log(3), color='darkseagreen')
                            # ax.set_xticks(list(ax.get_xticks()) + [np.log(2), np.log(3)])
                            trans = transforms.blended_transform_factory(ax.transData, ax.transAxes)
                            ax.text(np.log(2), 0.8, '$\ln 2$', fontsize='small', color='grey', ha='right', va='center', rotation='vertical', transform=trans)
                            ax.text(np.log(3), 0.8, '$\ln 3$', fontsize='small', color='darkseagreen', ha='right', va='center', rotation='vertical',
                                    transform=trans)
                            legendrow = get_legend_row(db=db, datanode=datanode, legend_col_keys=legend_col_keys)
                            for icol, (col, key) in enumerate(zip(legendrow, legend_col_keys)):
                                key, fmt = key.split(':') if ':' in key else [key, '']
                                f['legends'][idx][icol]['handle'].append(line)
                                f['legends'][idx][icol]['title'] = db['tex'][key]
                                f['legends'][idx][icol]['label'].append(col)

                            if not idx in f['axes_used']:
                                f['axes_used'].append(idx)

                    if dbval:
                        ax.set_title(get_title(dbval, sub3),
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
                    f['filename'] = "{}/{}_divg_fig({}_{}_{})_sub({}_{}_{}){}".format(meta['plotdir'], meta['plotprefix'],
                                                                                      str(key0), str(key1), str(key2), sub3[0], sub3[1], sub3[2],
                                                                                      suffix)

    return figs
