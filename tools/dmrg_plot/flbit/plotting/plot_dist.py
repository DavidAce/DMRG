from itertools import product
from pathlib import Path

import h5py
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


def plot_dist_fig3_sub3_line1(db, meta, figspec, subspec, linspec, algo_filter=None, state_filter=None, point_filter=None, figs=None, palette_name=None):
    if len(figspec) + len(subspec) + len(linspec) != 7:
        raise AssertionError("Must add to 7 elems: \n figspec {}\n subspec {}\n linespec {}")
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
            if not col in [l.split(':')[0] for l in subspec]:
                legend_col_keys.append(col)

    figprod = list(product(*get_keys(db, figspec)))  # All combinations of figspecs (names of parameters that iterate figures)
    subprod = list(product(*get_keys(db, subspec)))  # All combinations of subspecs (names of parameters that iterate subplots)
    linprod = list(product(*get_keys(db, linspec)))  # All combinations of linspecs (names of parameters that iterate lines)
    dirprod = list(product(db['keys']['algo'], db['keys']['state'], db['keys']['crono']))
    numfigs = len(figprod)
    numsubs = len(subprod)
    if figs is None:
        figs = [get_fig_meta(numsubs, meta=meta) for _ in range(numfigs)]

    for figkeys, f in zip(figprod, figs):
        print('- plotting figkeys: {}'.format(figkeys))
        dbval = None
        for idx, (subkeys, ax) in enumerate(zip(subprod, f['ax'])):
            popt = None
            pcov = None
            dbval = None
            print('-- plotting subkeys: {}'.format(subkeys))
            for dirkeys in dirprod:
                palette, lstyles = get_colored_lstyles(db, linspec, palette_name)
                for linkeys, color, lstyle in zip(linprod, palette, lstyles):
                    findlist = list(figkeys) + list(subkeys) + list(linkeys) + list(dirkeys) + [meta['groupname']]
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

                    tidxs = meta.get('tidx')
                    if isinstance(tidxs, str) and 'window' in tidxs:
                        tdata = mmntnode['avg']['physical_time']
                        ydata = mmntnode['avg']['entanglement_entropy']
                        widxs = find_loglog_window2(tdata=tdata, ydata=ydata, db=dbval)
                        # tidxs = [1, 5, widxs[0], len(tdata) - 1]
                        tidxs = [len(tdata) - 1]
                    # palette = sns.color_palette(palette=palette_name, n_colors=len(meta['tidx']))
                    for i, tidx in enumerate(tidxs):
                        data = datanode[meta['dsetname']][tidx, :]
                        dbval['vals']['t'] = mmntnode['avg']['physical_time'][tidx]
                        dbval['vals']['bavg'] = mmntnode['avg']['bond_mid'][tidx]
                        dbval['vals']['num'] = np.shape(datanode[meta['dsetname']])[1]

                        if meta.get('normpage') == True:
                            p = page_entropy(dbval['L'])
                            data /= page_entropy(dbval['L'])
                        if 'normalize' in meta:
                            data /= meta['normalize']

                        hist, edges = np.histogram(data, bins=meta['bins'], density=True)
                        bincentres = [(edges[j] + edges[j + 1]) / 2. for j in range(len(edges) - 1)]
                        line, = ax.step(x=bincentres, y=hist, where='mid', label=None,
                                        color=color, path_effects=path_effects, linestyle=lstyle)
                        if 'number_' in meta['dsetname']:
                            ax.axvline(x=np.log(2), color='grey')
                            ax.axvline(x=np.log(3), color='darkseagreen')
                            # ax.set_xticks(list(ax.get_xticks()) + [np.log(2), np.log(3)])
                            trans = transforms.blended_transform_factory(ax.transData, ax.transAxes)
                            ax.text(np.log(2), 0.8, '$\ln 2$', fontsize='small', color='grey', ha='right', va='center', rotation='vertical',
                                    transform=trans)
                            ax.text(np.log(3), 0.8, '$\ln 3$', fontsize='small', color='darkseagreen', ha='right', va='center', rotation='vertical',
                                    transform=trans)
                        if 'algorithm_time' in meta['dsetname']:
                            t_mean = np.mean(data)
                            trans = transforms.blended_transform_factory(ax.transData, ax.transAxes)
                            ax.axvline(x=t_mean, color=color)
                            ax.text(t_mean, 0.8, 'avg', fontsize='small', color='grey', ha='right', va='center', rotation='vertical',
                                    transform=trans)
                        legendrow = get_legend_row(db=db, datanode=datanode, legend_col_keys=legend_col_keys)
                        for icol, (col, key) in enumerate(zip(legendrow, legend_col_keys)):
                            key, fmt = key.split(':') if ':' in key else [key, '']
                            f['legends'][idx][icol]['handle'].append(line)
                            f['legends'][idx][icol]['title'] = db['tex'][key]
                            f['legends'][idx][icol]['label'].append(col)

                        if 'number' in meta['dsetname'] and tidx == tidxs[-1] and 'L_16' in findlist:
                            # Plot Luitz data
                            with h5py.File('external/raw_EE_NE_CE_distributions_random_XXX_chain.h5', 'r') as h5ext:
                                hist = h5ext['L16/W8.0']['hist[NE1][100]'][()]
                                edges = h5ext['L16/W8.0']['binedges[NE1][100]'][()]
                                bincentres = [(edges[j] + edges[j + 1]) / 2. for j in range(len(edges) - 1)]
                                line_ext, = ax.step(x=bincentres, y=hist, where='mid', label=None, color='black')
                                # for icol, (col, key) in enumerate(zip(legendrow, legend_col_keys)):
                                #     key, fmt = key.split(':') if ':' in key else [key, '']
                                #     f['legends'][idx][icol]['handle'].append(line_ext)
                                #     f['legends'][idx][icol]['title'] = db['tex'][key]
                                #     f['legends'][idx][icol]['label'].append('PhysRevB.102.100202')

                        if not idx in f['axes_used']:
                            f['axes_used'].append(idx)
            if dbval:
                ax.set_title(get_title(dbval, subspec),
                             horizontalalignment='left', x=0.05,
                             fontstretch="ultra-condensed",
                             # bbox=dict(boxstyle='square,pad=0.15', facecolor='white', alpha=0.6)
                             )

        if not prb_style and dbval:
            f['fig'].suptitle('{} distribution\n{}'.format(meta['titlename'], get_title(dbval, figspec)))

        # prettify_plot4(fmeta=f, lgnd_meta=axes_legends)
        if not f['filename']:
            suffix = ''
            suffix = suffix + '_normpage' if 'normpage' in meta and meta['normpage'] else suffix
            suffix = suffix + '_loglog' if 'timeloglevel' in meta and meta['timeloglevel'] >= 2 else suffix
            f['filename'] = "{}/{}_dist_fig({})_sub({}){}".format(meta['plotdir'], meta['plotprefix'],
                                                                  '-'.join(map(str, figkeys)),
                                                                  '-'.join(map(str, get_keys(db, subspec))),
                                                                  suffix)

    return figs
