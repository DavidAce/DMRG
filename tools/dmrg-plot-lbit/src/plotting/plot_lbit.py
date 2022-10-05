from itertools import product
from pathlib import Path
import seaborn as sns
from .tools import *


def plot_v2_lbit_fig3_sub3_line1(db, meta, fig3, sub3, l1, algo_filter=None, state_filter=None, f=None, palette_name=None):
    if len(fig3) != 3:
        raise AssertionError("fig must have length 3")
    if len(sub3) != 3:
        raise AssertionError("sub must have length 3")
    if len(l1) != 1:
        raise AssertionError("itr must have length 1")
    if 'mplstyle' in meta:
        plt.style.use(meta['mplstyle'])
    if 'plotdir' in meta and 'mplstyle' in meta:
        if Path(meta['plotdir']).stem != Path(meta['mplstyle']).stem:
            meta['plotdir'] = Path(meta['plotdir'], Path(meta['mplstyle']).stem)
            Path(meta['plotdir']).mkdir(parents=True, exist_ok=True)
            print("Setting plotdir: ", meta['plotdir'])
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
    legend_col_keys = l1.copy()
    if legendcols := meta['legendcols']:
        for col in legendcols:
            if not col in [l.split(':')[0] for l in sub3 + l1]:
                legend_col_keys.append(col)

    figprod = list(product(*get_keys(db, fig3)))
    for key0, key1, key2 in figprod:
        subprod = list(product(*get_keys(db, sub3)))
        numplots = len(subprod)
        if not f:
            f = get_fig_meta(numplots, meta=meta)
        for idx, ((key3, key4, key5), ax) in enumerate(zip(subprod, f['ax'])):
            axin = None
            popt = None
            pcov = None
            dbval = None
            for algokey, modelkey in product(db['keys']['algo'], db['keys']['model']):
                palette = sns.color_palette(palette=palette_name, n_colors=len(db['keys'][l1[0]]))
                for key6, color in zip(get_keys(db, l1[0]), palette):
                    findlist = [key0, key1, key2, key3, key4, key5, key6, algokey, modelkey, meta['dsetname']]
                    datanode = [value['node']['data'] for key, value in db['dsets'].items() if
                                all(k in key for k in findlist)]
                    if len(datanode) != 1:
                        print("ERROR: found", len(datanode), "datanodes: ", datanode, " | findlist: ", findlist)
                        continue
                        # raise LookupError("Found incorrect number of datanodes")
                    datanode = datanode[0]
                    dbval = db['dsets'][datanode.name]
                    ydata, _ = get_table_data(datanode['avg'])
                    edata, _ = get_table_data(datanode['ste'])
                    ndata = datanode['num'][()]
                    xdata = range(len(ydata))
                    if np.min(ndata) < 10:
                        continue
                    for i, (y, e) in enumerate(zip(ydata.T, edata.T)):
                        ax.fill_between(x=xdata, y1=y - e, y2=y + e, alpha=0.10, color=color)
                        line, = ax.plot(xdata, y, marker=None, color=color, path_effects=path_effects)
                        if i == 0:
                            legendrow = get_legend_row(db=db, datanode=datanode, legend_col_keys=legend_col_keys)
                            for icol, (col, key) in enumerate(zip(legendrow, legend_col_keys)):
                                key, fmt = key.split(':') if ':' in key else [key, '']
                                f['legends'][idx][icol]['handle'].append(line)
                                f['legends'][idx][icol]['label'].append(col)
                                f['legends'][idx][icol]['title'] = db['tex'][key]

                    if not idx in f['axes_used']:
                        f['axes_used'].append(idx)
            if dbval:
                ax.set_title(get_title(dbval, sub3), x=0, horizontalalignment='left')
            ax.set_xlabel("$|i-j|$")
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
            f['fig'].suptitle('{}\n{}'.format(meta['titlename'], get_title(dbval, fig3)))

        # prettify_plot4(fmeta=f, lgnd_meta=axes_legends)
        if not f['filename']:
            f['filename'] = "{}/{}(t)_fig({}_{}_{})_sub({}_{}_{})".format(meta['plotdir'], meta['plotprefix'],
                                                                          str(key0), str(key1), str(key2), sub3[0], sub3[1], sub3[2])

    return f
