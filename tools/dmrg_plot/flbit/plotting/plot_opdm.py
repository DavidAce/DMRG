from itertools import product
from pathlib import Path
import seaborn as sns
from .tools import *
import warnings
import logging
import matplotlib.patheffects as pe

logger = logging.getLogger('plot-opdm')

def format_err(val, err, min_decimals=None):
    err_exponent =  np.floor(np.log10(err))
    num_decimals =  int(np.abs(err_exponent))
    if min_decimals is not None:
        num_decimals = min_decimals
    return f'{val:.{num_decimals}f}({err*10**num_decimals:1.0f})'

def plot_v3_opdm_fig3_sub3_line1(db, meta, figspec, subspec, linspec, algo_filter=None, state_filter=None, figs=None,
                                 palette_name=None,dbidx=0,dbnum=1):
    if db['version'] != 3:
        raise ValueError("plot_v3_opdm_fig3_sub3_line1 requires db version 3")
    path_effects = [pe.SimpleLineShadow(offset=(-0.35, -0.35), alpha=0.4), pe.Normal()]
    path_effects_dashed = [pe.SimpleLineShadow(offset=(-0.1, -0.1), alpha=0.1), pe.Normal()]
    if 'mplstyle' in meta:
        plt.style.use(meta['mplstyle'])
        if 'slack' in meta['mplstyle']:
            path_effects = [pe.SimpleLineShadow(offset=(-0.5, -0.5), alpha=0.3), pe.Normal()]
            path_effects_dashed = [pe.SimpleLineShadow(offset=(-0.1, -0.1), alpha=0.1), pe.Normal()]

    legend_col_keys = [l for l in linspec if l in meta['legendcols']]
    if legendcols := meta['legendcols']:
        for col in legendcols:
            if not col in [l.split(':')[0] for l in subspec + linspec]:
                legend_col_keys.append(col)

    figprod = list(product(*get_vals(db=db, keyfmt=figspec, filter=meta.get('filter'))))  # All combinations of figspecs values
    subprod = list(product(*get_vals(db=db, keyfmt=subspec, filter=meta.get('filter'))))  # All combinations of subspecs values
    linprod = list(product(*get_vals(db=db, keyfmt=linspec, filter=meta.get('filter'))))  # All combinations of linspecs values
    numfigs = len(figprod)
    numsubs = len(subprod)
    if figs is None:
        figs = [get_fig_meta(numsubs, meta=meta) for _ in range(numfigs)]

    for figvals, f in zip(figprod, figs):
        logger.debug(f'-- plotting {figvals=}')
        dbval = None
        for idx, (subvals, ax, ix) in enumerate(zip(subprod, f['ax'], f['ix'])):
            logger.debug('-- plotting subkeys: {}'.format(subvals))
            palette, lstyles = get_colored_lstyles(db, linspec, palette_name, meta.get('filter'))
            for linvals, color, lstyle in zip(linprod, palette, lstyles):
                logger.debug('--- plotting linkeys: {}'.format(linvals))
                print('--- plotting linkeys: {}'.format(linvals))
                datanodes = match_datanodes(db=db, meta=meta, specs=figspec + subspec + linspec,
                                            vals=figvals + subvals + linvals)
                if len(datanodes) != 1:
                    logger.warning(f"match: \n"
                                   f"\tspec:{[figspec + subspec + linspec]}\n"
                                   f"\tvals:{[figvals + subvals + linvals]}")
                    logger.warning(f"found {len(datanodes)} datanodes: {datanodes=}")
                    continue
                for datanode in datanodes:
                    # datanode = datanode[0]
                    dbval = db['dsets'][datanode.name]
                    if dbval['vals']['f'] == 0.3:
                        continue
                    # if dbval['vals']['L'] != 20:
                    #     continue

                    # ndata = dbval['vals']['num']
                    Ldata = dbval['vals']['L']
                    # xdata = np.array(range(Ldata))
                    ndata = np.shape(datanode[()])[1]
                    ydata = np.mean(datanode[()], axis=1)
                    edata = np.std(datanode[()], axis=1)/np.sqrt(ndata)
                    xdata = np.arange(len(ydata))
                    # line = ax.errorbar(xdata, ydata, marker='o', color=color, path_effects=path_effects, zorder=1)
                    line = ax.errorbar(x=xdata, y=ydata, yerr=edata, linestyle='none',capsize=4.0,
                                       marker='o', color=color, path_effects=path_effects)

                    legendrow = get_legend_row(db=db, datanode=datanode, legend_col_keys=legend_col_keys)
                    for icol, (col, key) in enumerate(zip(legendrow, legend_col_keys)):
                        key, fmt = key.split(':') if ':' in key else [key, '']
                        f['legends'][idx][icol]['handle'].append(line)
                        f['legends'][idx][icol]['label'].append(col)
                        f['legends'][idx][icol]['title'] = db['tex'][key]

            if not idx in f['axes_used']:
                f['axes_used'].append(idx)

            if subspec_title := get_subspec_title(meta,dbval,subspec):
                ax.set_title(subspec_title,horizontalalignment='left', x=0.05,fontstretch="ultra-condensed")

        if figspec_title := get_figspec_title(meta,dbval,figspec):
            f['fig'].suptitle(figspec_title, y=0.66, x=0.85)

        f['filename'] = "{}/{}_fig({})_sub({})".format(meta['plotdir'], meta['plotprefix'],
                                                       get_specvals(db, figspec, figvals),
                                                       get_specvals(db, subspec))
    return figs



def plot_opdm_fig_sub_line(db, meta, figspec, subspec, linspec, algo_filter=None, state_filter=None, figs=None,
                           palette_name=None,dbidx=0,dbnum=1):
    if db['version'] == 3:
        return plot_v3_opdm_fig3_sub3_line1(db=db, meta=meta, figspec=figspec, subspec=subspec, linspec=linspec,
                                            algo_filter=algo_filter, state_filter=state_filter, figs=figs,
                                            palette_name=palette_name)
    else:
        raise NotImplementedError('database version not implemented:' + db['version'])
