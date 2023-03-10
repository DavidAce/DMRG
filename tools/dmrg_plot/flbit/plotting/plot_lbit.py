from itertools import product
from pathlib import Path
import seaborn as sns
from .tools import *
import warnings
import logging
import matplotlib.patheffects as pe

logger = logging.getLogger('plot-lbit')


def plot_v3_lbit_fig3_sub3_line1(db, meta, figspec, subspec, linspec, algo_filter=None, state_filter=None, figs=None,
                                 palette_name=None):
    if db['version'] != 3:
        raise ValueError("plot_v3_lbit_fig3_sub3_line1 requires db version 3")
    if 'mplstyle' in meta:
        plt.style.use(meta['mplstyle'])
    if 'mplstyle' in meta and 'slack' in meta['mplstyle']:
        # palette_name = "Spectral"
        if not palette_name:
            palette_name = "colorblind"
        # path_effects = None
        path_effects = [pe.SimpleLineShadow(offset=(0.5, -0.5), alpha=0.3), pe.Normal()]
    else:
        if not palette_name:
            palette_name = "colorblind"
        # path_effects = None
        path_effects = [pe.SimpleLineShadow(offset=(-0.5, -0.5), alpha=0.3), pe.Normal()]

    prb_style = 'prb' in meta['mplstyle'] if 'mplstyle' in meta else False

    legend_col_keys = linspec.copy()
    if legendcols := meta['legendcols']:
        for col in legendcols:
            if not col in [l.split(':')[0] for l in subspec + linspec]:
                legend_col_keys.append(col)

    figprod = list(product(*get_vals(db, figspec)))  # All combinations of figspecs values
    subprod = list(product(*get_vals(db, subspec)))  # All combinations of subspecs values
    linprod = list(product(*get_vals(db, linspec)))  # All combinations of linspecs values
    numfigs = len(figprod)
    numsubs = len(subprod)
    if figs is None:
        figs = [get_fig_meta(numsubs, meta=meta) for _ in range(numfigs)]

    for figvals, f in zip(figprod, figs):
        logger.debug('- plotting figkeys: {}'.format(figvals))
        dbval = None
        for idx, (subvals, ax, ix) in enumerate(zip(subprod, f['ax'], f['ix'])):
            logger.debug('-- plotting subkeys: {}'.format(subvals))
            palette, lstyles = get_colored_lstyles(db, linspec, palette_name)
            for linvals, color, lstyle in zip(linprod, palette, lstyles):
                logger.debug('--- plotting linkeys: {}'.format(linvals))
                datanodes = match_datanodes(db=db, meta=meta, specs=figspec + subspec + linspec,
                                            vals=figvals + subvals + linvals)
                logger.debug('Found {} datanodes'.format(len(datanodes)))
                for datanode in datanodes:
                    # datanode = datanode[0]
                    dbval = db['dsets'][datanode.name]
                    ndata, _ = get_table_data(datanode['num'], meta.get('colname'), 'i8')
                    ydata, _ = get_table_data(datanode['avg'])
                    edata, _ = get_table_data(datanode['ste'])
                    xdata = np.array(range(len(ydata)))
                    fdata = get_lbit_cls(xdata, ydata)  # Fit to get characteristic length-scale
                    if 'yscale' in meta and meta['yscale'] == 'log':
                        ydata = np.abs(ydata)
                    if np.min(ndata) < 10:
                        continue
                    for i, (y, e) in enumerate(zip(ydata.T, edata.T)):
                        # ax.fill_between(x=xdata, y1=y - e, y2=y + e, alpha=0.10, color=color)
                        line, = ax.plot(xdata, y, marker=None, color=color, path_effects=path_effects)
                        ax.scatter(xdata[1], y[1], marker='o', color=color, path_effects=path_effects)
                        ax.scatter(xdata[fdata[-1]], y[fdata[-1]], marker='o', color=color, path_effects=path_effects,
                                   zorder=50)
                        if i == 0:
                            legendrow = get_legend_row(db=db, datanode=datanode, legend_col_keys=legend_col_keys)
                            for icol, (col, key) in enumerate(zip(legendrow, legend_col_keys)):
                                key, fmt = key.split(':') if ':' in key else [key, '']
                                f['legends'][idx][icol]['handle'].append(line)
                                f['legends'][idx][icol]['label'].append(col)
                                f['legends'][idx][icol]['title'] = db['tex'][key]
                            icol = len(legendrow)
                            iuse = 0
                            if 'C' in meta['legendfits']:
                                f['legends'][idx][icol + iuse]['handle'].append(line)
                                f['legends'][idx][icol + iuse]['label'].append('{:.2f}'.format(fdata[0]))
                                f['legends'][idx][icol + iuse]['title'] = '$C$'
                                iuse += 1
                            if 'xi' in meta['legendfits']:
                                f['legends'][idx][icol + iuse]['handle'].append(line)
                                f['legends'][idx][icol + iuse]['label'].append('{:.2f}'.format(fdata[1]))
                                f['legends'][idx][icol + iuse]['title'] = '$\\xi_\\tau$'
                                iuse += 1
                            if 'beta' in meta['legendfits'] and fdata[2] is not None:
                                f['legends'][idx][icol + iuse]['handle'].append(line)
                                f['legends'][idx][icol + iuse]['label'].append('{:.2f}'.format(fdata[2]))
                                f['legends'][idx][icol + iuse]['title'] = '$\\beta$'

                        line, = ax.plot(xdata, regular_exp(xdata, fdata[0], fdata[1]), linewidth=1.0, marker=None,
                                        linestyle='dashed',
                                        color=color, path_effects=path_effects, zorder=0)
                        # if 'inset-cls' in meta:
                        #     if ix is None:
                        #         # pos tells where to put the inset, x0,y0, width, height in % units
                        #         ix = ax.inset_axes(meta['inset-cls']['pos'])
                        #         ix.set_ylabel('$\\xi_\\tau$')
                        #         ix.set_xlabel('$f$')
                        #     ix.scatter(dbval['vals']['f'], fdata[1], color=color)

            if not idx in f['axes_used']:
                f['axes_used'].append(idx)

            if dbval:
                ax.set_title(get_title(dbval, subspec),
                             horizontalalignment='left', x=0.05,
                             fontstretch="ultra-condensed",
                             # bbox=dict(boxstyle='square,pad=0.15', facecolor='white', alpha=0.6)
                             )

        if not prb_style and dbval:
            f['fig'].suptitle('{}\n{}'.format(meta['titlename'], get_title(dbval, figspec)))

        f['filename'] = "{}/{}_lbit_fig({})_sub({})".format(meta['plotdir'], meta['plotprefix'],
                                                            '-'.join(map(str, figvals)),
                                                            '-'.join(map(str, get_keys(db, subspec))))

    return figs


def plot_v2_lbit_fig3_sub3_line1(db, meta, figspec, subspec, linspec, algo_filter=None, state_filter=None, figs=None,
                                 palette_name=None):
    if db['version'] != 2:
        raise ValueError("plot_v2_lbit_fig3_sub3_line1 requires db version 2")
    if len(figspec) + len(subspec) + len(linspec) != 7:
        warnings.warn("plot_v2_lbit_fig3_sub3_line1: specs should add to 7")
        # raise AssertionError("Must add to 7 elems: \n figspec {}\n subspec {}\n linespec {}")

    if 'mplstyle' in meta:
        plt.style.use(meta['mplstyle'])
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

    legend_col_keys = linspec.copy()
    if legendcols := meta['legendcols']:
        for col in legendcols:
            if not col in [l.split(':')[0] for l in subspec + linspec]:
                legend_col_keys.append(col)

    figprod = list(product(*get_keys(db, figspec)))  # All combinations of figspecs (names of parameters that iterate figures)
    subprod = list(product(*get_keys(db, subspec)))  # All combinations of subspecs (names of parameters that iterate subplots)
    linprod = list(product(*get_keys(db, linspec)))  # All combinations of linspecs (names of parameters that iterate lines)
    dirprod = list(product(db['keys']['algo'], db['keys']['model']))
    numfigs = len(figprod)
    numsubs = len(subprod)
    if figs is None:
        figs = [get_fig_meta(numsubs, meta=meta) for _ in range(numfigs)]

    for figkeys, f in zip(figprod, figs):
        logger.debug('- plotting figkeys: {}'.format(figkeys))
        dbval = None
        for idx, (subkeys, ax, ix) in enumerate(zip(subprod, f['ax'], f['ix'])):
            logger.debug('-- plotting subkeys: {}'.format(subkeys))
            for dirkeys in dirprod:
                palette, lstyles = get_colored_lstyles(db, linspec, palette_name)
                for linkeys, color, lstyle in zip(linprod, palette, lstyles):
                    findlist = list(figkeys) + list(subkeys) + list(linkeys) + list(dirkeys) + [meta['groupname']]
                    datanode = [value['node']['data'] for key, value in db['dsets'].items() if
                                all(k in key for k in findlist)]
                    if len(datanode) != 1:
                        logger.error("found", len(datanode), "datanodes: ", datanode, " | findlist: ", findlist)
                        continue
                        # raise LookupError("Found incorrect number of datanodes")
                    datanode = datanode[0]
                    dbval = db['dsets'][datanode.name]
                    ydata, _ = get_table_data(datanode['avg'])
                    edata, _ = get_table_data(datanode['ste'])
                    ndata = datanode['num'][()]
                    xdata = np.array(range(len(ydata)))
                    C, xi, beta, yfit, res = get_lbit_cls(np.arange(0, len(ydata)),
                                                          ydata)  # Fit to get characteristic length-scale

                    if np.min(ndata) < 10:
                        continue
                    for i, (y, e) in enumerate(zip(ydata.T, edata.T)):
                        # ax.fill_between(x=xdata, y1=y - e, y2=y + e, alpha=0.10, color=color)
                        line, = ax.plot(xdata, y, marker=None, color=color, path_effects=path_effects)
                        if i == 0:
                            legendrow = get_legend_row(db=db, datanode=datanode, legend_col_keys=legend_col_keys)
                            for icol, (col, key) in enumerate(zip(legendrow, legend_col_keys)):
                                key, fmt = key.split(':') if ':' in key else [key, '']
                                f['legends'][idx][icol]['handle'].append(line)
                                f['legends'][idx][icol]['label'].append(col)
                                f['legends'][idx][icol]['title'] = db['tex'][key]
                            icol = len(legendrow)
                            iuse = 0
                            if 'C' in meta['legendfits']:
                                f['legends'][idx][icol + iuse]['handle'].append(line)
                                f['legends'][idx][icol + iuse]['label'].append('{:.2f}'.format(C))
                                f['legends'][idx][icol + iuse]['title'] = '$C$'
                                iuse += 1
                            if 'xi' in meta['legendfits']:
                                f['legends'][idx][icol + iuse]['handle'].append(line)
                                f['legends'][idx][icol + iuse]['label'].append('{:.2f}'.format(xi))
                                f['legends'][idx][icol + iuse]['title'] = '$\\xi_\\tau$'
                                iuse += 1
                            if 'beta' in meta['legendfits'] and beta is not None:
                                f['legends'][idx][icol + iuse]['handle'].append(line)
                                f['legends'][idx][icol + iuse]['label'].append('{:.2f}'.format(beta))
                                f['legends'][idx][icol + iuse]['title'] = '$\\beta$'

                        line, = ax.plot(xdata, yfit, marker=None, linestyle='dashed',
                                        color=color, path_effects=path_effects)
                        if 'inset-cls' in meta:
                            if ix is None:
                                # pos tells where to put the inset, x0,y0, width, height in % units
                                ix = ax.inset_axes(meta['inset-cls']['pos'])
                                ix.set_ylabel('$\\xi_\\tau$')
                                ix.set_xlabel('$f$')
                            ix.scatter(dbval['vals']['f'], fdata[1], color=color)

                if not idx in f['axes_used']:
                    f['axes_used'].append(idx)

            if dbval:
                ax.set_title(get_title(dbval, subspec),
                             horizontalalignment='left', x=0.05,
                             fontstretch="ultra-condensed",
                             # bbox=dict(boxstyle='square,pad=0.15', facecolor='white', alpha=0.6)
                             )

        if not prb_style and dbval:
            f['fig'].suptitle('{}\n{}'.format(meta['titlename'], get_title(dbval, figspec)))

        f['filename'] = "{}/{}_lbit_fig({})_sub({})".format(meta['plotdir'], meta['plotprefix'],
                                                            '-'.join(map(str, figkeys)),
                                                            '-'.join(map(str, get_keys(db, subspec))))

    return figs


def plot_lbit_fig_sub_line(db, meta, figspec, subspec, linspec, algo_filter=None, state_filter=None, figs=None,
                           palette_name=None):
    if db['version'] == 2:
        specs = figspec + subspec + linspec
        is_v3spec = lambda x: x in ['cstd', 'tstd', 'tgw8', 'cgw8']
        specs = list(filter(is_v3spec, specs))
        fig3 = specs[:3]
        sub3 = specs[3:6]
        lin1 = specs[6:7]
        return plot_v2_lbit_fig3_sub3_line1(db=db, meta=meta, figspec=fig3, subspec=sub3, linspec=lin1,
                                            algo_filter=algo_filter, state_filter=state_filter, figs=figs,
                                            palette_name=palette_name)
    elif db['version'] == 3:
        return plot_v3_lbit_fig3_sub3_line1(db=db, meta=meta, figspec=figspec, subspec=subspec, linspec=linspec,
                                            algo_filter=algo_filter, state_filter=state_filter, figs=figs,
                                            palette_name=palette_name)
    else:
        raise NotImplementedError('database version not implemented:' + db['version'])
