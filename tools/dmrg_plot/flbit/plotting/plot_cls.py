from itertools import product
from pathlib import Path
import seaborn as sns
import matplotlib.patheffects as pe
from .tools import *


def plot_v3_cls_fig3_sub3_line1(db, meta, figspec, subspec, linspec, xaxspec, algo_filter=None, state_filter=None,
                                point_filter=None, figs=None, palette_name=None):
    if db['version'] != 3:
        raise ValueError("plot_v3_dset_fig3_sub3_line1 requires db version 3")

    if 'mplstyle' in meta:
        plt.style.use(meta['mplstyle'])
    if 'plotdir' in meta and 'mplstyle' in meta:
        if Path(meta['plotdir']).stem != Path(meta['mplstyle']).stem:
            meta['plotdir'] = Path(meta['plotdir'], Path(meta['mplstyle']).stem)
            Path(meta['plotdir']).mkdir(parents=True, exist_ok=True)
    path_effects = [pe.SimpleLineShadow(offset=(0.5, -0.5), alpha=0.3), pe.Normal()]
    prb_style = 'prb' in meta['mplstyle'] if 'mplstyle' in meta else False

    # legend_col_keys = list(itertools.chain(l1, [col for col in meta['legendcols'] if 'legendcols' in meta]))
    legend_col_keys = linspec.copy()
    if legendcols := meta['legendcols']:
        for col in legendcols:
            if not col in [l.split(':')[0] for l in figspec + subspec + linspec + xaxspec]:
                legend_col_keys.append(col)

    figprod = list(product(*get_vals(db=db, keyfmt=figspec, filter=meta.get('filter'))))  # All combinations of figspecs values
    subprod = list(product(*get_vals(db=db, keyfmt=subspec, filter=meta.get('filter'))))  # All combinations of subspecs values
    linprod = list(product(*get_vals(db=db, keyfmt=linspec, filter=meta.get('filter'))))  # All combinations of linspecs values
    xaxprod = list(product(*get_vals(db=db, keyfmt=xaxspec, filter=meta.get('filter'))))  # All combinations of linspecs values
    # dirprod = list(product(db['keys']['algo'], db['keys']['state'], db['keys']['crono']))
    # print(dirprod)
    numfigs = len(figprod)
    numsubs = len(subprod)
    if figs is None:
        figs = [get_fig_meta(numsubs, meta=meta) for _ in range(numfigs)]

    for figvals, f in zip(figprod, figs):
        logger.debug('- plotting figs {}: {}'.format(figspec, figvals))
        dbval = None
        for idx, (subvals, ax, ix) in enumerate(zip(subprod, f['ax'], f['ix'])):
            logger.debug('-- plotting subs {}: {}'.format(subspec, subvals))
            palette, lstyles = get_colored_lstyles(db, linspec, palette_name)
            for linvals, color, lstyle in zip(linprod, palette, lstyles):
                logger.debug('--- plotting lins {}: {}'.format(linspec, linvals))
                xvals, yvals, evals, cvals, nvals = [], [], [], [], []
                legendrow = None
                for xaxvals in xaxprod:
                    logger.debug('--- plotting xaxs {}: {}'.format(xaxspec, xaxvals))
                    datanodes = match_datanodes(db=db, meta=meta, specs=figspec + subspec + linspec + xaxspec,
                                                vals=figvals + subvals + linvals + xaxvals)
                    if len(datanodes) != 1:
                        logger.warning(f"match: \n"
                                       f"\tspec:{[figspec + subspec + linspec]}\n"
                                       f"\tvals:{[figvals + subvals + linvals]}")
                        logger.warning(f"found {len(datanodes)} datanodes: {datanodes=}")

                    for datanode in datanodes:
                        dbval = db['dsets'][datanode.name]
                        ndata = dbval['vals']['num']
                        Ldata = dbval['vals']['L']
                        xdata = np.array(range(Ldata))
                        lbavg = get_lbit_avg(corrmat=datanode[()], site=meta.get('lbit-site'),
                                             mean=meta.get('lbit-mean'))

                        yfits = []
                        for yfull, ystdv in zip(lbavg.full.T, lbavg.stdv.T):
                            yfits.append(get_lbit_fit_data(x=xdata, y=np.atleast_2d(yfull).T, e=np.atleast_2d(ystdv).T,
                                                           beta=meta.get('fit-beta', True),
                                                           ymin=meta.get('fit-ymin', 1e-16),
                                                           ))  # Fit to get characteristic length-scale
                        print(
                            f"L{Ldata} | f{dbval['vals']['f']:.2f} | u{dbval['vals']['u']:2} | mean {meta.get('lbit-mean')} |  cls {[yfit.xi for yfit in yfits]}")
                        xvals.append(xaxvals)
                        yvals.append([yfit.xi for yfit in yfits])
                        evals.append([yfit.xierr for yfit in yfits])
                        nvals.append(ndata)
                        if legendrow is None:
                            legendrow = get_legend_row(db=db, datanode=datanode, legend_col_keys=legend_col_keys)
                xvals = np.array(xvals, ndmin=1)
                yvals = np.array(yvals, ndmin=2)
                evals = np.array(evals, ndmin=2)
                # cvals = np.array(cvals, ndmin=2).T
                if legendrow is not None:
                    for y, e in zip(yvals.T, evals.T):
                        if meta.get('fillerror'):
                            ax.fill_between(x=xvals[None], y1=y - e, y2=y + e, alpha=0.10, label=None, color=color)
                            line, = ax.plot(xvals, y, marker=None,
                                            # linestyle=linestyle,
                                            label=None, color=color, path_effects=path_effects)
                        else:
                            line = ax.errorbar(x=xvals, y=y, yerr=e, color=color, path_effects=path_effects)
                            # line = ax.errorbar(x=xvals, y=yvals, yerr=evals, color=color, path_effects=path_effects)

                        for icol, (col, key) in enumerate(zip(legendrow, legend_col_keys)):
                            key, fmt = key.split(':') if ':' in key else [key, '']
                            f['legends'][idx][icol]['handle'].append(line)
                            f['legends'][idx][icol]['label'].append(col)
                            f['legends'][idx][icol]['title'] = db['tex'][key]
                            f['legends'][idx][icol]['header'] = get_title(dbval, subspec, width=16)

                        if not idx in f['axes_used']:
                            f['axes_used'].append(idx)

            if dbval:
                ax.set_xlabel(get_tex(dbval, xaxspec))

            if axtitle := get_default(meta, 'axtitle'):
                if dbval and isinstance(axtitle, bool):
                    axtitle = get_title(dbval, subspec)
                ax.set_title(axtitle,horizontalalignment='left', x=0.05,fontstretch="ultra-condensed")

        if figspec_title := get_figspec_title(meta, dbval, figspec):
            f['fig'].suptitle(figspec_title)

        # prettify_plot4(fmeta=f, lgnd_meta=axes_legends)
        suffix = ''
        suffix = suffix + '_normpage' if 'normpage' in meta and meta['normpage'] else suffix
        suffix = suffix + '_loglog' if 'timeloglevel' in meta and meta['timeloglevel'] >= 2 else suffix
        f['filename'] = "{}/{}_fig({})_sub({}){}".format(meta['plotdir'], meta['plotprefix'],
                                                            '-'.join(map(str, figvals)),
                                                            '-'.join(map(str, get_keys(db, subspec))),
                                                            suffix)

    return figs


def plot_v2_cls_fig3_sub3_line1(db, meta, figspec, subspec, linspec, algo_filter=None, state_filter=None, figs=None,
                                palette_name=None):
    if len(figspec) + len(subspec) + len(linspec) != 7:
        raise AssertionError("Must add to 7 elems: \n figspec {}\n subspec {}\n linespec {}")

    if 'mplstyle' in meta:
        plt.style.use(meta['mplstyle'])
    if 'plotdir' in meta and 'mplstyle' in meta:
        if Path(meta['plotdir']).stem != Path(meta['mplstyle']).stem:
            meta['plotdir'] = Path(meta['plotdir'], Path(meta['mplstyle']).stem)
            Path(meta['plotdir']).mkdir(parents=True, exist_ok=True)
            print("Setting plotdir: ", meta['plotdir'])
    path_effects = [pe.SimpleLineShadow(offset=(0.5, -0.5), alpha=0.3), pe.Normal()]
    prb_style = 'prb' in meta['mplstyle'] if 'mplstyle' in meta else False

    legend_col_keys = linspec.copy()
    if legendcols := meta['legendcols']:
        for col in legendcols:
            if not col in [l.split(':')[0] for l in subspec + linspec]:
                legend_col_keys.append(col)

    figprod = list(
        product(*get_keys(db, figspec)))  # All combinations of figspecs (names of parameters that iterate figures)
    subprod = list(
        product(*get_keys(db, subspec)))  # All combinations of subspecs (names of parameters that iterate subplots)
    linprod = list(
        product(*get_keys(db, linspec)))  # All combinations of linspecs (names of parameters that iterate lines)
    dirprod = list(product(db['keys']['algo'], db['keys']['model']))
    numfigs = len(figprod)
    numsubs = len(subprod)
    if figs is None:
        figs = [get_fig_meta(numsubs, meta=meta) for _ in range(numfigs)]

    for figkeys, f in zip(figprod, figs):
        print('- plotting figkeys: {}'.format(figkeys))
        dbval = None
        for idx, (subkeys, ax) in enumerate(zip(subprod, f['ax'])):
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
                    ydata = []
                    xdata = []
                    for linkey in linkeys:
                        datanode = datanode[0]
                        dbval = db['dsets'][datanode.name]
                        adata = datanode['avg'][()]
                        popt = get_lbit_cls(adata)
                        if popt is not None:
                            ydata.append(popt[1])
                            xdata.append(dbval['vals']['f'])
                    line = ax.scatter(xdata, ydata, path_effects=path_effects)

                    # for i, (y, e) in enumerate(zip(ydata.T, edata.T)):
                    #     ax.fill_between(x=xdata, y1=y - e, y2=y + e, alpha=0.10, color=color)
                    #     line, = ax.plot(xdata, y, marker=None, color=color, path_effects=path_effects)
                    #
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

        if not f['filename']:
            f['filename'] = "{}/{}_lbit_fig({})_sub({})".format(meta['plotdir'], meta['plotprefix'],
                                                                '-'.join(map(str, figkeys)),
                                                                '-'.join(map(str, get_keys(db, subspec))))

    return figs


def plot_cls_fig_sub_line(db, meta, figspec, subspec, linspec, xaxspec, algo_filter=None, state_filter=None, figs=None,
                          palette_name=None):
    if db['version'] == 2:
        specs = figspec + subspec + linspec
        is_v3spec = lambda x: x in ['cstd', 'tstd', 'tgw8', 'cgw8']
        specs = list(filter(is_v3spec, specs))
        fig3 = specs[:3]
        sub3 = specs[3:6]
        lin1 = specs[6:7]
        xax1 = specs[7:8]
        return plot_v2_cls_fig3_sub3_line1(db=db, meta=meta, figspec=fig3, subspec=sub3, linspec=lin1,
                                           algo_filter=algo_filter, state_filter=state_filter, figs=figs,
                                           palette_name=palette_name)
    elif db['version'] == 3:
        return plot_v3_cls_fig3_sub3_line1(db=db, meta=meta, figspec=figspec, subspec=subspec, linspec=linspec,
                                           xaxspec=xaxspec,
                                           algo_filter=algo_filter, state_filter=state_filter, figs=figs,
                                           palette_name=palette_name)
    else:
        raise NotImplementedError('database version not implemented:' + db['version'])
