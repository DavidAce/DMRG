from itertools import product
from pathlib import Path
import matplotlib.patheffects as pe
import seaborn as sns
from scipy.optimize import curve_fit
import logging
from .tools import *

logger = logging.getLogger('plot-linearFit')


def plot_v3_slope_fig_sub_line(db, meta, figspec, subspec, linspec, xaxspec, algo_filter=None, state_filter=None,
                               point_filter=None, figs=None, palette_name=None):
    if db['version'] != 3:
        raise ValueError("plot_v3_time_fig3_sub3_line1 requires db version 3")

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
        # path_effects = None
        path_effects = [pe.SimpleLineShadow(offset=(0.5, -0.5), alpha=0.3), pe.Normal()]

    prb_style = 'prb' in meta['mplstyle'] if 'mplstyle' in meta else False

    # legend_col_keys = list(itertools.chain(l1, [col for col in meta['legendcols'] if 'legendcols' in meta]))
    legend_col_keys = linspec.copy()
    if legendcols := meta['legendcols']:
        for col in legendcols:
            if not col in [l.split(':')[0] for l in figspec + subspec + linspec + xaxspec]:
                legend_col_keys.append(col)

    figprod = list(product(*get_vals(db, figspec)))  # All combinations of figspecs values
    subprod = list(product(*get_vals(db, subspec)))  # All combinations of subspecs values
    linprod = list(product(*get_vals(db, linspec)))  # All combinations of linspecs values
    xaxprod = list(product(*get_vals(db, xaxspec)))  # All combinations of linspecs values
    # dirprod = list(product(db['keys']['algo'], db['keys']['state'], db['keys']['crono']))
    # print(dirprod)
    numfigs = len(figprod)
    numsubs = len(subprod)
    if figs is None:
        figs = [get_fig_meta(numsubs, meta=meta) for _ in range(numfigs)]

    for figvals, f in zip(figprod, figs):
        logger.debug('- plotting figs: {}'.format(figvals))
        dbval = None
        for idx, (subvals, ax, ix) in enumerate(zip(subprod, f['ax'], f['ix'])):
            popt = None
            pcov = None
            logger.debug('-- plotting subs: {}'.format(subvals))
            # for dirvals in dirprod:
            # palette = plt.rcParams['axes.prop_cycle'].by_key()['color']
            palette, lstyles = get_colored_lstyles(db, linspec, palette_name)
            for linvals, color, lstyle in zip(linprod, palette, lstyles):
                logger.debug('--- plotting lins: {}'.format(linvals))
                xvals, yvals, evals = [], [], []
                legendrow = None
                for xaxvals in xaxprod:
                    logger.debug('--- plotting xaxs: {}'.format(xaxvals))
                    datanodes = match_datanodes(db=db, meta=meta, specs=figspec + subspec + linspec + xaxspec,
                                                vals=figvals + subvals + linvals + xaxvals)
                    logger.debug('Found {} datanodes'.format(len(datanodes)))
                    for datanode in datanodes:
                        dbval = db['dsets'][datanode.name]
                        ydata = datanode['avg'][meta['colname']][()]
                        edata = datanode['std'][meta['colname']][()]  # Use standard deviation for curve_fit
                        tdata = datanode['avg']['physical_time'][()]
                        ndata = datanode['avg']['num'][()]

                        if np.min(ndata) < 10:
                            continue
                        if len(tdata) <= 1:
                            continue
                        if meta.get('findloglogwindow') is None:
                            raise ValueError('meta["findloglogwindow"] must be True')
                        if not 'entanglement_entropy' in datanode['avg'].dtype.fields:
                            raise LookupError('Could not find column "entanglement_entropy" in datanode["avg"]')

                        sdata = datanode['avg']['entanglement_entropy'][()]
                        idx1, idx2 = find_loglog_window2(tdata, sdata, dbval)

                        try:
                            if idx2 <= idx1:
                                raise IndexError("Invalid index order: idx1 {} | idx2 {}".format(idx1, idx2))
                            if idx1 + 10 > idx2:
                                raise IndexError("Too few datapoints for a fit: idx1 {} | idx2 {}".format(idx1, idx2))

                            bounds_v2 = ([-np.inf, 0], [np.inf, np.inf])
                            popt, pcov = curve_fit(f=flinear, xdata=np.log(np.log(tdata[idx1:idx2])),
                                                   ydata=ydata[idx1:idx2],
                                                   sigma=edata[idx1:idx2], bounds=bounds_v2)
                            pstd = np.sqrt(np.diag(pcov))
                            xvals.append(xaxvals)
                            yvals.append(popt[1])
                            evals.append(pstd[1])
                        except IndexError as e:
                            logger.error("Fit failed: {}".format(e))
                            pass
                        except ValueError as e:
                            logger.error("Fit failed: {}".format(e))
                            pass

                        if legendrow is None:
                            legendrow = get_legend_row(db=db, datanode=datanode, legend_col_keys=legend_col_keys)
                if legendrow is not None:
                    line = ax.errorbar(x=xvals, y=yvals, yerr=evals, color=color, path_effects=path_effects)
                    for icol, (col, key) in enumerate(zip(legendrow, legend_col_keys)):
                        key, fmt = key.split(':') if ':' in key else [key, '']
                        f['legends'][idx][icol]['handle'].append(line)
                        f['legends'][idx][icol]['label'].append(col)
                        f['legends'][idx][icol]['title'] = db['tex'][key]
                        f['legends'][idx][icol]['header'] = get_title(dbval, subspec, width=16)
                    if not idx in f['axes_used']:
                        f['axes_used'].append(idx)
            if dbval:
                ax.set_title(get_title(dbval, subspec, width=16),
                             horizontalalignment='left', x=0.05,
                             fontstretch="ultra-condensed",
                             # bbox=dict(boxstyle='square,pad=0.15', facecolor='white', alpha=0.6)
                             )
                # ax.set_xlabel(dbval['tex']['keys'][xaxspec[0].split(':')[0]])

                ax.set_xlabel(get_tex(dbval, xaxspec))
        if not prb_style and dbval:
            f['fig'].suptitle('{}\n{}'.format(meta['titlename'], get_title(dbval, figspec)))

        # prettify_plot4(fmeta=f, lgnd_meta=axes_legends)
        suffix = ''
        suffix = suffix + '_normpage' if 'normpage' in meta and meta['normpage'] else suffix
        f['filename'] = "{}/{}-slope_fig({})_sub({}){}".format(meta['plotdir'], meta['plotprefix'],
                                                       get_specvals(db, figspec, figvals),
                                                       get_specvals(db, subspec),suffix)
    return figs


def plot_v2_slope_fig3_sub3_line1(db, meta, figspec, subspec, linspec, xaxspec, algo_filter=None, state_filter=None,
                                  point_filter=None, figs=None, palette_name=None):
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
        # path_effects = None
        path_effects = [pe.SimpleLineShadow(offset=(0.5, -0.5), alpha=0.3), pe.Normal()]

    prb_style = 'prb' in meta['mplstyle'] if 'mplstyle' in meta else False

    # legend_col_keys = list(itertools.chain(l1, [col for col in meta['legendcols'] if 'legendcols' in meta]))
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
    dirprod = list(product(db['keys']['algo'], db['keys']['state'], db['keys']['crono']))
    numfigs = len(figprod)
    numsubs = len(subprod)
    if figs is None:
        figs = [get_fig_meta(numsubs, meta=meta) for _ in range(numfigs)]

    for figkeys, f in zip(figprod, figs):
        logger.debug(f'-- plotting {figkeys=}')
        dbval = None
        for idx, (subkeys, ax, ix) in enumerate(zip(subprod, f['ax'], f['ix'])):
            popt = None
            pcov = None
            logger.debug(f'-- plotting {subkeys=}')
            for dirkeys in dirprod:
                # palette = plt.rcParams['axes.prop_cycle'].by_key()['color']
                palette, lstyles = get_colored_lstyles(db, linspec, palette_name)
                for linkeys, color, lstyle in zip(linprod, palette, lstyles):
                    findlist = list(figkeys) + list(subkeys) + list(linkeys) + list(dirkeys) + [meta['groupname']]
                    datanode = [value['node']['data'] for key, value in db['dsets'].items() if
                                all(k in key for k in findlist)]
                    if len(datanode) != 1:
                        logger.warning(f"found {len(datanode)} datanodes: {datanode=} | {findlist=}")
                        continue
                        # raise LookupError("Found incorrect number of datanodes")
                    datanode = datanode[0]
                    dbval = db['dsets'][datanode.name]
                    ydata, colnames = get_table_data(datanode['avg'], meta['colname'],
                                                     'f8')  # Supports multiple columns
                    edata, colnames = get_table_data(datanode['ste'], meta['colname'],
                                                     'f8')  # Supports multiple columns
                    tdata = datanode['avg']['physical_time'][()]
                    ndata = datanode['avg']['num'][()]

                    if np.min(ndata) < 10:
                        continue
                    if len(tdata) <= 1:
                        logger.warning(f"too short: {tdata=}")
                        continue
                        # raise AssertionError("tdata is too short");
                    if meta.get('ydiff') == True:
                        ydata = np.diff(ydata.T, n=1, prepend=0).T

                    if meta.get('normpage') == True:
                        for i, (y, e) in enumerate(zip(ydata, edata)):
                            p = page_entropy(dbval['L'])
                            ydata[i] = y / p
                            edata[i] = e / p
                    if 'normalize' in meta:
                        for i, (y, e) in enumerate(zip(ydata, edata)):
                            ydata[i] = y / meta['normalize']
                            edata[i] = e / meta['normalize']

                    for i, (y, e, colname) in enumerate(zip(ydata.T, edata.T, colnames)):
                        linestyle = meta['linestyle'][i] if 'linestyle' in meta and len(
                            meta['linestyle']) == len(ydata) else '-'

                        label = None
                        if 'entanglement' in colname or 'Neumann' in colname:
                            label = '$S_\mathrm{vN}$'
                        elif 'number' in colname:
                            label = '$S_\mathrm{N}$'

                        if meta.get('plotsatapproach'):
                            sdata = datanode['avg']['entanglement_entropy'][()]
                            idx1, idx2 = find_loglog_window2(tdata, sdata, dbval)
                            ysat = np.mean(ydata[idx2:])  # Saturation value
                            y = np.abs(y - ysat)
                            try:
                                if idx2 <= idx1:
                                    raise IndexError("Invalid index order: idx1 {} | idx2 {}".format(idx1, idx2))
                                if idx1 + 10 > idx2:
                                    raise IndexError(
                                        "Too few datapoints for a fit: idx1 {} | idx2 {}".format(idx1, idx2))
                                bounds = ([-np.inf, -np.inf], [np.inf, 0])
                                with np.errstate(invalid='ignore'):
                                    try:
                                        ylog = np.log10(y)
                                        tlog = np.log10(tdata)
                                        popt, pcov = curve_fit(f=flinear, xdata=tlog[idx1:idx2],
                                                               ydata=ylog[idx1:idx2], bounds=bounds)
                                        yfit = 10 ** flinear(tlog, *popt)
                                        ax.plot(xdata, yfit, marker=None, linewidth=0.8,
                                                linestyle='--', label='fit', color=color,
                                                path_effects=path_effects)
                                        sep = 0.05 * len(f['legends'][idx][0]['handle'])
                                        xmid = 10 ** ((0.4 + sep) * (np.log10(xdata[idx2]) + np.log10(xdata[idx1])))
                                        ymid = 10 ** ((0.4 + sep) * (np.log10(yfit[idx2]) + np.log10(yfit[idx1])))
                                        xtxt = 10 ** ((0.75 + sep) * (np.log10(xdata[idx2]) + np.log10(xdata[idx1])))
                                        ytxt = 10 ** ((0.35 + sep) * (np.log10(yfit[idx2]) + np.log10(yfit[idx1])))
                                        ax.annotate('$\sim t^{{{:.2f}}}$'.format(popt[1]), xy=(xmid, ymid),
                                                    xytext=(xtxt, ytxt),
                                                    arrowprops=dict(arrowstyle="->", color=color))
                                    except:
                                        pass

                            except IndexError as e:
                                pass


                        elif meta.get('fillerror'):
                            ax.fill_between(x=xdata, y1=y - e, y2=y + e, alpha=0.10, label=None, color=color)

                        line, = ax.plot(xdata, y, marker=None, linestyle=linestyle, label=label, color=color,
                                        path_effects=path_effects)

                        if i == 0:
                            legendrow = get_legend_row(db=db, datanode=datanode, legend_col_keys=legend_col_keys)
                            for icol, (col, key) in enumerate(zip(legendrow, legend_col_keys)):
                                key, fmt = key.split(':') if ':' in key else [key, '']

                                f['legends'][idx][icol]['handle'].append(line)
                                f['legends'][idx][icol]['label'].append(col)
                                f['legends'][idx][icol]['title'] = db['tex'][key]
                                f['legends'][idx][icol]['header'] = get_title(dbval, subspec, width=16)

                        if meta.get('findloglogwindow') and 'entanglement_entropy' in datanode['avg'].dtype.fields:
                            sdata = datanode['avg']['entanglement_entropy'][()]
                            idx1, idx2 = find_loglog_window2(tdata, sdata, dbval)
                            f['ymax'] = np.max([f['ymax'], np.max(y)]) if f['ymax'] else np.max(y)
                            f['ymin'] = np.min([f['ymin'], y[idx1]]) if f['ymin'] else y[idx1]
                            if meta.get('markloglogwindow'):
                                mark, = ax.plot([xdata[idx1], xdata[idx2]], [y[idx1], y[idx2]],
                                                color=color,
                                                marker='o', markersize=6, linestyle='None', markeredgecolor='w',
                                                path_effects=path_effects)

                            if meta.get('fitloglogwindow') and not meta.get('zoomloglogwindow'):
                                # bounds = ([0., 0., 1.], [np.inf, np.inf, np.exp(1)])
                                try:
                                    # popt, pcov = curve_fit(f=floglog, xdata=tdata[idx1:idx2], ydata=y[idx1:idx2], bounds=bounds)
                                    if idx2 <= idx1:
                                        raise IndexError("Invalid index order: idx1 {} | idx2 {}".format(idx1, idx2))
                                    if idx1 + 10 > idx2:
                                        raise IndexError(
                                            "Too few datapoints for a fit: idx1 {} | idx2 {}".format(idx1, idx2))

                                    bounds_v2 = ([-np.inf, 0], [np.inf, np.inf])
                                    with np.errstate(invalid='ignore'):
                                        tloglog = np.log(np.log(tdata))
                                        popt, pcov = curve_fit(f=flinear, xdata=tloglog[idx1:idx2], ydata=y[idx1:idx2],
                                                               bounds=bounds_v2)
                                        idx_min = int(idx1 * 0.75)
                                        idx_max = int(np.min([len(tdata) - 1, idx2 * 1.25]))
                                        tfit = xdata[idx_min:idx_max]
                                        yfit = flinear(tloglog[idx_min:idx_max], *popt)
                                        ax.plot(tfit, yfit, marker=None, linewidth=0.8, alpha=1.0, zorder=50,
                                                linestyle='--', label='fit', color=color,
                                                path_effects=path_effects)
                                except IndexError as err:
                                    pass
                                except ValueError as err:
                                    logger.error(f"Fit failed: {err}")

                            if meta.get('zoomloglogwindow') and i in meta['zoomloglogwindow']['colnum']:
                                if ix is None:
                                    # pos tells where to put the inset, x0,y0, width, height in % units
                                    ix = ax.inset_axes(meta['zoomloglogwindow']['pos'])
                                ix.fill_between(x=xdata, y1=y - e, y2=y + e, alpha=0.15, label=None,
                                                color=color)
                                ix.plot(xdata, y, marker=None, linestyle=linestyle,
                                        label=None, color=color, path_effects=path_effects)
                                x1, x2, y1, y2 = meta['zoomloglogwindow']['coords']
                                x1 = np.min([x1, xdata[idx1]]) if x1 else xdata[idx1]
                                x2 = np.max([x2, xdata[idx2]]) if x2 else xdata[idx2]
                                y1 = np.min([y1, y[idx1]]) if y1 else y[idx1]
                                y2 = np.max([y2, y[idx2]]) if y2 else y[idx2]
                                meta['zoomloglogwindow']['coords'] = [x1, x2, y1, y2]
                                if meta.get('fitloglogwindow'):
                                    if idx + 5 < idx2:
                                        try:
                                            if idx2 <= idx1:
                                                raise RuntimeError("Invalid values: idx1 {} | idx2 {}", idx1,
                                                                   idx2)
                                            bounds_v2 = ([-np.inf, 0], [np.inf, np.inf])
                                            with np.errstate(invalid='ignore'):
                                                tloglog = np.log(np.log(tdata))
                                                popt, pcov = curve_fit(f=flinear, xdata=tloglog[idx1:idx2],
                                                                       ydata=y[idx1:idx2], bounds=bounds_v2)
                                                ix.plot(xdata, flinear(tloglog, *popt), marker=None,
                                                        linewidth=0.8,
                                                        linestyle='--', label=None, color=color,
                                                        path_effects=path_effects)
                                        except RuntimeError as e:
                                            pass
                                if meta.get('markloglogwindow'):
                                    ix.plot([xdata[idx1], xdata[idx2]], [y[idx1], y[idx2]],
                                            color=color,
                                            marker='o', markersize=6, linestyle='None', markeredgecolor='w',
                                            path_effects=path_effects)

                    if not idx in f['axes_used']:
                        f['axes_used'].append(idx)
            if dbval:
                ax.set_title(get_title(dbval, subspec, width=16),
                             horizontalalignment='left', x=0.05,
                             fontstretch="ultra-condensed",
                             # bbox=dict(boxstyle='square,pad=0.15', facecolor='white', alpha=0.6)
                             )

            ax.set_xlabel("$t$")
            if meta.get('timeselection') == 'lnt':
                ax.set_xscale('log')
                ymin = None
                ymax = None
                if ix is not None:
                    ix.set_xscale('log')

            elif meta.get('timeselection') == 'lnlnt':
                ax.set_xlabel("$\ln\ln t$")

            if meta.get('zoomloglogwindow') and ix is not None:
                x1, x2, y1, y2 = meta['zoomloglogwindow']['coords']  # sub region of the original image
                ix.set_xlim(xmin=0.1 * x1, xmax=10 * x2)
                ix.set_ylim(ymin=0.975 * y1, ymax=1.025 * y2)
                # ix.set_xticklabels('')
                # ix.set_yticklabels('')
                ix.tick_params(axis='both', which='both', labelsize='x-small')
                # ix.xaxis.set_major_locator(plt.MaxNLocator(5))
                ix.xaxis.set_major_locator(plt.LogLocator(base=10, numticks=6))
                if meta.get('timeselection') == 'lnlnt':
                    ix.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
                    ix.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))

                ix.yaxis.tick_right()
                ix.patch.set_linewidth('2')
                ix.patch.set_edgecolor('black')
                ix.legend(title=meta['zoomloglogwindow']['legendtitle'], fontsize='x-small',
                          loc='lower center', framealpha=0.9,
                          bbox_to_anchor=(0.49, -0.01),
                          bbox_transform=ix.transAxes,
                          borderaxespad=0.5
                          )
                ax.indicate_inset_zoom(ix, edgecolor="black")

        if f['ymin']:
            f['ymin'] = 0.9 * f['ymin']
        if f['ymax']:
            f['ymax'] = 1.1 * f['ymax']

        if not prb_style and dbval:
            f['fig'].suptitle('{} vs Time\n{}'.format(meta['titlename'], get_title(dbval, figspec)))

        # prettify_plot4(fmeta=f, lgnd_meta=axes_legends)
        suffix = ''
        suffix = suffix + '_normpage' if 'normpage' in meta and meta['normpage'] else suffix
        suffix = suffix + '_loglog' if meta.get('timeselection') == 'lnlnt' else suffix
        f['filename'] = "{}/{}(t)_fig({})_sub({}){}".format(meta['plotdir'], meta['plotprefix'],
                                                            '-'.join(map(str, figkeys)),
                                                            '-'.join(map(str, get_keys(db, subspec))),
                                                            suffix)

    return figs


def plot_slope_fig_sub_line(db, meta, figspec, subspec, linspec, xaxspec, algo_filter=None, state_filter=None,
                            point_filter=None, figs=None, palette_name=None):
    if db['version'] == 2:
        specs = figspec + subspec + linspec
        nonv3spec = lambda x: x not in ['cstd', 'tstd', 'tgw8', 'cgw8']
        specs = list(filter(nonv3spec, specs))
        fig3 = specs[:3]
        sub2 = specs[3:5]
        lin1 = specs[5:6]
        xax1 = specs[6:7]
        return plot_v2_slope_fig3_sub3_line1(db=db, meta=meta, figspec=fig3, subspec=sub2, linspec=lin1, xaxspec=xax1,
                                             algo_filter=algo_filter, state_filter=state_filter,
                                             point_filter=point_filter, figs=figs, palette_name=palette_name)
    elif db['version'] == 3:
        return plot_v3_slope_fig_sub_line(db=db, meta=meta, figspec=figspec, subspec=subspec, linspec=linspec,
                                          xaxspec=xaxspec,
                                          algo_filter=algo_filter, state_filter=state_filter, point_filter=point_filter,
                                          figs=figs, palette_name=palette_name)
    else:
        raise NotImplementedError('database version not implemented:' + db['version'])
