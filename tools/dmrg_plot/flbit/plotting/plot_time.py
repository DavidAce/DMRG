from itertools import product
from pathlib import Path
import matplotlib.patheffects as pe
import scipy.stats
import seaborn as sns
import logging
from scipy.optimize import curve_fit
from .tools import *
import matplotlib.transforms as transforms

logger = logging.getLogger('plot-time')

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

@njit(cache=True,parallel=True)
def get_particle_number_fluctuation(p):
    n_len, t_len, r_len = np.shape(p)
    n_vals = np.arange(n_len)
    n_var = np.zeros(shape=(t_len,r_len))
    for r in prange(r_len):
        for t in range(t_len):
            n_mean = np.sum(p[:, t, r] * n_vals)
            n_var[t,r] = np.sum(p[:,t,r] * (n_vals - n_mean)**2)
    return n_var

@njit(cache=True,parallel=True)
def get_renyi_number_entropy(p,alpha):
    n_len, t_len, r_len = np.shape(p)
    n_vals = np.arange(n_len)
    nrenyi = np.zeros(shape=(t_len,r_len))
    if alpha != 1:
        alpha_inv = (1.0/(1-alpha))
        for r in prange(r_len):
            for t in range(t_len):
                nrenyi[t,r] = alpha_inv * np.log(np.sum(p[:, t, r]** alpha ))
    else:
        for r in prange(r_len):
            for t in range(t_len):
                nrenyi[t,r] = - np.sum(p[:, t, r] * np.log( p[:, t, r] ))
    return nrenyi


def plot_v3_time_fig_sub_line(db, meta, figspec, subspec, linspec, algo_filter=None, state_filter=None,
                              point_filter=None, figs=None, palette_name=None, debug=False,dbidx=0, dbnum=1):
    if db['version'] != 3:
        raise ValueError("plot_v3_time_fig3_sub3_line1 requires db version 3")
    path_effects = [pe.SimpleLineShadow(offset=(0.35, -0.35), alpha=0.3), pe.Normal()]
    # path_effects = [pe.Stroke(linewidth=1.4, foreground='black'), pe.Normal()]
    # yfit_effects = [pe.SimpleLineShadow(offset=(-0.35, -0.35), alpha=0.5, shadow_color='gray'), pe.Normal()]
    # yfit_effects = [pe.Stroke(linewidth=1.2, foreground='white'), pe.Normal()]
    # yfit_effects = [pe.Stroke(linewidth=1.2, foreground='white'), pe.Normal()]
    yfit_effects = [pe.Stroke(linewidth=0.50, foreground='gray'), pe.Normal()]
    # mark_effects = [pe.Stroke(linewidth=0.90, foreground='black'), pe.Normal()]
    mark_effects = [pe.SimpleLineShadow(offset=(0.35, -0.35), alpha=0.4), pe.Normal()]

    patch_effects = [pe.SimplePatchShadow(offset=(0.35, -0.35), alpha=0.3), pe.Normal()]
    if 'mplstyle' in meta:
        plt.style.use(meta['mplstyle'])
        if 'slack' in meta['mplstyle']:
            path_effects = [pe.SimpleLineShadow(offset=(0.5, -0.5), alpha=0.3), pe.Normal()]


    # legend_col_keys = list(itertools.chain(l1, [col for col in meta['legendcols'] if 'legendcols' in meta]))
    legend_col_keys = linspec.copy()
    if legendcols := meta.get('legendcols'):
        for col in legendcols:
            if not col in [l.split(':')[0] for l in figspec + subspec + linspec]:
                legend_col_keys.append(col)

    figprod = list(product(*get_vals(db=db, keyfmt=figspec, filter=meta.get('filter'))))  # All combinations of figspecs values
    subprod = list(product(*get_vals(db=db, keyfmt=subspec, filter=meta.get('filter'))))  # All combinations of subspecs values
    linprod = list(product(*get_vals(db=db, keyfmt=linspec, filter=meta.get('filter'))))  # All combinations of linspecs values

    # dirprod = list(product(db['keys']['algo'], db['keys']['state'], db['keys']['crono']))
    numfigs = len(figprod)
    numsubs = len(subprod)

    numsubs *=dbnum
    subprod *=dbnum

    if figs is None:
        figs = [get_fig_meta(numsubs, meta=meta) for _ in range(numfigs)]

    for figvals, f in zip(figprod, figs):
        logger.debug('- plotting figs: {}'.format(figvals))
        dbval = None
        for idx, (subvals, ax, ix) in enumerate(zip(subprod, f['ax'], f['ix'])):
            idx_palette = idx % len(palette_name)
            if dbidx != int(idx/dbnum):
                continue

            popt = None
            pcov = None
            logger.debug('-- plotting subs: {}'.format(subvals))
            # for dirvals in dirprod:
            # palette = plt.rcParams['axes.prop_cycle'].by_key()['color']
            palette, lstyles = get_colored_lstyles(db, linspec, palette_name, filter=None, idx=idx_palette)
            xm1,ym1,cm1,xm2,ym2,cm2,xm3,ym3,cm3,xm4,ym4,cm4 = [],[],[],[],[],[],[],[],[],[],[],[]

            for lidx, (linvals, color, lstyle) in enumerate(zip(linprod, palette, lstyles)):
                logger.debug('--- plotting lins: {}'.format(linvals))
                datanodes = match_datanodes(db=db, meta=meta, specs=figspec + subspec + linspec,
                                            vals=figvals + subvals + linvals)
                logger.debug('Found {} datanodes'.format(len(datanodes)))
                for datanode in datanodes:
                    dbval = db['dsets'][datanode.name]

                    # if dbval['vals']['f'] < 0.4:
                    #     continue
                    ystat = meta['ystat'] if 'ystat' in meta else 'avg'
                    if ystat == 'typ' and not 'typ' in datanode:
                        data = datanode[meta.get('colname')][()]
                        ydata = np.atleast_2d(scipy.stats.mstats.gmean(np.where(data < 1e-10, 1e-10, data), axis=1)).T
                    else:
                        ydata, colnames = get_table_data(datanode[ystat], meta.get('colname'),'f8')  # Supports multiple columns
                    edata, colnames = get_table_data(datanode['ste'], meta.get('colname'),'f8')  # Supports multiple columns
                    tdata = datanode['avg']['physical_time'][()].astype(float)
                    ndata = datanode['avg']['num'][()]
                    if meta.get('plotrenyi2'):
                        # p = datanode.parent.parent['number_probabilities'][()] # indexed by n, time and realization
                        # nren2 = get_renyi_number_entropy(p,alpha=2)
                        # ydata = np.atleast_2d(np.mean(nren2, axis=1)).T
                        ydata = np.atleast_2d(datanode.parent.parent['renyi2_number_entropies/avg'][()]).T
                    if meta.get('plothartley'):
                        ydata = np.atleast_2d(datanode.parent.parent['hartley_number_entropies/avg'][()]).T
                    if meta.get('plotVarN'):
                        # Compute the particle number fluctuation from the number probabilities
                        p = datanode.parent.parent['number_probabilities'][()] # indexed by n, time and realization
                        n_var = get_particle_number_fluctuation(p)
                        vdata = np.atleast_2d(np.mean(n_var, axis=1)).T
                        ydata = 0.5*np.log(2*np.pi * np.e *(vdata + 1.0/12))

                    if np.min(ndata) < 10:
                        continue
                    if len(tdata) <= 1:
                        logger.warning("tdata is too short: ", tdata)
                        continue
                        # raise AssertionError("tdata is too short");
                    if meta.get('ydiff') == True:
                        ydata = np.diff(ydata.T, n=1, prepend=0).T

                    if meta.get('normpage') == True:
                        for i, (y, e) in enumerate(zip(ydata.T, edata.T)):
                            p = page_entropy(dbval['vals']['L'])
                            ydata.T[i] = y / p
                            edata.T[i] = e / p
                    if meta.get('normsinf') is True:
                        for i, (y, e) in enumerate(zip(ydata.T, edata.T)):
                            t = get_timepoints(tdata,dbval)
                            idx_num, idx_ent = t.idx_num_saturated, t.idx_ent_saturated
                            idx_sat = idx_num if 'number' in meta['colname'] else idx_ent
                            # sdata = datanode['avg']['entanglement_entropy'][()]
                            ytavg = np.mean(y[idx_sat:])
                            ydata.T[i] = y / ytavg
                            edata.T[i] = e / ytavg

                    if 'normalize' in meta:
                        for i, (y, e) in enumerate(zip(ydata, edata)):
                            ydata[i] = y / meta['normalize']
                            edata[i] = e / meta['normalize']


                    if meta.get('timeselection') == "t":
                        xdata = tdata
                    elif meta.get('timeselection') == "lnt":
                        xdata = tdata
                        ax.set_xscale('log')
                    elif meta.get('timeselection') == "lnlnt":
                        with np.errstate(divide='ignore'):
                            xdata = np.log(np.log(tdata))
                    elif meta.get('timeselection') == "SE":
                        xdata = datanode['avg']['entanglement_entropy'][()]
                    elif meta.get('timeselection') == "lnSE":
                        xdata = np.log(datanode['avg']['entanglement_entropy'][()])
                    else:
                        print(f'Unexpected value: {meta.get("timeselection")=}')
                        xdata = tdata

                    if timenormalization := meta.get('timenormalization'):
                        t = get_timepoints(tdata, dbval)
                        idx_ent = t.idx_ent_saturated
                        sdata = datanode['avg']['entanglement_entropy'][()]
                        stavg = np.mean(sdata[idx_ent:], axis=0)
                        L = dbval['vals']['L']
                        if '/SEinfty' in timenormalization:
                            xdata = xdata / stavg
                        if '*SEinfty' in timenormalization:
                            xdata = xdata * stavg
                        if '/Spage' in timenormalization:
                            xdata = xdata / page_entropy(L)
                        if '*Spage' in timenormalization:
                            xdata = xdata * page_entropy(L)
                        if '/exp(L/4)' in timenormalization:
                            xdata = xdata / np.exp(L/4)
                        if '/exp(L/2)' in timenormalization:
                            xdata = xdata / np.exp((L)/2)
                        if '/exp(L)' in timenormalization:
                            xdata = xdata / np.exp((L))
                    for i, (y, e, colname) in enumerate(zip(ydata.T, edata.T, colnames)):
                        linestyle = meta['linestyle'][i] if 'linestyle' in meta and len(
                            meta['linestyle']) == len(ydata) else '-'

                        label = None
                        if 'entanglement' in colname or 'Neumann' in colname:
                            label = '$S_\mathrm{vN}$'
                        elif 'number' in colname:
                            label = '$S_\mathrm{N}$'
                        kappa = None # Slope of loglog fit of approach
                        alpha = None # Exponent to power fit of approach


                        if meta.get('plotsatapproach'):
                            sdata = datanode['avg']['entanglement_entropy'][()]
                            # idx1, idx2 = find_loglog_window2(tdata, sdata, dbval)
                            t = get_timepoints(tdata, dbval)
                            idx_sat = t.idx_num_saturated if 'number' in meta['colname'] else t.idx_ent_saturated
                            ytavg = np.mean(y[idx_sat:])
                            y = np.abs(y - ytavg)

                            if meta.get('marksatapproach'):
                                x1 = xdata[t.idx_num_lnlnt_begin]
                                x2 = xdata[t.idx_num_lnlnt_cease]
                                x3 = xdata[t.idx_num_saturated]
                                x4 = xdata[t.idx_ent_saturated]
                                y1 = y[t.idx_num_lnlnt_begin]
                                y2 = y[t.idx_num_lnlnt_cease]
                                y3 = y[t.idx_num_saturated]
                                y4 = y[t.idx_ent_saturated]
                                ax.plot([x1], [y1], color=color, marker='>',
                                        markersize=5.0,
                                        linestyle='None',
                                        markeredgecolor=color,
                                        # markeredgewidth=0.1,
                                        markerfacecolor='none',
                                        markeredgewidth=0.6,
                                        path_effects=mark_effects,
                                        zorder=10,
                                        )
                                ax.plot([x2], [y2], color=color, marker='s',
                                        markersize=5.0,
                                        linestyle='None',
                                        markeredgecolor=color,
                                        markerfacecolor='none',
                                        markeredgewidth=0.6,
                                        path_effects=mark_effects,
                                        zorder=10,
                                        )
                                # ax.plot([x3], [y3], color=color, marker='<',
                                #         markersize=4.0,
                                #         linestyle='None',
                                #         markeredgecolor='black',
                                #         markeredgewidth=0.1, path_effects=patch_effects, zorder=10,
                                #         )
                                # ax.plot([x4], [y4], color=color, marker='o',
                                #         markersize=4.0,
                                #         linestyle='None',
                                #         markeredgecolor='black',
                                #         markeredgewidth=0.1, path_effects=patch_effects, zorder=10,
                                #         )
                            if meta.get('fitsatapproachlnlnt'):
                                try:
                                    idx1 = t.idx_num_lnlnt_begin
                                    idx2 = t.idx_num_lnlnt_cease
                                    if idx2 <= idx1:
                                        raise IndexError("Invalid index order: idx1 {} | idx2 {}".format(idx1, idx2))
                                    if idx1 + 10 > idx2:
                                        raise IndexError("Too few datapoints for a fit: idx1 {} | idx2 {}".format(idx1, idx2))
                                    try:
                                        with np.errstate(divide='ignore', invalid='ignore'):
                                            tloglog = np.log(np.log(tdata))
                                        idx_min = int(idx1 * 0.5)
                                        idx_max = int(np.min([len(tdata) - 1, idx2 * 2.0]))
                                        popt, pcov = curve_fit(f=flinear, xdata=tloglog[idx1:idx2],
                                                                          ydata=y[idx1:idx2])
                                        tfit = xdata[idx_min:idx_max]
                                        yfit = flinear(tloglog[idx_min:idx_max], *popt)
                                        ax.plot(tfit, yfit, marker=None,
                                                linewidth=0.4,
                                                linestyle=(0, (4.1, 1.5)), label='fit', color=color,
                                                path_effects=yfit_effects,
                                                zorder=9,)

                                        kappa = popt[1]
                                    except IndexError as err:
                                        print(f'failed to plot approach: {e}')
                                        pass
                                    except ValueError as err:
                                        logger.error("Fit to loglog approach failed: {}".format(err))

                                except IndexError as e:
                                    pass
                            if meta.get('fitsatapproachpower'):
                                try:
                                    idx1 = t.idx_num_lnlnt_begin
                                    idx2 = t.idx_num_lnlnt_cease
                                    if idx2 <= idx1:
                                        raise IndexError("Invalid index order: idx1 {} | idx2 {}".format(idx1, idx2))
                                    if idx1 + 10 > idx2:
                                        raise IndexError("Too few datapoints for a fit: idx1 {} | idx2 {}".format(idx1, idx2))
                                    bounds = ([-np.inf, -np.inf], [np.inf, 0])
                                    with np.errstate(divide='ignore'):
                                        try:
                                            ylog = np.log10(y)
                                            tlog = np.log10(tdata)
                                            idx_min = int(idx1 * 0.75)
                                            idx_max = int(np.min([len(tdata) - 1, idx2 * 1.5]))
                                            # If y ~ At^alpha, then log(y) = log(A) + alpha * log(t), meaning we can do a linear
                                            # fit to the logx-logy data.
                                            popt_linear, pcov_linear = curve_fit(f=flinear, xdata=tlog[idx1:idx2],
                                                                   ydata=ylog[idx1:idx2], bounds=bounds)
                                            yfit_linear = 10 ** flinear(tlog, *popt_linear)
                                            ax.plot(xdata[idx_min:idx_max], yfit_linear[idx_min:idx_max], marker=None,
                                                    linewidth=0.4,
                                                    linestyle=(0, (1, 1.5)), label='fit', color=color,
                                                    path_effects=yfit_effects,
                                                    zorder=9,)
                                            alpha = popt_linear[1]
                                        except Exception as e:
                                            print(f'failed to plot approach: {e}')
                                            pass

                                except IndexError as e:
                                    pass


                        if meta.get('fillerror'):
                            line = ax.fill_between(x=xdata, y1=y - e, y2=y + e, alpha=0.10, label=None, color=color)


                        if meta.get('usescatterplot'):
                            line = ax.scatter(xdata, y, marker='o', s=7.0, edgecolors='gray', linewidths=0.2,
                                            # linestyle=linestyle,
                                            label=label, color=color,
                                            # linewidth=1.8,
                                            # path_effects=path_effects,
                                               zorder=8)
                        else:
                            line, = ax.plot(xdata, y, marker=None, linestyle=linestyle, label=label, color=color,
                                            # linewidth=1.2,
                                            path_effects=path_effects,zorder=8)


                        if i == 0:
                            legendrow = get_legend_row(db=db, datanode=datanode, legend_col_keys=legend_col_keys)
                            for icol, (col, key) in enumerate(zip(legendrow, legend_col_keys)):
                                key, fmt = key.split(':') if ':' in key else [key, '']

                                f['legends'][idx][icol]['handle'].append(line)
                                f['legends'][idx][icol]['label'].append(col)
                                f['legends'][idx][icol]['title'] = db['tex'][key]
                            if alpha is not None:
                                f['legends'][idx][icol+1]['handle'].append(line)
                                f['legends'][idx][icol+1]['label'].append(f'{alpha:.2f}')
                                f['legends'][idx][icol+1]['title'] = '$\\alpha$'
                                # f['legends'][idx][icol]['header'] = get_title(dbval, subspec, width=16)
                            if kappa is not None:
                                f['legends'][idx][icol+1]['handle'].append(line)
                                f['legends'][idx][icol+1]['label'].append(f'{100*abs(kappa):.2f}')
                                f['legends'][idx][icol+1]['title'] = '$100\kappa$'
                                # f['legends'][idx][icol]['header'] = get_title(dbval, subspec, width=16)
                        if meta.get('marksaturation') is True:
                            t = get_timepoints(tdata,dbval)
                            idx3, idx4 = t.idx_num_saturated, t.idx_ent_saturated
                            xm3.append(xdata[idx3])
                            ym3.append(y[idx3])
                            cm3.append(color)
                            xm4.append(xdata[idx4])
                            ym4.append(y[idx4])
                            cm4.append(color)

                        if meta.get('findloglogwindow') and 'entanglement_entropy' in datanode['avg'].dtype.fields:
                            sdata = datanode['avg']['entanglement_entropy'][()]
                            t = get_timepoints(tdata, dbval)
                            idx_sat = t.idx_num_saturated if 'number' in meta['colname'] else t.idx_ent_saturated
                            idx1, idx2 = t.idx_num_lnlnt_begin, t.idx_num_lnlnt_cease
                            print(f'L={dbval["vals"]["L"]:2} '
                                  f't▶={t.time_num_lnlnt_begin:.3e} '
                                  f't◼={t.time_num_lnlnt_cease:.3e} '
                                  f't◀={t.time_num_saturated:.3e} '
                                  f't●{t.time_ent_saturated:.3e} '
                                  f'num{len(tdata[idx_sat:])}')
                            f['ymax'] = np.max([f['ymax'], np.max(y)]) if f['ymax'] else np.max(y)
                            f['ymin'] = np.min([f['ymin'], y[idx1]]) if f['ymin'] else y[idx1]
                            if meta.get('markloglogwindow') is True:
                                # if linvals != linprod[-1]:
                                #     # This makes sure that we only do this once, at the last iteration
                                #     # Otherwise this gets too messy, full of dots
                                #     continue
                                # the x coords of this transformation are data, and the y coord are axes
                                yrange = np.abs(ax.get_ylim()[1]-ax.get_ylim()[0])
                                if len(f['ymarkoffset'])  <= idx:
                                    yreloffset = 0.02 * yrange
                                    f['ymarkoffset'].append(y[idx2] + yreloffset)
                                if meta.get('colname') == 'entanglement_entropy':
                                    idx3, idx4 = t.idx_num_saturated, t.idx_ent_saturated
                                    f['ymarkoffset'][idx] = y[idx4]
                                xm1.append(xdata[idx1])
                                ym1.append(y[idx1]),
                                cm1.append(color)
                                xm2.append(xdata[idx2])
                                ym2.append(y[idx2]),
                                cm2.append(color)
                                #
                                # ax.plot([xm1], [ym1], color=color, marker='>', markersize=4.5,
                                #         linestyle='None', markeredgecolor='w',
                                #                 markeredgewidth=0.0, path_effects=path_effects,zorder=10)
                                # ax.plot([xm2], [ym2], color=color, marker='s', markersize=4.5,
                                #         linestyle='None', markeredgecolor='w',
                                #         markeredgewidth=0.0, path_effects=path_effects,zorder=10)

                            if meta.get('fitloglogwindow') and not meta.get('zoomloglogwindow'):
                                # bounds = ([0., 0., 1.], [np.inf, np.inf, np.exp(1)])
                                try:
                                    # popt, pcov = curve_fit(f=floglog, xdata=tdata[idx1:idx2], ydata=y[idx1:idx2], bounds=bounds)
                                    if idx2 <= idx1:
                                        raise IndexError("Invalid index order: idx1 {} | idx2 {}".format(idx1, idx2))
                                    if idx1 + 10 > idx2:
                                        raise IndexError("Too few datapoints for a fit: idx1 {} | idx2 {}".format(idx1, idx2))

                                    bounds_v2 = ([-np.inf, 0], [np.inf, np.inf])
                                    with np.errstate(invalid='ignore', divide='ignore'):
                                        tloglog = np.log(np.log(tdata))
                                        popt, pcov = curve_fit(f=flinear, xdata=tloglog[idx1:idx2], ydata=y[idx1:idx2],
                                                               bounds=bounds_v2)
                                        idx_min = int(idx1 * 0.5)
                                        idx_max = int(np.min([len(tdata) - 1, idx2 * 1.5]))
                                        tfit = xdata[idx_min:idx_max]
                                        yfit = flinear(tloglog[idx_min:idx_max], *popt)
                                        ax.plot(tfit, yfit, marker=None,
                                                #linewidth=0.4,
                                                alpha=0.85,
                                                zorder=0,
                                                linestyle='--', label='fit', color=color,
                                                path_effects=path_effects)
                                except IndexError as err:
                                    pass
                                except ValueError as err:
                                    logger.error("Fit to loglog failed: {}".format(err))

                            if meta.get('zoomloglogwindow') and i in meta['zoomloglogwindow']['colnum']:
                                if ix is None:
                                    # pos tells where to put the inset, x0,y0, width, height in % units
                                    ix = ax.inset_axes(meta['zoomloglogwindow']['pos'], zorder=20)
                                ix.fill_between(x=xdata, y1=y - e, y2=y + e, alpha=0.15, label=None,
                                                color=color)
                                ix.plot(xdata, y, marker=None, linestyle=linestyle,
                                        label=None, color=color, path_effects=path_effects)
                                x1 = xdata[idx1]
                                x2 = xdata[idx2]
                                y1 = y[idx1]
                                y2 = y[idx2]

                                meta['zoomloglogwindow']['coords'] = [0.9*x1, 1.2*x2, y1-0.01, y2+0.01]
                                if meta.get('fitloglogwindow'):
                                    if idx1 + 5 < idx2:
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
                                                        linewidth=0.4,
                                                        linestyle='--', label=None, color=color, alpha=0.85,
                                                        path_effects=path_effects)
                                        except RuntimeError as e:
                                            pass
                                # if meta.get('markloglogwindow') is True:
                                    # yrange = np.abs(ix.get_ylim()[1] - ix.get_ylim()[0])
                                    # xm1 = x1
                                    # ym1 = f['ymarkoffset'][idx] + lidx/70 * yrange,
                                    # xm2 = x2
                                    # ym2 = f['ymarkoffset'][idx] + lidx/70 * yrange,
                                    # ix.plot([xm1, xm2],
                                    #         [ym1, ym2],
                                    #         color=color,
                                    #         marker=None, linestyle='dotted',
                                    #         path_effects=None, zorder=10)
                                    # ix.plot([xm1],[ym1],
                                    #         color=color,
                                    #         marker='>', markersize=4.5, linestyle='None', markeredgecolor='b',
                                    #         markeredgewidth=0.0, path_effects=None, zorder=10)
                                    # ix.plot([xm2],[ym2],
                                    #         color=color,
                                    #         marker='s', markersize=4.5, linestyle='None', markeredgecolor='b',
                                    #         markeredgewidth=0.0, path_effects=path_effects, zorder=10)


                    if not idx in f['axes_used']:
                        f['axes_used'].append(idx)
            if axtitle := get_default(meta, 'axtitle'):
                if dbidx == 0:
                    if dbval and isinstance(axtitle, bool):
                        axtitle = get_title(dbval, subspec, width=16)
                    ax.set_title(axtitle,
                                 horizontalalignment='center', x=0.5,fontstretch="ultra-condensed",
                                 )

            if xlabel := meta.get('xlabel'):
                ax.set_xlabel(xlabel)
            if ix is not None:
                ix.set_xscale(ax.get_xscale())

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
                # ix.patch.set_linewidth('1')
                # ix.patch.set_edgecolor('black')
                if 'legendtitle' in meta['zoomloglogwindow']:
                    ix.legend(title=meta['zoomloglogwindow']['legendtitle'], fontsize='x-small',
                              loc='lower center', framealpha=0.9,
                              bbox_to_anchor=(0.49, -0.01),
                              bbox_transform=ix.transAxes,
                              borderaxespad=0.5
                              )


                ax.indicate_inset_zoom(ix, edgecolor="black")

            for axs in [ax, ix]:
                if axs is None:
                    continue
                if meta.get('markloglogwindow') is not True:
                    continue
                yrange = np.abs(axs.get_ylim()[1] - axs.get_ylim()[0])
                ymmax = np.max(ym1 + ym2)
                ymmin = np.min(ym1 + ym2)
                for lidx,(x1, y1, c1, x2, y2, c2, x3, y3, c3, x4, y4, c4) in enumerate(zip(xm1, ym1, cm1, xm2,ym2,cm2, xm3,ym3,cm3, xm4,ym4,cm4)):
                    yl = np.max([y4, ]) if meta['colname'] == 'entanglement_entropy' else ymmin - (lidx) / 30 * yrange
                    if meta.get('markloglogwindow') is True:
                        axs.plot([x1, x2], [yl, yl], color=c1, marker=None, linestyle='dotted', path_effects=None, zorder=8) # Between rtriangle and square
                        # axs.plot([x1, x1], [y1, yl], color=c1, marker=None, linestyle='dotted', path_effects=None, zorder=8) # From rtriangle to line
                        # axs.plot([x2, x2], [y2, yl], color=c1, marker=None, linestyle='dotted', path_effects=None, zorder=8) # From square to line

                        axs.plot([x1], [yl], color=c1, marker='>',
                                # markersize=4.5,
                                linestyle='None',
                                markeredgecolor='lightgrey',
                                markeredgewidth=0.0, path_effects=patch_effects, zorder=10,
                                )
                        axs.plot([x2], [yl], color=c2, marker='s',
                                # markersize=4.5,
                                linestyle='None',
                                markeredgecolor='lightgrey',
                                markeredgewidth=0.0, path_effects=patch_effects, zorder=9,
                                )
                    if meta.get('marklogwindow') is True:
                        axs.plot([x2, x3], [yl, yl], color=c1, marker=None, linestyle='dotted', path_effects=None,
                                 zorder=8)  # Between square and rtriangle
                        # axs.plot([x3, x3], [y3, yl], color=c1, marker=None, linestyle='dotted', path_effects=None,
                        #          zorder=8)  # From rtriangle to line
                        axs.plot([x3], [yl], color=c3, marker='<',
                                # markersize=5,
                                linestyle='None',
                                markeredgecolor='lightgrey',
                                markeredgewidth=0.0, path_effects=patch_effects, zorder=10)
                    if meta.get('marksaturation') is True:
                        ax.plot([x4], [yl], color=c4, marker='o',
                                # markersize=5,
                                linestyle='None',
                                markeredgecolor='lightgrey',
                                markeredgewidth=0.0, path_effects=patch_effects, zorder=10)
                        ax.plot([x3, x4], [yl, yl], color=c4, marker=None, linestyle='dotted', path_effects=None,
                                 zorder=8)  # Between rtriangle and circle

            # if meta.get('marklogwindow') is True:
            #     for x3, y3, c3, x4, y4, c4 in zip(xm3,ym3,cm3, xm4,ym4,cm4):
            #         # yl = 0.5*(y3+y4)
            #         # yl = y4
            #         yl = np.max([y4, ])
            #         ax.plot([x3], [yl], color=c3, marker='<', markersize=5, linestyle='None', markeredgecolor='lightgrey',
            #                         markeredgewidth=0.35, path_effects=path_effects, zorder=10)

            # if meta.get('marksaturation') is True:
            #     for x3, y3, c3, x4, y4, c4 in zip(xm3,ym3,cm3, xm4,ym4,cm4):
            #         ax.plot([x4], [y4], color=c4, marker='o', markersize=5, linestyle='None', markeredgecolor='lightgrey',
            #                         markeredgewidth=0.35, path_effects=path_effects, zorder=10)
            #         ax.plot([x3, x4], [y4, y4], color=c4, marker=None, linestyle='dotted', path_effects=None,
            #                  zorder=8)  # Between rtriangle and circle

        if f['ymin']:
            f['ymin'] = 0.9 * f['ymin']
        if f['ymax']:
            f['ymax'] = 1.1 * f['ymax']

        if figspec_title := get_figspec_title(meta, dbval, figspec):
            f['fig'].suptitle(figspec_title)

        # prettify_plot4(fmeta=f, lgnd_meta=axes_legends)
        suffix = ''
        suffix = suffix + '_normpage' if 'normpage' in meta and meta['normpage'] else suffix
        suffix = suffix + '_loglog' if meta.get('timeselection') == 'lnlnt' else suffix
        f['filename'] = "{}/{}(t)_fig({})_sub({}){}".format(meta['plotdir'], meta['plotprefix'],
                                                       get_specvals(db, figspec, figvals),
                                                       get_specvals(db, subspec), suffix)
    return figs


def plot_v2_time_fig3_sub3_line1(db, meta, figspec, subspec, linspec, algo_filter=None, state_filter=None,
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
        logger.debug('- plotting figkeys: {}'.format(figkeys))
        dbval = None
        for idx, (subkeys, ax, ix) in enumerate(zip(subprod, f['ax'], f['ix'])):
            popt = None
            pcov = None
            logger.debug('-- plotting subkeys: {}'.format(subkeys))
            for dirkeys in dirprod:
                # palette = plt.rcParams['axes.prop_cycle'].by_key()['color']
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
                    ydata, colnames = get_table_data(datanode['avg'], meta['colname'],'f8')  # Supports multiple columns
                    edata, colnames = get_table_data(datanode['ste'], meta['colname'],'f8')  # Supports multiple columns
                    tdata = datanode['avg']['physical_time'][()]
                    ndata = datanode['avg']['num'][()]



                    if np.min(ndata) < 10:
                        continue
                    if len(tdata) <= 1:
                        logger.warning("tdata is too short: {}".format(tdata))
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

                    if meta.get('timeselection') == 'lnlnt':
                        with np.errstate(invalid='ignore'):
                            xdata = np.log(np.log(tdata))
                            if not 'xmin' in meta:
                                ax.set_xlim(xmin=-1)
                            if not 'xmax' in meta:
                                ax.set_xlim(xmax=1.05 * xdata[-1])
                    else:
                        xdata = tdata
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
                                        ax.plot(xdata, yfit, marker=None, linewidth=0.4,
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
                                # f['legends'][idx][icol]['header'] = get_title(dbval, subspec, width=16)

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
                                        ax.plot(tfit, yfit, marker=None, linewidth=0.4, alpha=1.0, zorder=50,
                                                linestyle='--', label='fit', color=color,
                                                path_effects=path_effects)
                                except IndexError as err:
                                    pass
                                except ValueError as err:
                                    logger.error(f"Fit failed: {err}")

                            if meta.get('zoomloglogwindow') and i in meta['zoomloglogwindow']['colnum']:
                                if ix is None:
                                    # pos tells where to put the inset, x0,y0, width, height in % units
                                    ix = ax.inset_axes(meta['zoomloglogwindow']['pos'], zorder=20)
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
                                    if idx1 + 5 < idx2:
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
                                                        linewidth=0.4,
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
                             horizontalalignment='center', x=0.45,
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


def plot_time_fig_sub_line(db, meta, figspec, subspec, linspec, algo_filter=None, state_filter=None, point_filter=None,
                           figs=None, palette_name=None, dbidx=0, dbnum=1):
    if db['version'] == 2:
        specs = figspec + subspec + linspec
        nonv3spec = lambda x: x not in ['cstd', 'tstd', 'tgw8', 'cgw8']
        specs = list(filter(nonv3spec, specs))
        fig3 = specs[:3]
        sub3 = specs[3:6]
        lin1 = specs[6:7]
        return plot_v2_time_fig3_sub3_line1(db=db, meta=meta, figspec=fig3, subspec=sub3, linspec=lin1,
                                            algo_filter=algo_filter, state_filter=state_filter,
                                            point_filter=point_filter, figs=figs, palette_name=palette_name)
    elif db['version'] == 3:
        return plot_v3_time_fig_sub_line(db=db, meta=meta, figspec=figspec, subspec=subspec, linspec=linspec,
                                         algo_filter=algo_filter, state_filter=state_filter, point_filter=point_filter,
                                         figs=figs, palette_name=palette_name,dbidx=dbidx,dbnum=dbnum)
    else:
        raise NotImplementedError('database version not implemented:' + db['version'])
