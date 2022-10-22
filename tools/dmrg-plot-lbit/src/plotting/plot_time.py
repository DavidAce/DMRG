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


def plot_v2_time_fig3_sub3_line1(db, meta, fig3, sub3, l1, algo_filter=None, state_filter=None, point_filter=None, f=None, palette_name=None):
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

    for key0 in get_keys(db, fig3[0]):
        for key1 in get_keys(db, fig3[1]):
            for key2 in get_keys(db, fig3[2]):
                keyprod = list(product(*get_keys(db, sub3)))
                numplots = len(keyprod)
                if not f:
                    f = get_fig_meta(numplots, meta=meta)
                for idx, ((key3, key4, key5), ax) in enumerate(zip(keyprod, f['ax'])):
                    popt = None
                    pcov = None
                    dbval = None
                    for algokey, statekey, cronokey in product(db['keys']['algo'], db['keys']['state'], db['keys']['crono']):
                        # palette = plt.rcParams['axes.prop_cycle'].by_key()['color']
                        palette = sns.color_palette(palette=palette_name, n_colors=len(db['keys'][l1[0]]))
                        for key6, color in zip(get_keys(db, l1[0]), palette):
                            findlist = [key0, key1, key2, key3, key4, key5, key6, algokey, statekey, cronokey,
                                        meta['groupname']]
                            datanode = [value['node']['data'] for key, value in db['dsets'].items() if
                                        all(k in key for k in findlist)]
                            if len(datanode) != 1:
                                print("ERROR: found", len(datanode), "datanodes: ", datanode, " | findlist: ", findlist)
                                continue
                                # raise LookupError("Found incorrect number of datanodes")
                            datanode = datanode[0]
                            dbval = db['dsets'][datanode.name]
                            ydata, colnames = get_table_data(datanode['avg'], meta['colname'], 'f8')  # Supports multiple columns
                            edata, colnames = get_table_data(datanode['ste'], meta['colname'], 'f8')  # Supports multiple columns
                            tdata = datanode['avg']['physical_time'][()]
                            ndata = datanode['avg']['num'][()]

                            if np.min(ndata) < 10:
                                continue
                            if len(tdata) <= 1:
                                print("WARNING: tdata is too short: ", tdata)
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

                            if meta.get('timeloglevel') == 2:
                                with np.errstate(invalid='ignore'):
                                    xdata = np.log(np.log(tdata))
                                    ax.set_xlim(-1, 1.05 * xdata[-1])

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
                                    sdata = datanode['avg']['entanglement_entropy_midchain'][()]
                                    idx1, idx2 = find_loglog_window2(tdata, sdata, dbval)
                                    ysat = np.mean(ydata[idx2:])  # Saturation value
                                    y = np.abs(y - ysat)
                                    try:
                                        if idx2 <= idx1:
                                            raise RuntimeError("Invalid values: idx1 {} | idx2 {}", idx1, idx2)
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
                                                ax.annotate('$\sim t^{{{:.2f}}}$'.format(popt[1]), xy=(xmid, ymid), xytext=(xtxt, ytxt),
                                                            arrowprops=dict(arrowstyle="->", color=color))
                                            except:
                                                pass


                                    except RuntimeError as e:
                                        raise RuntimeError(e)


                                else:
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

                                if meta.get('findloglogwindow') and 'entanglement_entropy_midchain' in datanode['avg'].dtype.fields:
                                    sdata = datanode['avg']['entanglement_entropy_midchain'][()]
                                    idx1, idx2 = find_loglog_window2(tdata, sdata, dbval)
                                    f['ymax'] = np.max([f['ymax'], np.max(y)]) if f['ymax'] else np.max(y)
                                    f['ymin'] = np.min([f['ymin'], y[idx1]]) if f['ymin'] else y[idx1]
                                    if meta.get('markloglogwindow'):
                                        mark, = ax.plot([xdata[idx1], xdata[idx2]], [y[idx1], y[idx2]],
                                                        color=color,
                                                        marker='o', markersize=6, linestyle='None', markeredgecolor='w',
                                                        path_effects=path_effects)

                                    if meta.get('fitloglogwindow') and not meta.get('zoomloglogwindow'):
                                        if idx + 5 < idx2:
                                            # bounds = ([0., 0., 1.], [np.inf, np.inf, np.exp(1)])
                                            try:
                                                # popt, pcov = curve_fit(f=floglog, xdata=tdata[idx1:idx2], ydata=y[idx1:idx2], bounds=bounds)
                                                if idx2 <= idx1:
                                                    raise RuntimeError("Invalid values: idx1 {} | idx2 {}", idx1, idx2)
                                                bounds_v2 = ([-np.inf, 0], [np.inf, np.inf])
                                                with np.errstate(invalid='ignore'):
                                                    tloglog = np.log(np.log(tdata))
                                                    popt, pcov = curve_fit(f=flinear, xdata=tloglog[idx1:idx2], ydata=y[idx1:idx2], bounds=bounds_v2)
                                                    print("{} | idx {} {} | {}".format(popt, idx1, idx2, findlist))
                                                    ax.plot(xdata, flinear(tloglog, *popt), marker=None, linewidth=0.8,
                                                            linestyle='--', label='fit', color=color,
                                                            path_effects=path_effects)
                                            except RuntimeError as e:
                                                pass

                                    if meta.get('zoomloglogwindow') and i in meta['zoomloglogwindow']['colnum']:
                                        if f['axin'] is None:
                                            # pos tells where to put the inset, x0,y0, width, height in % units
                                            f['axin'] = ax.inset_axes(meta['zoomloglogwindow']['pos'])
                                        f['axin'].fill_between(x=xdata, y1=y - e, y2=y + e, alpha=0.15, label=None,
                                                               color=color)
                                        f['axin'].plot(xdata, y, marker=None, linestyle=linestyle,
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
                                                        f['axin'].plot(xdata, flinear(tloglog, *popt), marker=None,
                                                                       linewidth=0.8,
                                                                       linestyle='--', label=None, color=color,
                                                                       path_effects=path_effects)
                                                except RuntimeError as e:
                                                    pass
                                        if meta.get('markloglogwindow'):
                                            f['axin'].plot([xdata[idx1], xdata[idx2]], [y[idx1], y[idx2]],
                                                           color=color,
                                                           marker='o', markersize=6, linestyle='None', markeredgecolor='w',
                                                           path_effects=path_effects)

                            if not idx in f['axes_used']:
                                f['axes_used'].append(idx)
                    if dbval:
                        ax.set_title(get_title(dbval, sub3, width=16), fontstretch="ultra-condensed", bbox=dict(facecolor='white', alpha=1.0))

                    ax.set_xlabel("$t$")
                    if meta.get('timeloglevel'):
                        if meta['timeloglevel'] == 1:
                            ax.set_xscale('log')
                            ymin = None
                            ymax = None
                            if f['axin'] is not None:
                                f['axin'].set_xscale('log')

                        if meta['timeloglevel'] == 2:
                            ax.set_xlabel("$\ln\ln t$")

                    if meta.get('zoomloglogwindow') and f['axin'] is not None:
                        x1, x2, y1, y2 = meta['zoomloglogwindow']['coords']  # sub region of the original image
                        f['axin'].set_xlim(xmin=0.1 * x1, xmax=10 * x2)
                        f['axin'].set_ylim(ymin=0.975 * y1, ymax=1.025 * y2)
                        # f['axin'].set_xticklabels('')
                        # f['axin'].set_yticklabels('')
                        f['axin'].tick_params(axis='both', which='both', labelsize='x-small')
                        # f['axin'].xaxis.set_major_locator(plt.MaxNLocator(5))
                        f['axin'].xaxis.set_major_locator(plt.LogLocator(base=10, numticks=6))
                        if 'timeloglevel' in meta and meta['timeloglevel'] == 2:
                            f['axin'].yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
                            f['axin'].xaxis.set_major_formatter(FormatStrFormatter('%.1f'))

                        f['axin'].yaxis.tick_right()
                        f['axin'].patch.set_linewidth('2')
                        f['axin'].patch.set_edgecolor('black')
                        f['axin'].legend(title=meta['zoomloglogwindow']['legendtitle'], fontsize='x-small',
                                         loc='lower center', framealpha=0.9,
                                         bbox_to_anchor=(0.49, -0.01),
                                         bbox_transform=f['axin'].transAxes,
                                         borderaxespad=0.5
                                         )
                        ax.indicate_inset_zoom(f['axin'], edgecolor="black")


                if f['ymin']:
                    f['ymin'] = 0.9 * f['ymin']
                if f['ymax']:
                    f['ymax'] = 1.1 * f['ymax']

                if not prb_style and dbval:
                    f['fig'].suptitle('{} vs Time\n{}'.format(meta['titlename'], get_title(dbval, fig3)))

                # prettify_plot4(fmeta=f, lgnd_meta=axes_legends)
                if not f['filename']:
                    suffix = ''
                    suffix = suffix + '_normpage' if 'normpage' in meta and meta['normpage'] else suffix
                    suffix = suffix + '_loglog' if 'timeloglevel' in meta and meta['timeloglevel'] >= 2 else suffix
                    f['filename'] = "{}/{}(t)_fig({}_{}_{})_sub({}_{}_{}){}".format(meta['plotdir'], meta['plotprefix'],
                                                                                    str(key0), str(key1), str(key2), sub3[0], sub3[1], sub3[2], suffix)

    return f


def plot_v2_time_fig3_sub2_line2(db, meta, fig3, sub2, l2, algo_filter=None, state_filter=None, point_filter=None):
    if len(fig3) != 3:
        raise AssertionError("fig must have length 3")
    if len(sub2) != 2:
        raise AssertionError("sub must have length 2")
    if len(l2) != 2:
        raise AssertionError("itr must have length 2")

    if 'mplstyle' in meta:
        plt.style.use(meta['mplstyle'])

    l2_legend = {'handle': [], 'label': [], 'ncol': 1, 'unique': True, 'loc': 'upper left', 'insubfig': False,
                 'title': []}
    m1_legend = {'handle': [], 'label': [], 'ncol': 1, 'unique': True, 'loc': 'upper left', 'insubfig': False,
                 'title': ["$t_{}$".format('{\ln\ln}')]}

    for key in l2:
        l2_legend['title'].append(db['tex'][key])

    l2_legend['title'].extend(['$n$', '$\\bar \chi$', '$\\bar t_\mathrm{sim}$'])

    for key0 in db['keys'][fig3[0]]:
        for key1 in db['keys'][fig3[1]]:
            for key2 in db['keys'][fig3[2]]:
                figrows, figcols = get_optimal_subplot_num(len(db['keys'][sub2[0]]) * len(db['keys'][sub2[1]]))
                fig, axes = plt.subplots(nrows=figrows, ncols=figcols, figsize=(1.25 * 5 * figcols, 5 * figrows),
                                         sharey='all')
                axes_used = []

                axes_legends = []
                ymax = None
                ymin = None

                # if l1[0] == 'b':
                #     m1_legend['insubfig'] = True

                for idx, ((key3, key4), ax) in enumerate(
                        zip(product(db['keys'][sub2[0]], db['keys'][sub2[1]]), np.ravel(axes))):
                    l = deepcopy(l2_legend)
                    m = deepcopy(m1_legend)
                    axin = None
                    if len(l['label']) != 0:
                        raise AssertionError('l is not fresh')
                    for algokey, statekey, cronokey in product(db['keys']['algo'], db['keys']['state'],
                                                               db['keys']['crono']):
                        palette_names = get_uniform_palette_names(len(db['keys'][l2[0]]))
                        for key5, palette_name in zip(db['keys'][l2[0]], palette_names):
                            # palette = plt.rcParams['axes.prop_cycle'].by_key()['color']
                            palette = sns.color_palette(palette_name, len(db['keys'][l2[1]]))
                            for key6, color in zip(db['keys'][l2[1]], palette):
                                findlist = [key0, key1, key2, key3, key4, key5, key6, algokey, statekey, cronokey,
                                            meta['groupname']]
                                datanode = [value['datanode'] for key, value in db['dsets'].items() if
                                            all(k in key for k in findlist)]
                                if len(datanode) != 1:
                                    print("found", len(datanode), "datanodes: ", datanode, " | findlist: ", findlist)
                                    continue
                                    # raise LookupError("Found incorrect number of datanodes")

                                datanode = datanode[0]
                                dbval = db['dsets'][datanode.name]
                                if isinstance(meta['colname'], list):  # Support for plotting multiple quantities
                                    ydata = [datanode['avg'][col][()] for col in meta['colname']]
                                    edata = [datanode['ste'][col][()] for col in meta['colname']]
                                else:
                                    ydata = [datanode['avg'][meta['colname']][()]]
                                    edata = [datanode['ste'][meta['colname']][()]]

                                tdata = datanode['avg']['physical_time'][()]
                                adata = datanode['avg']['algorithm_time'][()]
                                ndata = datanode['avg']['num'][()]
                                bdata = datanode['avg']['bond_dimension_midchain'][()]

                                if np.min(ndata) < 10:
                                    continue
                                print(np.min(ndata), datanode.name)

                                if meta['normpage']:
                                    for i, (y, e) in enumerate(zip(ydata, edata)):
                                        p = page_entropy(dbval['L'])
                                        ydata[i] = y / p
                                        edata[i] = e / p

                                if meta['timeloglevel'] == 2:
                                    xdata = tdata
                                    xdata = np.log(xdata + np.exp(1))
                                    xdata = np.log(xdata)

                                else:
                                    xdata = tdata

                                for i, (y, e) in enumerate(zip(ydata, edata)):
                                    linestyle = meta['linestyle'][i] if 'linestyle' in meta and len(
                                        meta['linestyle']) == len(ydata) else '-'

                                    ax.fill_between(x=xdata, y1=y - e, y2=y + e, alpha=0.15, label=None, color=color)
                                    line, = ax.plot(xdata, y, marker=None, linewidth=1.2, linestyle=linestyle,
                                                    label=None, alpha=1.0, color=color,
                                                    path_effects=[pe.SimpleLineShadow(offset=(0.6, -0.6), alpha=0.2),
                                                                  pe.Normal()])

                                    if i == 0:
                                        l['handle'].append(line)
                                        l['label'].append([])
                                        for key in l2:
                                            l['label'][-1].extend(['{}'.format(dbval[key])])
                                        l['label'][-1].extend(
                                            [
                                                '{}'.format(np.min(ndata)),
                                                '{:>.1f}'.format(np.max(bdata)),
                                                '{:>.1f}m'.format(np.max(adata) / 60)
                                            ])

                                    if 'findloglogwindow' in meta and meta.get(
                                            'findloglogwindow') and 'entanglement_entropy_midchain' in datanode[
                                        'avg'].dtype.fields:
                                        sdata = datanode['avg']['entanglement_entropy_midchain'][()]
                                        idx1, idx2 = find_loglog_window2(tdata, sdata)
                                        ymax = np.max([ymax, np.max(y)]) if ymax else np.max(y)
                                        ymin = np.min([ymin, y[idx1]]) if ymin else y[idx1]
                                        if meta.get('markloglogwindow'):
                                            mark, = ax.plot([xdata[idx1], xdata[idx2]], [y[idx1], y[idx2]], color=color,
                                                            marker='o', markersize=4, linestyle='None',
                                                            path_effects=[pe.Stroke(linewidth=2, foreground='black'),
                                                                          pe.Normal()])
                                            if i == 0:
                                                m['handle'].append(mark)
                                                m['label'].append(["{:.1e}".format(tdata[idx2] - tdata[idx1])])

                                        if meta.get('zoomloglogwindow') and i in meta['zoomloglogwindow']['colnum']:
                                            if not axin:
                                                # pos tells where to put the inset, x0,y0, width, height in % units
                                                axin = ax.inset_axes(meta['zoomloglogwindow']['pos'])
                                            ax.fill_between(x=xdata, y1=y - e, y2=y + e, alpha=0.15, label=None,
                                                            color=color)
                                            axin.plot(xdata, y, marker=None, linewidth=1.2, linestyle=linestyle,
                                                      label=None, alpha=1.0, color=color,
                                                      path_effects=[pe.SimpleLineShadow(offset=(0.6, -0.6), alpha=0.2),
                                                                    pe.Normal()])
                                            x1, x2, y1, y2 = meta['zoomloglogwindow']['coords']
                                            x1 = np.min([x1, xdata[idx1]]) if x1 else xdata[idx1]
                                            x2 = np.max([x2, xdata[-1]]) if x2 else xdata[-1]
                                            y1 = np.min([y1, y[idx1]]) if y1 else y[idx1]
                                            y2 = np.max([y2, y[idx2]]) if y2 else y[idx2]
                                            meta['zoomloglogwindow']['coords'] = [x1, x2, y1, y2]

                    axes_legends.append(dict({'ax': ax, 'legends': ['l', 'm'], 'l': l.copy(), 'm': m.copy()}))
                    axes_used.append(idx)

                    if dbval:
                        ax.set_title(get_title(dbval, sub2), x=0.5, horizontalalignment='left')
                    # ax.set_title("{}, {}".format(dbval['texeq'][sub2[0]], dbval['texeq'][sub2[1]]))
                    ax.set_xlabel("$t$")

                    ax.xaxis.set_tick_params(labelbottom=True)
                    ax.yaxis.set_tick_params(labelleft=True)
                    if meta['timeloglevel'] == 1:
                        ax.set_xscale('log')
                        ymin = None
                        ymax = None
                    if meta['timeloglevel'] == 2:
                        ax.set_xlabel("$\ln\ln(t+e)$")
                    if 'zoomloglogwindow' in meta and axin:
                        x1, x2, y1, y2 = meta['zoomloglogwindow']['coords']  # sub region of the original image
                        axin.set_xlim(x1, x2)
                        axin.set_ylim(y1 * 0.95, y2 * 1.05)
                        # axin.set_xticklabels('')
                        # axin.set_yticklabels('')
                        axin.tick_params(axis='both', which='both', labelsize='x-small')
                        axin.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
                        axin.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
                        axin.yaxis.tick_right()
                        axin.patch.set_linewidth('2')
                        axin.patch.set_edgecolor('black')
                        axin.legend(title=meta['zoomloglogwindow']['legendtitle'], fontsize='x-small',
                                    loc='lower center', framealpha=1.0,
                                    bbox_to_anchor=(0.49, -0.01),
                                    bbox_transform=axin.transAxes,
                                    borderaxespad=0.
                                    )

                        ax.indicate_inset_zoom(axin, edgecolor="black")

                if ymin:
                    ymin = 0.9 * ymin
                if ymax:
                    ymax = 1.1 * ymax
                fig.suptitle('{} vs Time\n{}, {}, {}'.format(
                    meta['titlename'],
                    dbval['texeq'][fig3[0]],
                    dbval['texeq'][fig3[1]],
                    dbval['texeq'][fig3[2]]))
                plt.subplots_adjust(top=0.9)
                prettify_plot3(fig=fig, axes=axes, cols=figcols, rows=figrows, axes_used=axes_used, ymin=ymin,
                               ymax=ymax, extra_legend=axes_legends)
                # fig.tight_layout(w_pad=10.5, rect=[0, 0, 0.85, 1])
                # fig.subplots_adjust(wspace=0.5, hspace=0)
                # plt.subplots_adjust(right=0.5)
                suffix = ''
                suffix = suffix + '_normpage' if 'normpage' in meta and meta['normpage'] else suffix
                suffix = suffix + '_loglog' if 'timeloglevel' in meta and meta['timeloglevel'] >= 2 else suffix

                plt.savefig(
                    "{}/{}(t)_fig({}_{}_{})_sub({}_{}){}.pdf".format(meta['plotdir'], meta['plotprefix'], str(key0),
                                                                     str(key1), str(key2), sub2[0], sub2[1], suffix),
                    bbox_inches="tight", format='pdf')
                plt.savefig(
                    "{}/{}(t)_fig({}_{}_{})_sub({}_{}){}.png".format(meta['plotdir'], meta['plotprefix'], str(key0),
                                                                     str(key1), str(key2), sub2[0], sub2[1], suffix),
                    bbox_inches="tight", format='png')
