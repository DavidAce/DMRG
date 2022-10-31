from itertools import product
from src.database.database import *
from src.plotting.tools import *


def find_saturation_idx(ydata, std_threshold):
    sdata = []
    for i in range(len(ydata)):
        sdata.append(np.std(ydata[i:]))
    sdiff = -np.log(np.abs(np.diff(sdata)))
    return np.argmax(sdiff)


def find_saturation_idx2(ydata, threshold=1e-2):
    if len(ydata) <= 2:
        return
    ylog = -np.log10(ydata)
    ylog = ylog / ylog[-1]
    sdata = []
    w = 2
    for i, yl in enumerate(ylog):
        min_idx = np.min([len(ylog) - w, i])
        min_idx = np.max([min_idx, 0])
        s = np.std(ylog[min_idx:])
        sdata.append(s)
    idx = np.argwhere(np.asarray(sdata) < threshold)[0, 0]
    return idx


# def find_loglog_window(tdata, ydata, threshold1=0.6, threshold2=1e-2):
#     if len(ydata) <= 2:
#         return
#     ylog = -np.log10(ydata)
#     ylog = ylog / ylog[-1]
#     sdata = []
#     w = 2
#     for i, yl in enumerate(ylog):
#         min_idx = np.min([len(ylog) - w, i])
#         min_idx = np.max([min_idx, 0])
#         s = np.std(ylog[min_idx:])
#         sdata.append(s)
#     tdx = np.argwhere(np.asarray(tdata) >= 1.0)
#     idx1 = np.argwhere(np.asarray(sdata) < threshold1)
#     idx2 = np.argwhere(np.asarray(sdata) < threshold2)
#
#     tdx = tdx[0, 0] if len(tdx) > 0 else 0
#     idx1 = idx1[0, 0] if len(idx1) > 0 else 0
#     idx2 = idx2[0, 0] if len(idx2) > 0 else 0
#
#     idx1 = np.max([tdx, idx1])
#     return idx1, idx2


def page_entropy(L):
    n = int(2 ** (L / 2))
    S = - float(n - 1) / (2 * n)
    for k in range(1 + n, 1 + n ** 2):  # Include last
        S = S + 1.0 / k
    return S


def gather_time_series(mmntnode, yname=None, tname=None, aname=None):
    ydata = []
    tdata = []
    edata = []
    ndata = []
    adata = []
    for dsetkey, dsetpath, dsetnode in h5py_group_finder(node=mmntnode, keypattern=yname, num=1, dep=1):
        for iterkey, iterpath, iternode in h5py_group_iterator(node=dsetnode, dep=1):
            ydata.append(iternode['avg'][()])
            edata.append(iternode['ste'][()])
            ndata.append(iternode['num'][()])

    for dsetkey, dsetpath, dsetnode in h5py_group_finder(node=mmntnode, keypattern=tname, num=1, dep=1):
        for iterkey, iterpath, iternode in h5py_group_iterator(node=dsetnode, dep=1):
            tdata.append(iternode['avg'][()])
    for dsetkey, dsetpath, dsetnode in h5py_group_finder(node=mmntnode, keypattern=aname, num=1, dep=1):
        for iterkey, iterpath, iternode in h5py_group_iterator(node=dsetnode, dep=1):
            adata.append(iternode['avg'][()])
    if not ydata:
        raise LookupError("ydata is empty")
    if not tdata:
        raise LookupError("tdata is empty")
    if not adata:
        raise LookupError("tdata is empty")
    return np.array(ydata), np.array(tdata), np.array(edata), np.array(ndata), np.array(adata)


def plot_S_vs_Time_fig3_sub3_line1(db, meta, fig3, sub3, l1, algo_filter=None, state_filter=None, point_filter=None):
    if len(fig3) != 3:
        raise AssertionError("fig must have length 3")
    if len(sub3) != 3:
        raise AssertionError("sub must have length 3")
    if len(l1) != 1:
        raise AssertionError("itr must have length 1")

    for key0 in db['keys'][fig3[0]]:
        for key1 in db['keys'][fig3[1]]:
            for key2 in db['keys'][fig3[2]]:
                figrows, figcols = get_optimal_subplot_num(1 + len(db['keys'][sub3[0]]) * len(db['keys'][sub3[1]] * len(db['keys'][sub3[2]])))
                fig, axes = plt.subplots(nrows=figrows, ncols=figcols, figsize=(5 * figcols, 5 * figrows))
                fig.tight_layout(pad=5, w_pad=1.0, h_pad=1.0)
                fig.subplots_adjust(wspace=0.2, hspace=0.45)
                axes_used = []
                l1_legend = {'ax': [], 'handle': [], 'label': [], 'ncol': 1, 'unique': True, 'loc': 'upper center', 'insubfig': False}
                m1_legend = {'ax': [], 'handle': [], 'label': [], 'ncol': 1, 'unique': True, 'loc': 'lower center', 'insubfig': False}
                if l1[0] == 'b':
                    m1_legend['insubfig'] = True

                for idx, ((key3, key4, key5), ax) in enumerate(zip(product(db['keys'][sub3[0]], db['keys'][sub3[1]], db['keys'][sub3[2]]), np.ravel(axes))):
                    axes_used.append(idx)
                    for algokey, statekey, pointkey in product(db['keys']['algo'], db['keys']['state'], db['keys']['point']):
                        palette = sns.color_palette("Spectral", len(db['keys'][l1[0]]))
                        ymin = None
                        for key6, color in zip(db['keys'][l1[0]], palette):
                            findlist = [key0, key1, key2, key3, key4, key5, key6, algokey, statekey, pointkey, meta['dsetname']]
                            datanode = [value['datanode'] for key, value in db['dsets'].items() if all(k in key for k in findlist)]
                            if len(datanode) != 1:
                                print("found", len(datanode), "datanodes: ", datanode, " | findlist: ", findlist)
                                continue
                                raise LookupError("Found incorrect number of datanodes")

                            datanode = datanode[0]
                            mmntnode = datanode.parent
                            db_vals = db['dsets'][datanode.name]
                            ydata, tdata, edata, ndata, adata = gather_time_series(mmntnode, yname=meta['dsetname'],
                                                                                   tname="physical_time", aname="algorithm_time")

                            if np.min(ndata) < 10:
                                continue
                            print(np.min(ndata), datanode.name)

                            if meta['dsetname'] == 'number_entropy_midchain':
                                sdata, _, _, _, _ = gather_time_series(mmntnode, yname='entanglement_entropy_midchain',
                                                                       tname="physical_time", aname="algorithm_time")
                            else:
                                sdata = ydata

                            if meta['normpage']:
                                ydata = ydata / page_entropy(db_vals['L'])
                            if meta['timeloglevel'] == 2:
                                xdata = tdata
                                xdata = np.log(xdata + np.exp(1))
                                xdata = np.log(xdata)
                                # newidx = np.argwhere(tdata > 0)[:,0]
                                # xdata = tdata[newidx]
                                # tdata = tdata[newidx]
                                # ydata = ydata[newidx]
                                # sdata = sdata[newidx]
                                # edata = edata[newidx]
                                # ndata = ndata[newidx]

                            else:
                                xdata = tdata

                            # print(tdata)
                            ax.fill_between(x=xdata, y1=ydata - edata, y2=ydata + edata, alpha=0.15, label=None,
                                            color=color)
                            line, = ax.plot(xdata, ydata, marker=None, linewidth=1.2, linestyle='-', label=None, alpha=1.0, color=color,
                                            path_effects=[pe.SimpleLineShadow(offset=(0.6, -0.6), alpha=0.2), pe.Normal()])
                            l1_legend['ax'].append(ax)
                            l1_legend['handle'].append(line)
                            l1_legend['label'].append(
                                "${}={}$ (n:{}, ${}:{:>4.1f}$min)".format(l1[0], db_vals[l1[0]], np.min(ndata), "\\bar t_\mathrm{sim}", np.mean(adata) / 60))
                            idx1, idx2 = find_loglog_window2(tdata, sdata)
                            mark, = ax.plot([xdata[idx1], xdata[idx2]], [ydata[idx1], ydata[idx2]], color=color, marker='o', markersize=4, linestyle='None',
                                            path_effects=[pe.Stroke(linewidth=2, foreground='black'), pe.Normal()])
                            ymin = np.min([ymin, ydata[idx1]]) if ymin else ydata[idx1]
                            m1_legend['ax'].append(ax)
                            m1_legend['handle'].append(mark)
                            m1_legend['label'].append("$t_{}=$ {:.1e} $\\rightarrow$ {:.1e}".format('{\ln\ln}', tdata[idx1], tdata[idx2]))

                    ax.set_title("${}={}$, ${}={}$, ${}={}$ ".format(sub3[0], db_vals[sub3[0]], sub3[1], db_vals[sub3[1]], sub3[2], db_vals[sub3[2]]))
                    ax.set_xlabel("$t$")
                    ax.set_ylabel("$S(L/2)$")
                    ax.xaxis.set_tick_params(labelbottom=True)
                    ax.yaxis.set_tick_params(labelleft=True)
                    if meta['timeloglevel'] == 1:
                        ax.set_xscale('log')
                    if meta['timeloglevel'] == 2:
                        ax.set_xlabel("$\ln\ln(t+e)$")
                        if ymin:
                            ax.set_ylim(0.9 * ymin)

                fig.suptitle('{} entropy vs Time: ${}={}$,  ${}={}$, ${}={}$, {} realizations'.format(
                    meta['titlename'],
                    fig3[0], db_vals[fig3[0]],
                    fig3[1], db_vals[fig3[1]],
                    fig3[2], db_vals[fig3[2]], meta['realizations']))
                prettify_plot2(fig=fig, axes=axes, cols=figcols, rows=figrows, axes_used=axes_used, extra_legend=[l1_legend, m1_legend])
                suffix = ''
                suffix = suffix + '_normpage' if meta['normpage'] else suffix
                suffix = suffix + '_loglog' if meta['timeloglevel'] >= 2 else suffix
                plt.savefig(
                    "{}/{}(t)_fig({}_{}_{})_sub({}_{}_{}){}.pdf".format(meta['plotdir'], meta['plotprefix'], str(key0), str(key1), str(key2), sub3[0], sub3[1],
                                                                        sub3[2], suffix), format='pdf')


def plot_S_vs_Time_fig2_sub2_line2(db, meta, fig2, sub2, l2, algo_filter=None, state_filter=None, point_filter=None):
    if len(fig2) != 2:
        raise AssertionError("fig must have length 2")
    if len(sub2) != 2:
        raise AssertionError("sub must have length 2")
    if len(l2) != 2:
        raise AssertionError("itr must have length 2")

    for key0 in db['keys'][fig2[0]]:
        for key1 in db['keys'][fig2[1]]:
            figrows, figcols = get_optimal_subplot_num(1 + len(db['keys'][sub2[0]]) * len(db['keys'][sub2[1]]))
            fig, axes = plt.subplots(nrows=figrows, ncols=figcols, figsize=(3.5 * figcols, 3.5 * figrows))
            fig.tight_layout(pad=5, w_pad=1.0, h_pad=1.0)
            fig.subplots_adjust(wspace=0.25, hspace=0.4)
            axes_used = []
            l1_legend = {'handle': [], 'label': [], 'ncol': 1, 'unique': True, 'loc': 'upper left', 'color': None}
            l2_legend = {'handle': [], 'label': [], 'ncol': 1, 'unique': True, 'loc': 'upper right', 'color': 'dimgray'}

            for idx, ((key2, key3), ax) in enumerate(zip(product(db['keys'][sub2[0]], db['keys'][sub2[1]]), np.ravel(axes))):
                axes_used.append(idx)
                for algokey, statekey, pointkey in product(db['keys']['algo'], db['keys']['state'], db['keys']['point']):
                    palette = sns.color_palette("Spectral", len(db['keys'][l2[0]]), as_cmap=False, desat=1.0)
                    for key4, color in zip(db['keys'][l2[0]], palette):
                        lstyles = itertools.cycle(['-', ':', '--', '-.'])
                        for key5, lstyle in zip(db['keys'][l2[1]], lstyles):
                            findlist = [key0, key1, key2, key3, key4, key5, algokey, statekey, pointkey, meta['dsetname']]
                            datanode = [value['datanode'] for key, value in db['dsets'].items() if all(k in key for k in findlist)]
                            if len(datanode) != 1:
                                # print("found", len(datanode),"datanodes: ",datanode," | findlist: ", findlist)
                                # continue
                                raise LookupError("Found incorrect number of datanodes\nFound: {}\nFindlist: {}\n{}".format(len(datanode), findlist, datanode))

                            datanode = datanode[0]
                            mmntnode = datanode.parent
                            db_vals = db['dsets'][datanode.name]

                            ydata, xdata, edata, ndata = gather_time_series(mmntnode, yname=meta['dsetname'],
                                                                            tname="physical_time")
                            if meta['normpage']:
                                ydata = ydata / page_entropy(db_vals['L'])

                            if meta['timeloglevel'] >= 2:
                                xdata = np.log10(xdata)

                            ax.fill_between(x=xdata, y1=ydata - edata, y2=ydata + edata, alpha=0.2, label=None,
                                            color=color)
                            line, = ax.plot(xdata, ydata, marker=None, linewidth=1.8, linestyle=lstyle, label=None, alpha=1.0, color=color,
                                            path_effects=[pe.SimpleLineShadow(offset=(0.6, -0.6), alpha=0.2), pe.Normal()]
                                            )
                            l1_legend['handle'].append(line)
                            l1_legend['label'].append("${}={}$".format(l2[0], db_vals[l2[0]]))
                            l2_legend['handle'].append(
                                plt.Line2D([0], [0], color=l2_legend['color'], markersize=None, marker=None, linewidth=1.4, linestyle=lstyle))
                            l2_legend['label'].append("${}={}$".format(l2[1], db_vals[l2[1]]))

                ax.set_title("${}={}$, ${}={}$ ".format(sub2[0], db_vals[sub2[0]], sub2[1], db_vals[sub2[1]]))
                ax.set_xlabel("$t$")
                ax.set_ylabel("$S(L/2)$")
                ax.xaxis.set_tick_params(labelbottom=True)
                ax.yaxis.set_tick_params(labelleft=True)
                if meta['timeloglevel'] >= 1:
                    ax.set_xscale('log')

            fig.suptitle('{} entropy vs Time: ${}={}$, ${}={}$, {} realizations'.format(
                meta['titlename'],
                fig2[0], db_vals[fig2[0]],
                fig2[1], db_vals[fig2[1]], meta['realizations']))
            prettify_plot2(fig=fig, axes=axes, cols=figcols, rows=figrows, axes_used=axes_used, extra_legend=[l1_legend, l2_legend])
            suffix = ''
            suffix = suffix + '_normpage' if meta['normpage'] else suffix
            suffix = suffix + '_loglog' if meta['timeloglevel'] >= 2 else suffix
            plt.savefig("{}/{}(t)_fig({}_{})_sub({}_{}){}.pdf".format(meta['plotdir'], meta['plotprefix'], str(key0), str(key1), sub2[0], sub2[1], suffix),
                        format='pdf')


def plot_S_vs_Time_fig3_sub1_line2(db, meta, fig3, sub1, l2, algo_filter=None, state_filter=None, point_filter=None):
    if len(fig3) != 3:
        raise AssertionError("fig must have length 3")
    if len(sub1) != 1:
        raise AssertionError("sub must have length 1")
    if len(l2) != 2:
        raise AssertionError("itr must have length 2")

    for key0 in db['keys'][fig3[0]]:
        for key1 in db['keys'][fig3[1]]:
            for key2 in db['keys'][fig3[2]]:
                figrows, figcols = get_optimal_subplot_num(1 + len(db['keys'][sub1[0]]))
                fig, axes = plt.subplots(nrows=figrows, ncols=figcols, figsize=(3.5 * figcols, 3.5 * figrows))
                fig.tight_layout(pad=5, w_pad=1.0, h_pad=1.0)
                fig.subplots_adjust(wspace=0.25, hspace=0.4)
                axes_used = []
                l1_legend = {'handle': [], 'label': [], 'ncol': 1, 'unique': True, 'loc': 'upper left', 'color': None}
                l2_legend = {'handle': [], 'label': [], 'ncol': 1, 'unique': True, 'loc': 'upper right', 'color': 'dimgray'}

                for idx, (key3, ax) in enumerate(zip(db['keys'][sub1[0]], np.ravel(axes))):
                    axes_used.append(idx)
                    for algokey, statekey, pointkey in product(db['keys']['algo'], db['keys']['state'], db['keys']['point']):
                        palette = sns.color_palette("Spectral", len(db['keys'][l2[0]]), as_cmap=False, desat=1.0)
                        palette.reverse()
                        for key4, color in zip(db['keys'][l2[0]], palette):
                            lstyles = itertools.cycle(['-', ':', '--', '-.'])
                            for key5, lstyle in zip(db['keys'][l2[1]], lstyles):
                                findlist = [key0, key1, key2, key3, key4, key5, algokey, statekey, pointkey, meta['dsetname']]
                                datanode = [value['datanode'] for key, value in db['dsets'].items() if all(k in key for k in findlist)]
                                if len(datanode) != 1:
                                    print("found", len(datanode), "datanodes: ", datanode, " | findlist: ", findlist)
                                    continue
                                    raise LookupError("Found incorrect number of datanodes")

                                datanode = datanode[0]
                                mmntnode = datanode.parent
                                db_vals = db['dsets'][datanode.name]

                                ydata, xdata, edata, ndata = gather_time_series(mmntnode, yname=meta['dsetname'],
                                                                                tname="physical_time")
                                if meta['normpage']:
                                    ydata = ydata / page_entropy(db_vals['L'])

                                ax.fill_between(x=xdata, y1=ydata - edata, y2=ydata + edata, alpha=0.2, label=None,
                                                color=color)
                                line, = ax.plot(xdata, ydata, marker=None, linewidth=1.8, linestyle=lstyle, label=None, alpha=1.0, color=color,
                                                path_effects=[pe.SimpleLineShadow(offset=(0.6, -0.6), alpha=0.2), pe.Normal()]
                                                )
                                l1_legend['handle'].append(line)
                                l1_legend['label'].append("${}={}$".format(l2[0], db_vals[l2[0]]))
                                l2_legend['handle'].append(
                                    plt.Line2D([0], [0], color=l2_legend['color'], markersize=None, marker=None, linewidth=1.4, linestyle=lstyle))
                                l2_legend['label'].append("${}={}$".format(l2[1], db_vals[l2[1]]))

                    ax.set_title("${}={}$".format(sub1[0], db_vals[sub1[0]]))
                    ax.set_xlabel("$t$")
                    ax.set_ylabel("$S(L/2)$")
                    ax.xaxis.set_tick_params(labelbottom=True)
                    ax.yaxis.set_tick_params(labelleft=True)
                    ax.set_xscale('log')

                fig.suptitle('{} entropy vs Time: ${}={}$, ${}={}$, ${}={}$, {} realizations'.format(
                    meta['titlename'],
                    fig3[0], db_vals[fig3[0]],
                    fig3[1], db_vals[fig3[1]],
                    fig3[2], db_vals[fig3[2]],
                    meta['realizations']))
                prettify_plot2(fig=fig, axes=axes, cols=figcols, rows=figrows, axes_used=axes_used, extra_legend=[l1_legend, l2_legend])
                suffix = '_normpage' if meta['normpage'] else ''
                plt.savefig(
                    "{}/{}(t)_fig({}_{}_{})_sub({}){}.pdf".format(meta['plotdir'], meta['plotprefix'], str(key0), str(key1), str(key2), sub1[0], suffix),
                    format='pdf')
