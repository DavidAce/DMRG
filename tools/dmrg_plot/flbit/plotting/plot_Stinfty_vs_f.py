from itertools import product
import numpy as np
from database.database import *
from plotting.tools import *
import seaborn as sns
import matplotlib.patheffects as pe


def find_saturation_idx(ydata, std_threshold):
    # Consider Y vs X: a noisy signal decaying in the shape of a hockey-club, say.
    # We want to identify the point at which the signal stabilizes. We use the fact that the
    # standard deviation is high if it includes parts of the non-stable signal, and low if
    # it includes only the stable part.
    # Here we monitor the standard deviation of the signal between [start_point, end_point],
    # and move "start_point" towards the end. If the standard deviation goes below a certain
    # threshold, i.e. threshold < max_std, then we have found the stabilization point.
    sdata = []
    for i in range(len(ydata)):
        sdata.append(np.std(ydata[i:]))
    sdiff = -np.log(np.abs(np.diff(sdata)))
    return np.argmax(sdiff)


def gather_at_tinfty(mmntnode, yname=None):
    ydata = []
    # tdata = []
    # edata = []
    # ndata = []
    for dsetkey, dsetpath, dsetnode in h5py_group_finder(node=mmntnode, keypattern=yname, num=1, dep=1):
        for iterkey, iterpath, iternode in h5py_group_iterator(node=dsetnode, dep=1):
            ydata.append(iternode['avg'][()])
            # edata.append(iternode['ste'][()])
            # ndata.append(iternode['num'][()])

    # for dsetkey,dsetpath,dsetnode in h5py_group_finder(node=mmntnode, keypattern=tname, num=1, dep=1):
    #     for iterkey,iterpath,iternode in h5py_group_iterator(node=dsetnode, dep=1):
    #         tdata.append(iternode['avg'][()])
    if not ydata:
        raise LookupError("ydata is empty")
    # if not tdata:
    #     raise LookupError("tdata is empty")

    # Instead of returning the whole time series, we want the average saturation value.
    # So first we need to find the index where y saturates
    idx = find_saturation_idx(ydata, np.max(ydata) / 10)
    ydata = ydata[idx:]
    return np.mean(ydata), np.std(ydata) / np.sqrt(len(ydata)), idx


def page_entropy(L):
    n = int(2 ** (L / 2))
    S = - float(n - 1) / (2 * n)
    for k in range(1 + n, 1 + n ** 2):  # Include last
        S = S + 1.0 / k
    return S


# For plots S(t -> infty) vs f, each line has unique u

def plot_Stinfty_fig2_sub2_line1(db, meta, fig2, sub2, l1, x1, algo_filter=None, state_filter=None, point_filter=None):
    if len(fig2) != 2:
        raise AssertionError("fig must have length 2")
    if len(sub2) != 2:
        raise AssertionError("sub must have length 2")
    if len(l1) != 1:
        raise AssertionError("itr must have length 1")

    for key0 in db['keys'][fig2[0]]:
        for key1 in db['keys'][fig2[1]]:
            figrows, figcols = get_optimal_subplot_num(1 + len(db['keys'][sub2[0]]) * len(db['keys'][sub2[1]]))
            fig, axes = plt.subplots(nrows=figrows, ncols=figcols, figsize=(3.5 * figcols, 3.5 * figrows), sharey='all', sharex='all')
            fig.tight_layout(pad=5, w_pad=1.0, h_pad=1.0)
            fig.subplots_adjust(wspace=0.25, hspace=0.40)
            axes_used = []
            l1_legend = {'handle': [], 'label': [], 'ncol': 1, 'unique': True, 'loc': 'upper center'}

            for idx, ((key2, key3), ax) in enumerate(zip(product(db['keys'][sub2[0]], db['keys'][sub2[1]]), np.ravel(axes))):
                print(idx, key2, key3)
                axes_used.append(idx)
                for algokey, statekey, pointkey in product(db['keys']['algo'], db['keys']['state'], db['keys']['point']):
                    palette = sns.color_palette("husl", len(db['keys'][l1[0]]))
                    for key4, color in zip(db['keys'][l1[0]], palette):
                        ydata = []
                        xdata = []
                        edata = []
                        for key5 in db['keys'][x1[0]]:
                            findlist = [key0, key1, key2, key3, key4, key5, algokey, statekey, pointkey, meta['dsetname']]
                            datanode = [value['datanode'] for key, value in db['dsets'].items() if all(k in key for k in findlist)]
                            if len(datanode) != 1:
                                print("found", len(datanode), "datanodes: ", datanode, " | findlist: ", findlist)
                                continue
                                raise LookupError("Found incorrect number of datanodes")

                            datanode = datanode[0]
                            mmntnode = datanode.parent
                            db_vals = db['dsets'][datanode.name]
                            x = db_vals[x1[0]]
                            cachestr = '{}_tinfty'.format(meta['plotprefix'])
                            if cachestr in db_vals:
                                y, e, i = db_vals[cachestr]
                            else:
                                y, e, i = gather_at_tinfty(mmntnode, yname=meta['dsetname'])
                                db_vals[cachestr] = y, e, i
                            ydata.append(y)
                            xdata.append(x)
                            edata.append(e)

                        ydata = np.asarray(ydata)
                        edata = np.asarray(edata)
                        xdata = np.asarray(xdata)
                        if meta['normpage']:
                            ydata = ydata / page_entropy(db_vals['L'])

                        line, = ax.plot(xdata, ydata, marker='.', markersize=5.0, linewidth=1.2, linestyle='-', label=None, alpha=1.0, color=color,
                                        path_effects=[pe.SimpleLineShadow(offset=(0.6, -0.6), alpha=0.2), pe.Normal()])
                        ax.fill_between(x=xdata, y1=ydata - edata, y2=ydata + edata, alpha=0.15, label=None, color=color)

                        l1_legend['handle'].append(line)
                        l1_legend['label'].append("${}={}$".format(l1[0], db_vals[l1[0]]))
                ax.set_title("${}={}$, ${}={}$ ".format(sub2[0], db_vals[sub2[0]], sub2[1], db_vals[sub2[1]]))
                ax.set_xlabel("${}$".format(x1[0]))
                ax.set_ylabel("$S(L/2)$")
                ax.xaxis.set_tick_params(labelbottom=True)
                ax.yaxis.set_tick_params(labelleft=True)
            fig.suptitle('{} entropy: ${}={}$,  ${}={}$, {} realizations'.format(
                meta['titlename'],
                fig2[0], db_vals[fig2[0]],
                fig2[1], db_vals[fig2[1]], meta['realizations']))
            prettify_plot2(fig=fig, axes=axes, cols=figcols, rows=figrows, axes_used=axes_used, extra_legend=[l1_legend])
            suffix = '_normpage' if meta['normpage'] else ''
            plt.savefig(
                "{}/{}_tinfty({})_fig({}_{})_sub({}_{}){}.pdf".format(meta['plotdir'], meta['plotprefix'], x1[0], str(key0), str(key1), sub2[0], sub2[1],
                                                                      suffix), format='pdf')


def plot_Stinfty_fig2_sub1_line2(db, meta, fig2, sub1, l2, x1, algo_filter=None, state_filter=None, point_filter=None):
    if len(fig2) != 2:
        raise AssertionError("fix must have length 2")
    if len(sub1) != 1:
        raise AssertionError("sub must have length 2")
    if len(l2) != 2:
        raise AssertionError("itr must have length 1")

    for key0 in db['keys'][fig2[0]]:
        for key1 in db['keys'][fig2[1]]:
            figrows, figcols = get_optimal_subplot_num(1 + len(db['keys'][sub1[0]]))
            fig, axes = plt.subplots(nrows=figrows, ncols=figcols, figsize=(3.5 * figcols, 3.5 * figrows), sharey='all', sharex='all')
            fig.tight_layout(pad=5, w_pad=1.0, h_pad=1.0)
            fig.subplots_adjust(wspace=0.25, hspace=0.40)
            axes_used = []
            l1_legend = {'handle': [], 'label': [], 'ncol': 1, 'unique': True, 'loc': 'upper left'}
            l2_legend = {'handle': [], 'label': [], 'ncol': 1, 'unique': True, 'loc': 'upper right'}

            for idx, (key2, ax) in enumerate(zip(db['keys'][sub1[0]], np.ravel(axes))):
                print(idx, key2)
                axes_used.append(idx)
                for algokey, statekey, pointkey in product(db['keys']['algo'], db['keys']['state'], db['keys']['point']):
                    palette = sns.color_palette("husl", len(db['keys'][l2[0]]))
                    for key3, color in zip(db['keys'][l2[0]], palette):
                        lstyles = itertools.cycle(['-', ':', '--', '-.'])
                        for key4, lstyle in zip(db['keys'][l2[1]], lstyles):
                            ydata = []
                            xdata = []
                            edata = []
                            for key5 in db['keys'][x1[0]]:
                                findlist = [key0, key1, key2, key3, key4, key5, algokey, statekey, pointkey, meta['dsetname']]
                                datanode = [value['datanode'] for key, value in db['dsets'].items() if all(k in key for k in findlist)]
                                if len(datanode) != 1:
                                    print("found", len(datanode), "datanodes: ", datanode, " | findlist: ", findlist)
                                    continue
                                    raise LookupError("Found incorrect number of datanodes")

                                datanode = datanode[0]
                                mmntnode = datanode.parent
                                db_vals = db['dsets'][datanode.name]
                                x = db_vals[x1[0]]
                                cachestr = '{}_tinfty'.format(meta['plotprefix'])
                                if cachestr in db_vals:
                                    y, e, i = db_vals[cachestr]
                                else:
                                    y, e, i = gather_at_tinfty(mmntnode, yname=meta['dsetname'])
                                    db_vals[cachestr] = y, e, i
                                ydata.append(y)
                                xdata.append(x)
                                edata.append(e)

                            ydata = np.asarray(ydata)
                            edata = np.asarray(edata)
                            xdata = np.asarray(xdata)
                            if meta['normpage']:
                                ydata = ydata / page_entropy(db_vals['L'])

                            line, = ax.plot(xdata, ydata, marker='.', markersize=5.0, linewidth=1.2, linestyle=lstyle, label=None, alpha=1.0, color=color,
                                            path_effects=[pe.SimpleLineShadow(offset=(0.6, -0.6), alpha=0.2), pe.Normal()])
                            ax.fill_between(x=xdata, y1=ydata - edata, y2=ydata + edata, alpha=0.15, label=None, color=color)

                            l1_legend['handle'].append(line)
                            l1_legend['label'].append("${}={}$".format(l2[0], db_vals[l2[0]]))
                            l2_legend['handle'].append(line)
                            l2_legend['label'].append("${}={}$".format(l2[1], db_vals[l2[1]]))

                ax.set_title("${}={}$ ".format(sub1[0], db_vals[sub1[0]]))
                ax.set_xlabel("${}$".format(x1[0]))
                ax.set_ylabel("$S(L/2)$")
                ax.xaxis.set_tick_params(labelbottom=True)
                ax.yaxis.set_tick_params(labelleft=True)

            fig.suptitle('{} entropy: ${}={}$,  ${}={}$, {} realizations'.format(
                meta['titlename'],
                fig2[0], db_vals[fig2[0]],
                fig2[1], db_vals[fig2[1]], meta['realizations']))
            prettify_plot2(fig=fig, axes=axes, cols=figcols, rows=figrows, axes_used=axes_used, extra_legend=[l1_legend, l2_legend])
            suffix = '_normpage' if meta['normpage'] else ''
            plt.savefig("{}/{}_tinfty({})_fig({}_{})_sub({}){}.pdf".format(meta['plotdir'], meta['plotprefix'], x1[0], str(key0), str(key1), sub1[0], suffix),
                        format='pdf')


def plot_Stinfty_fig1_sub1_line3(db, meta, fig1, sub1, l3, x1, algo_filter=None, state_filter=None, point_filter=None):
    if len(fig1) != 1:
        raise AssertionError("fix must have length 2")
    if len(sub1) != 1:
        raise AssertionError("sub must have length 2")
    if len(l3) != 3:
        raise AssertionError("itr must have length 1")

    for key0 in db['keys'][fig1[0]]:
        figrows, figcols = get_optimal_subplot_num(1 + len(db['keys'][sub1[0]]))
        fig, axes = plt.subplots(nrows=figrows, ncols=figcols, figsize=(3.5 * figcols, 3.5 * figrows), sharey='all', sharex='all')
        fig.tight_layout(pad=5, w_pad=1.0, h_pad=1.0)
        fig.subplots_adjust(wspace=0.25, hspace=0.40)
        axes_used = []
        l1_legend = {'handle': [], 'label': [], 'ncol': 1, 'unique': True, 'loc': 'upper left', 'color': None}
        l2_legend = {'handle': [], 'label': [], 'ncol': 1, 'unique': True, 'loc': 'upper right', 'color': 'dimgray'}
        l3_legend = {'handle': [], 'label': [], 'ncol': 1, 'unique': True, 'loc': 'lower right', 'color': 'dimgray'}

        for idx, (key1, ax) in enumerate(zip(db['keys'][sub1[0]], np.ravel(axes))):
            print(idx, key0, key1)
            axes_used.append(idx)
            for algokey, statekey, pointkey in product(db['keys']['algo'], db['keys']['state'], db['keys']['point']):
                palette = sns.color_palette("husl", len(db['keys'][l3[0]]))
                for key2, color in zip(db['keys'][l3[0]], palette):
                    mstyles = itertools.cycle(['o', 's', '*', 'h', '+', 'v', '^', '<', '>'])
                    for key3, mstyle in zip(db['keys'][l3[1]], mstyles):
                        lstyles = itertools.cycle(['-', ':', '--', '-.'])
                        for key4, lstyle in zip(db['keys'][l3[2]], lstyles):
                            ydata = []
                            xdata = []
                            edata = []
                            for key5 in db['keys'][x1[0]]:
                                findlist = [key0, key1, key2, key3, key4, key5, algokey, statekey, pointkey, meta['dsetname']]
                                datanode = [value['datanode'] for key, value in db['dsets'].items() if all(k in key for k in findlist)]
                                if len(datanode) != 1:
                                    print("found", len(datanode), "datanodes: ", datanode, " | findlist: ", findlist)
                                    continue
                                    raise LookupError("Found incorrect number of datanodes")

                                datanode = datanode[0]
                                mmntnode = datanode.parent
                                db_vals = db['dsets'][datanode.name]
                                x = db_vals[x1[0]]
                                cachestr = '{}_tinfty'.format(meta['plotprefix'])
                                if cachestr in db_vals:
                                    y, e, i = db_vals[cachestr]
                                else:
                                    y, e, i = gather_at_tinfty(mmntnode, yname=meta['dsetname'])
                                    db_vals[cachestr] = y, e, i
                                ydata.append(y)
                                xdata.append(x)
                                edata.append(e)

                            ydata = np.asarray(ydata)
                            edata = np.asarray(edata)
                            xdata = np.asarray(xdata)
                            if meta['normpage']:
                                ydata = ydata / page_entropy(db_vals['L'])

                            line, = ax.plot(xdata, ydata, markersize=4.0, linewidth=1.1, marker=mstyle, linestyle=lstyle, label=None, alpha=1.0, color=color,
                                            path_effects=[pe.SimpleLineShadow(offset=(0.6, -0.6), alpha=0.2), pe.Normal()])
                            ax.fill_between(x=xdata, y1=ydata - edata, y2=ydata + edata, alpha=0.15, label=None, color=color)

                            l1_legend['handle'].append(line)
                            l1_legend['label'].append("${}={}$".format(l3[0], db_vals[l3[0]]))
                            l2_legend['handle'].append(
                                plt.Line2D([0], [0], color=l2_legend['color'], markersize=4.0, marker=mstyle, linewidth=None, linestyle="None"))
                            l2_legend['label'].append("${}={}$".format(l3[1], db_vals[l3[1]]))
                            l3_legend['handle'].append(
                                plt.Line2D([0], [0], color=l3_legend['color'], markersize=None, marker=None, linewidth=1.1, linestyle=lstyle))
                            l3_legend['label'].append("${}={}$".format(l3[2], db_vals[l3[2]]))

            ax.set_title("${}={}$ ".format(sub1[0], db_vals[sub1[0]]))
            ax.set_xlabel("${}$".format(x1[0]))
            ax.set_ylabel("$S(L/2)$")
            ax.xaxis.set_tick_params(labelbottom=True)
            ax.yaxis.set_tick_params(labelleft=True)

        fig.suptitle('{} entropy: ${}={}$, {} realizations'.format(meta['titlename'], fig1[0], db_vals[fig1[0]], meta['realizations']))
        prettify_plot2(fig=fig, axes=axes, cols=figcols, rows=figrows, axes_used=axes_used, extra_legend=[l1_legend, l2_legend, l3_legend])
        suffix = '_normpage' if meta['normpage'] else ''
        plt.savefig("{}/{}_tinfty({})_fig({})_sub({}){}.pdf".format(meta['plotdir'], meta['plotprefix'], x1[0], str(key0), sub1[0], suffix), format='pdf')
