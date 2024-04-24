from database.database import *
from plotting.tools import *
from itertools import product
from copy import deepcopy
import numpy as np
import seaborn as sns
import matplotlib.patheffects as pe
from pathlib import Path


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


@njit(parallel=True, cache=True)
def gather_at_tinfty_v2(ydata, findsaturation=True):
    if len(ydata) == 0:
        raise LookupError("ydata is empty")

    if findsaturation:
        # Instead of returning the whole time series, we want the average saturation value.
        # So first we need to find the index where y saturates
        idx = find_saturation_idx(ydata, np.max(ydata) / 10)
    else:
        idx = len(ydata) - 1

    y = ydata[idx:]
    return np.mean(y), np.std(y) / np.sqrt(len(y)), idx


def plot_v2_at_tinfty_fig3_sub2_line1(db, meta, fig3, sub2, l1, x1, algo_filter=None, state_filter=None, point_filter=None, f=None, palette_name=None):
    if len(fig3) != 3:
        raise AssertionError("fig must have length 3")
    if len(sub2) != 2:
        raise AssertionError("sub must have length 2")
    if len(l1) != 1:
        raise AssertionError("l must have length 1")
    if len(x1) != 1:
        raise AssertionError("x must have length 1")

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

    legend_col_keys = l1.copy()
    if legendcols := meta['legendcols']:
        for col in legendcols:
            if not col in [l.split(':')[0] for l in sub2 + l1]:
                legend_col_keys.append(col)

    for key0 in db['keys'][fig3[0]]:
        for key1 in db['keys'][fig3[1]]:
            for key2 in db['keys'][fig3[2]]:
                keyprod = list(product(*get_keys(db, sub2)))
                numplots = len(keyprod)
                if not f:
                    f = get_fig_meta(numplots, meta=meta)

                for idx, ((key3, key4), ax) in enumerate(zip(keyprod, f['ax'])):
                    for algokey, statekey, cronokey in product(db['keys']['algo'], db['keys']['state'], db['keys']['crono']):
                        palette = sns.color_palette(palette=palette_name, n_colors=len(db['keys'][l1[0]]))
                        for idx5, (key5, color) in enumerate(zip(get_keys(db, l1[0]), palette)):
                            ydata = []
                            xdata = []
                            edata = []
                            for key6 in get_keys(db, x1[0]):
                                findlist = [key0, key1, key2, key3, key4, key5, key6, algokey, statekey, cronokey, meta['groupname']]
                                datanode = [value['datanode'] for key, value in db['dsets'].items() if all(k in key for k in findlist)]
                                if len(datanode) != 1:
                                    print("found", len(datanode), "datanodes: ", datanode, " | findlist: ", findlist)
                                    continue

                                    # raise LookupError("Found incorrect number of datanodes")
                                print('datanode:', datanode)

                                datanode = datanode[0]
                                dbval = db['dsets'][datanode.name]
                                dbtex = db['dsets'][datanode.name]['tex']
                                x = dbval[x1[0]]
                                cachestr = '{}_tinfty'.format(meta['plotprefix'])
                                if cachestr in dbval:
                                    y, e, i = dbval[cachestr]
                                else:
                                    y, e, i = gather_at_tinfty_v2(datanode['avg'][meta['colname']][()], meta['findsaturation'])
                                    dbval[cachestr] = y, e, i
                                ydata.append(y)
                                xdata.append(x)
                                edata.append(e)

                            ydata = np.asarray(ydata)
                            edata = np.asarray(edata)
                            xdata = np.asarray(xdata)

                            if 'normalize' in meta:
                                ydata = ydata / meta['normalize']
                                edata = edata / meta['normalize']
                            ax.fill_between(x=xdata, y1=ydata - edata, y2=ydata + edata, alpha=0.15, label=None, color=color)
                            line, = ax.plot(xdata, ydata, marker='.', linestyle='-', label=None, color=color, path_effects=path_effects)

                            label_str = '{}'.format(db['vals'][l1[0]][idx5])
                            f['legends'][idx][0]['handle'].append(line)
                            f['legends'][idx][0]['label'].append(label_str)
                            f['legends'][idx][0]['title'] = db['tex'][l1[0]]

                    f['axes_used'].append(idx)
                    ax.set_title(get_title(dbval, sub2), x=0, horizontalalignment='left')
                    ax.set_xlabel(db['tex'][x1[0]])
                prettify_plot5(fmeta=f)

                suffix = ''
                suffix = suffix + '_normpage' if 'normpage' in meta and meta['normpage'] else suffix
                suffix = suffix + '_loglog' if meta.get('timeselection') == 'lnlnt' else suffix
                f['fig'].savefig(
                    "{}/{}_tinfty({})_fig({}_{}_{})_sub({}_{}){}.pdf"
                    .format(meta['plotdir'], meta['plotprefix'],
                            x1[0], str(key0), str(key1), str(key2), sub2[0], sub2[1], suffix), format='pdf')
                f['fig'].savefig(
                    "{}/{}_tinfty({})_fig({}_{}_{})_sub({}_{}).png"
                    .format(meta['plotdir'], meta['plotprefix'],
                            x1[0], str(key0), str(key1), str(key2), sub2[0], sub2[1], suffix), format='png')


def plot_v2_at_tinfty_fig2_sub2_line2(db, meta, fig2, sub2, l2, x1, algo_filter=None, state_filter=None, point_filter=None):
    if len(fig2) != 2:
        raise AssertionError("fix must have length 2")
    if len(sub2) != 2:
        raise AssertionError("sub must have length 2")
    if len(l2) != 2:
        raise AssertionError("itr must have length 2")
    if len(x1) != 1:
        raise AssertionError("x must have length 1")

    if 'mplstyle' in meta:
        plt.style.use(meta['mplstyle'])
    if 'plotdir' in meta and 'mplstyle' in meta:
        meta['plotdir'] = Path(meta['plotdir'], Path(meta['mplstyle']).stem)
        Path(meta['plotdir']).mkdir(parents=True, exist_ok=True)
        print("Setting plotdir: ", meta['plotdir'])
    if 'mplstyle' in meta and 'slack' in meta['mplstyle']:
        palette_name = "Spectral"
        path_effects = [pe.SimpleLineShadow(offset=(0.6, -0.6), alpha=0.2), pe.Normal()]
    else:
        palette_name = "colorblind"
        path_effects = None

    prb_style = 'prb' in meta['mplstyle'] if 'mplstyle' in meta else False

    legend_col_keys = l2
    if 'legendcols' in meta:
        legend_col_keys.extend([col for col in meta['legendcols']])

    for key0 in db['keys'][fig2[0]]:
        for key1 in db['keys'][fig2[1]]:
            numplots = len(db['keys'][sub2[0]]) * len(db['keys'][sub2[1]])
            f = get_fig_meta(numplots, meta=meta)

            for idx, ((key2, key3), ax) in enumerate(
                    zip(product(db['keys'][sub2[0]], db['keys'][sub2[1]]), f['ax'])):
                for algokey, statekey, cronokey in product(db['keys']['algo'], db['keys']['state'], db['keys']['crono']):
                    palette_names = get_uniform_palette_names(len(db['keys'][l2[0]]))
                    mstyles = itertools.cycle(['o', 's', '*', 'h', '+', 'v', '^', '<', '>'])
                    for idx4, (key4, palette_name, mstyle) in enumerate(zip(db['keys'][l2[0]], palette_names, mstyles)):
                        palette = sns.color_palette(palette_name, len(db['keys'][l2[1]]))
                        # lstyles = itertools.cycle(['-', ':', '--', '-.'])
                        for idx5, (key5, color) in enumerate(zip(db['keys'][l2[1]], palette)):
                            ydata = []
                            xdata = []
                            edata = []
                            for key6 in db['keys'][x1[0]]:
                                findlist = [key0, key1, key2, key3, key4, key5, key6, algokey, statekey, cronokey, meta['groupname']]
                                datanode = [value['datanode'] for key, value in db['dsets'].items() if all(k in key for k in findlist)]
                                if len(datanode) != 1:
                                    print("found", len(datanode), "datanodes: ", datanode, " | findlist: ", findlist)
                                    continue
                                    raise LookupError("Found incorrect number of datanodes")

                                datanode = datanode[0]
                                dbval = db['dsets'][datanode.name]
                                x = dbval[x1[0]]
                                cachestr = '{}_tinfty'.format(meta['plotprefix'])
                                if cachestr in dbval:
                                    y, e, i = dbval[cachestr]
                                else:
                                    y, e, i = gather_at_tinfty_v2(datanode['avg'][meta['colname']][()], meta['findsaturation'])
                                    dbval[cachestr] = y, e, i
                                ydata.append(y)
                                xdata.append(x)
                                edata.append(e)
                            ydata = np.asarray(ydata)
                            edata = np.asarray(edata)
                            xdata = np.asarray(xdata)
                            if meta['normpage']:
                                ydata = ydata / page_entropy(dbval['L'])

                            line, = ax.plot(xdata, ydata, markersize=6.0, linewidth=1.1, marker=mstyle, markerfacecolor='none', linestyle='-', label=None,
                                            alpha=1.0, color=color,
                                            path_effects=[pe.SimpleLineShadow(offset=(0.6, -0.6), alpha=0.2), pe.Normal()])
                            ax.fill_between(x=xdata, y1=ydata - edata, y2=ydata + edata, alpha=0.15, label=None, color=color)

                            f['legends'][idx][0]['handle'].append(line)
                            f['legends'][idx][0]['label'].append(db['vals'][l2[0]][idx4])
                            f['legends'][idx][0]['title'] = db['tex'][l2[0]]

                            f['legends'][idx][1]['handle'].append(line)
                            f['legends'][idx][1]['label'].append(db['vals'][l2[1]][idx5])
                            f['legends'][idx][1]['title'] = db['tex'][l2[1]]

                            # legendrow = get_legend_row(dbval=dbval, avgnode=datanode['avg'], meta=meta,
                            #                            legend_col_keys=legend_col_keys)
                            # for icol, (col, key) in enumerate(zip(legendrow, legend_col_keys)):
                            #     f['legends'][idx][icol]['handle'].append(line)
                            #     f['legends'][idx][icol]['label'].append(col)
                            #     f['legends'][idx][icol]['title'] = db['tex'][key]

                            # l['handle'].append(line)
                            # l['label'].append([])
                            # for key in l2:
                            #     l['label'][-1].extend(['{}'.format(dbval[key])])

                f['axes_used'].append(idx)
                ax.set_title(get_title(dbval, sub2), x=0, horizontalalignment='left')
                ax.set_xlabel(db['tex'][x1[0]])
            prettify_plot5(fmeta=f)
            suffix = ''
            suffix = suffix + '_normpage' if 'normpage' in meta and meta['normpage'] else suffix
            # suffix = suffix + '_loglog' if meta.get('timeselection') == 'lnlnt' else suffix
            plt.savefig(
                "{}/{}_tinfty({})_fig({}_{})_sub({}_{}){}.pdf".format(meta['plotdir'], meta['plotprefix'], x1[0], str(key0),
                                                                      str(key1), sub2[0], sub2[1], suffix),
                bbox_inches="tight", format='pdf')
            plt.savefig(
                "{}/{}_tinfty({})_fig({}_{})_sub({}_{}){}.png".format(meta['plotdir'], meta['plotprefix'], x1[0], str(key0),
                                                                      str(key1), sub2[0], sub2[1], suffix),
                bbox_inches="tight", format='png')


def plot_v2_at_tinfty_fig2_sub1_line3(db, meta, fig2, sub1, l3, x1, algo_filter=None, state_filter=None, point_filter=None):
    if len(fig2) != 2:
        raise AssertionError("fix must have length 2")
    if len(sub1) != 1:
        raise AssertionError("sub must have length 1")
    if len(l3) != 3:
        raise AssertionError("itr must have length 3")
    if len(x1) != 1:
        raise AssertionError("x must have length 1")

    if 'mplstyle' in meta:
        plt.style.use(meta['mplstyle'])

    if 'plotdir' in meta and 'mplstyle' in meta:
        meta['plotdir'] = Path(meta['plotdir'], Path(meta['mplstyle']).stem)
        Path(meta['plotdir']).mkdir(parents=True, exist_ok=True)
        print("Setting plotdir: ", meta['plotdir'])
    if 'mplstyle' in meta and 'slack' in meta['mplstyle']:
        palette_name = "Spectral"
        path_effects = [pe.SimpleLineShadow(offset=(0.6, -0.6), alpha=0.2), pe.Normal()]
    else:
        palette_name = "colorblind"
        path_effects = None

    prb_style = 'prb' in meta['mplstyle'] if 'mplstyle' in meta else False

    legend_col_keys = l3
    if 'legendcols' in meta:
        legend_col_keys.extend([col for col in meta['legendcols']])

    for key0 in db['keys'][fig2[0]]:
        for key1 in db['keys'][fig2[1]]:
            numplots = len(db['keys'][sub1[0]])
            f = get_fig_meta(numplots, meta=meta)

            for idx, (key2, ax) in enumerate(zip(db['keys'][sub1[0]], f['ax'])):
                for algokey, statekey, cronokey in product(db['keys']['algo'], db['keys']['state'], db['keys']['crono']):
                    mstyles = itertools.cycle(['o', 's', '*', 'h', '+', 'v', '^', '<', '>'])
                    # lstyles = itertools.cycle(['-', ':', '--', '-.'])
                    for idx3, (key3, mstyle) in enumerate(zip(db['keys'][l3[0]], mstyles)):
                        palette_names = get_uniform_palette_names(len(db['keys'][l3[1]]))
                        for idx4, (key4, palette_name) in enumerate(zip(db['keys'][l3[1]], palette_names)):
                            palette = sns.color_palette(palette_name, len(db['keys'][l3[2]]))
                            for idx5, (key5, color) in enumerate(zip(db['keys'][l3[2]], palette)):
                                ydata = []
                                xdata = []
                                edata = []
                                for key6 in db['keys'][x1[0]]:
                                    findlist = [key0, key1, key2, key3, key4, key5, key6, algokey, statekey, cronokey, meta['groupname']]
                                    datanode = [value['datanode'] for key, value in db['dsets'].items() if all(k in key for k in findlist)]
                                    if len(datanode) != 1:
                                        print("found", len(datanode), "datanodes: ", datanode, " | findlist: ", findlist)
                                        continue
                                        raise LookupError("Found incorrect number of datanodes")

                                    datanode = datanode[0]
                                    dbval = db['dsets'][datanode.name]
                                    x = dbval[x1[0]]
                                    cachestr = '{}_tinfty'.format(meta['plotprefix'])
                                    if cachestr in dbval:
                                        y, e, i = dbval[cachestr]
                                    else:
                                        y, e, i = gather_at_tinfty_v2(datanode['avg'][meta['colname']][()], meta['findsaturation'])
                                        dbval[cachestr] = y, e, i
                                    ydata.append(y)
                                    xdata.append(x)
                                    edata.append(e)
                                ydata = np.asarray(ydata)
                                edata = np.asarray(edata)
                                xdata = np.asarray(xdata)
                                if meta['normpage']:
                                    ydata = ydata / page_entropy(dbval['L'])

                                line, = ax.plot(xdata, ydata, marker=mstyle, markerfacecolor='none', linestyle='-', label=None, alpha=1.0, color=color,
                                                path_effects=path_effects)
                                ax.fill_between(x=xdata, y1=ydata - edata, y2=ydata + edata, alpha=0.15, label=None, color=color)

                                legendrow = [db['vals'][l3[0]][idx3], db['vals'][l3[1]][idx4], db['vals'][l3[2]][idx5]]
                                for icol, (col, key) in enumerate(zip(legendrow, legend_col_keys)):
                                    f['legends'][idx][icol]['handle'].append(line)
                                    f['legends'][idx][icol]['label'].append(str(col))
                                    f['legends'][idx][icol]['title'] = db['tex'][key]

                                # f['legends'][idx][0]['handle'].append(line)
                                # f['legends'][idx][0]['title'] = [db['tex'][l3[0]], db['tex'][l3[1]], db['tex'][l3[2]]]
                                # f['legends'][idx][0]['label'].append(db['vals'][l3[0]][idx3])
                                # f['legends'][idx][0]['label'].append(db['vals'][l3[1]][idx4])
                                # f['legends'][idx][0]['label'].append(db['vals'][l3[2]][idx5])

                                # legendrow = get_legend_row(dbval=dbval, avgnode=datanode['avg'], meta=meta,
                                #                            legend_col_keys=legend_col_keys)
                                # for icol, (col, key) in enumerate(zip(legendrow, legend_col_keys)):
                                #     f['legends'][idx][icol]['handle'].append(line)
                                #     f['legends'][idx][icol]['label'].append(col)
                                #     f['legends'][idx][icol]['title'] = db['tex'][key]
                                #

                                # l['handle'].append(line)
                                # l['label'].append([])
                                # for key in l3:
                                #     l['label'][-1].extend(['{}'.format(dbval[key])])

                f['axes_used'].append(idx)
                ax.set_title(get_title(dbval, sub1), x=0, horizontalalignment='left')
                ax.set_xlabel(db['tex'][x1[0]])
            prettify_plot5(fmeta=f)

            suffix = ''
            suffix = suffix + '_normpage' if 'normpage' in meta and meta['normpage'] else suffix
            # suffix = suffix + '_loglog' if meta.get('timeselection') == 'lnlnt' else suffix
            plt.savefig(
                "{}/{}_tinfty({})_fig({}_{})_sub({}){}.pdf".format(meta['plotdir'], meta['plotprefix'], x1[0], str(key0),
                                                                   str(key1), sub1[0], suffix), format='pdf')
            plt.savefig(
                "{}/{}_tinfty({})_fig({}_{})_sub({}){}.png".format(meta['plotdir'], meta['plotprefix'], x1[0], str(key0),
                                                                   str(key1), sub1[0], suffix), format='png')
