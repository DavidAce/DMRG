from src.plotting.tools import *
import matplotlib.pyplot as plt
from src.database.database import *
from src.general.filter import *
from matplotlib.ticker import MaxNLocator
import itertools


def multiplot_Smid_sub1_l1_x1(h5_avg, db, meta, plotdir='', g3=['algo', 'state', 'point'], algo_inc='', state_inc='', sub1=['L'], l1=['l'], x1=['d'],
                              vwin=None):
    print('Plotting: S vs Site for: ', algo_inc, state_inc)
    if len(sub1) != 1:
        raise AssertionError("sub2 must have length 1")
    if len(l1) != 1:
        raise AssertionError("l1 must have length 1")
    if len(x1) != 1:
        raise AssertionError("x1 must have length 1")
    # One figure per L
    # One subplot per delta
    # One line per lambda per state
    # One figure per L
    # One subplot per delta
    # One line per lambda per state
    numsub = len(db['keys'][sub1[0]])
    figrows, figcols = get_optimal_subplot_num(numsub)
    fig, axes = plt.subplots(nrows=figrows, ncols=figcols, figsize=(5 * figcols, 5 * figrows), sharey='all')
    fig.tight_layout(pad=5, w_pad=1.0, h_pad=1.0)
    fig.subplots_adjust(wspace=0.2, hspace=0.2)
    axes_used = []
    for subidx, (key0, ax) in enumerate(zip(db['keys'][sub1[0]], np.ravel(axes))):
        ed_palette = itertools.cycle(sns.color_palette("Set2"))
        numcolors = len(db['keys'][g3[0]]) * len(db['keys'][g3[1]]) * len(db['keys'][g3[2]])
        current_palette = itertools.cycle(sns.color_palette("colorblind", numcolors))
        lstyles = itertools.cycle(['-.', '-', '--', ':', ])
        mstyles = itertools.cycle(('.', ',', '+', 'o', '*'))
        delt = None
        size = None
        for lidx, key1 in enumerate(db['keys'][l1[0]]):
            for gidx, (key2, key3, key4) in enumerate(product(db['keys'][g3[0]], db['keys'][g3[1]], db['keys'][g3[2]])):
                ydata = []
                xdata = []
                edata = []
                ndata = []
                for key5 in db['keys'][x1[0]]:
                    findlist = [key0, key1, key2, key3, key4, key5, meta['dsetname']]
                    datanode = [value['datanode'] for key, value in db['dsets'].items() if all(k in key for k in findlist)]
                    if len(datanode) != 1:
                        continue
                    if not any(k in datanode[0].name for k in algo_inc):
                        continue
                    if not any(k in datanode[0].name for k in state_inc):
                        continue

                        # raise LookupError("Found incorrect number of datanodes")
                    datanode = datanode[0]
                    dbval = db['dsets'][datanode.name]
                    midx = dbval['midx']
                    lamb = dbval['l']
                    delt = dbval['d']
                    size = dbval['L']
                    ydata.append(get_data(datanode['avg'], 'L_', 'f8')[midx])
                    edata.append(get_data(datanode['ste'], 'L_', 'f8')[midx])
                    ndata.append(get_data(datanode['num'], 'L_', 'f8')[()])
                    xdata.append(dbval[x1[0]])
                    print("ndata: {}".format(ndata))
                if not xdata or not ydata:
                    continue
                sort = np.argsort(xdata)
                ydata = np.array(ydata)[sort]
                edata = np.array(edata)[sort]
                xdata = np.array(xdata)[sort]
                print("ydata after:", np.shape(ydata))
                print("edata after:", np.shape(edata))
                print("xdata after:", np.shape(xdata))
                print("ndata after:", np.shape(ndata))
                ndata = np.min(ndata)

                statekey = dbval['keys']['state']
                algokey = dbval['keys']['algo']

                style = dbval['style']

                lstyle = next(lstyles)
                mstyle = next(mstyles)

                if vwin:
                    ydata, xdata, edata, ndata = get_v_filtered_edata(ydata, xdata, edata, ndata, dbval, vwin)

                if "states" in statekey:
                    color = next(ed_palette)
                    mstyle = None
                    lstyle = 'solid'
                    lgnd = "ED ({})".format(ndata)
                    # nicename = "ED e=[" + statenode.attrs["efmt"] + "]"
                else:
                    color = next(current_palette)
                    lgnd = "{} {} ({})".format(re.sub(r'[\W_]', ' ', "{} {}".format(algokey, statekey)), db['tex'][key1], ndata)
                print("Plotting xdata {}, legend {}".format(np.shape(xdata), np.shape(lgnd)))
                ax.errorbar(x=xdata, y=ydata, yerr=edata, label=lgnd, capsize=2,
                            color=color, elinewidth=0.3, markeredgewidth=0.8,
                            marker=mstyle, markersize=4.0, linestyle=lstyle,
                            linewidth=style['lwidth'], alpha=style['lalpha'])

        ax.set_xlabel(db['tex'][x1[0]])
        ax.set_ylabel(meta['ylabel'] + "(L/2)")
        ax.set_title('{}'.format(db['tex'][key0]))
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        ax.legend(loc='lower center', fontsize='small', labelspacing=0.25)
        fig.suptitle('Entanglement entropy vs {} @ {}'.format(db['tex'][x1[0]], db['tex'][key0]))

        # $\Delta = $' + str(delt) + '$\lambda = $' + str(lamb))
        axes_used.append(subidx) if not subidx in axes_used else axes_used

        remove_empty_subplots(fig=fig, axes=axes, axes_used=axes_used)
        if plotdir != '':
            plt.savefig('{}/Smid_vs_Delta_{}.pdf'.format(plotdir, key0), format='pdf')
            plt.savefig('{}/Smid_vs_Delta_{}.png'.format(plotdir, key0), format='png')

    return

    # One figure per unique_l, unique_J and unique_h
    for l in path_l:
        for d in path_d:
            # In each figure we want one subplot per unique_L
            rows, cols = get_optimal_subplot_num(len(path_L))
            fig, axes = plt.subplots(nrows=rows, ncols=cols, figsize=(7 * cols, 7 * rows))
            fig.tight_layout(pad=5, w_pad=1.0, h_pad=1.0)
            fig.subplots_adjust(wspace=0.3, hspace=0.3)
            used_ax = 0
            delt = 0
            lamb = 0
            max_S = 0
            for ax, L in zip(np.ravel(axes), path_L):
                ed_palette = itertools.cycle(sns.color_palette("Set2"))
                current_palette = itertools.cycle(sns.color_palette())
                basenode = h5_src[L][l][d]
                size = basenode.attrs['L']
                delt = basenode.attrs['d']
                lamb = basenode.attrs['l']
                for algokey, algopath, algonode in h5py_group_iterator(g=basenode, filter=algo_inc, dep=1):
                    for statekey, statepath, statenode in h5py_group_iterator(g=algonode, filter=state_inc,
                                                                              dep=1):
                        for v_win_idx, v_win in enumerate(variance_window_limits):
                            for e_win_idx, e_win in enumerate(energy_window_limits):
                                idx = get_v_e_filtered_index_list(statenode, v_win, e_win)
                                for datakey, datapath, datanode in h5py_node_finder(g=statenode,
                                                                                    filter='entanglement_entropies',
                                                                                    dep=8):
                                    ndata = datanode['num'][()]
                                    ydata = np.array(datanode['avg'])
                                    edata = np.array(datanode['ste'])
                                    xdata = range(len(ydata))

                                    if np.any(np.isnan(ydata)):
                                        raise ValueError("Data contains nan's")
                                    if np.any(np.isnan(ydata)):
                                        raise ValueError("Standard error contains nan's")

                                    if "states" in statekey:
                                        color = next(ed_palette)
                                        nicename = "ED e=[" + statenode.attrs["efmt"] + "]"
                                        lwidth = 3.0
                                        lalpha = 0.9
                                        mstyle = None
                                        lstyle = 'solid'
                                        # Test taking a subset of states
                                        # edata = bootstrap_sterr(data=datanode['data'],chunksize=2000)

                                        # idx = random.sample(range(0,ndata-1),2000)
                                        # idx.sort()
                                        # ydata = np.nanmean(np.array(datanode['data'][:,idx]),axis=1)


                                    else:
                                        color = next(current_palette)
                                        nicename = re.sub(r'[\W_]', ' ', str(algokey + " " + statekey))
                                        lwidth = 1.4
                                        lalpha = 1.0
                                        mstyle = '.'
                                        lstyle = 'dotted'

                                    nicename = nicename + ' (' + str(ndata) + ')'
                                    ax.errorbar(x=xdata, y=ydata, yerr=edata, label=nicename, capsize=2,
                                                color=color,
                                                elinewidth=0.3, markeredgewidth=0.8, marker=mstyle,
                                                linestyle=lstyle,
                                                linewidth=lwidth, alpha=lalpha)
                                    max_S = np.max([max_S, np.max(ydata)])

                                    # Select a random subset of the data
                                    if "states" in statekey:
                                        continue
                                    data = np.array(datanode['data'])
                                    ndata = int(data.shape[1] * 0.5)
                                    for rep in range(0):
                                        idx = np.random.choice(data.shape[1], ndata, replace=False)
                                        subdata = data[:, idx]
                                        std = np.nanstd(subdata, axis=1)
                                        ydata = np.nanmean(subdata, axis=1)
                                        edata = std / np.sqrt(ndata)
                                        nicename = re.sub(r'[\W_]', ' ',
                                                          str(algokey + " " + statekey)) + ' (random ' + str(
                                            ndata) + ')'
                                        ax.errorbar(x=xdata, y=ydata, yerr=edata, label=nicename, capsize=2,
                                                    elinewidth=0.3, markeredgewidth=0.8, linewidth=0.3)
                                        max_S = np.max([max_S, np.max(ydata)])

                ax.set_xlabel('Site $l$')
                ax.set_ylabel('$S_E(l)$')
                ax.set_title('$L = ' + str(size) + '$')
                used_ax = used_ax + 1
                ax.legend()

            for ax in np.ravel(axes)[used_ax:]:
                fig.delaxes(ax)
            for ax in np.ravel(axes):
                ax.set_ylim(ymin=0, ymax=max_S * 1.2)
            fig.suptitle('Entanglement entropy vs site @ $\Delta = ' + str(delt) + '\quad \lambda = ' + str(
                lamb) + '$')
            if plotdir != '':
                plt.savefig(plotdir + '/S_vs_Site_' + l + '_' + d + '.pdf', format='pdf')
                plt.savefig(plotdir + '/S_vs_Site_' + l + '_' + d + '.png', format='png')

    h5close(h5_src)


def multiplot_Smid_vs_Delta(h5_src, db=None, plotdir='', algo_inc='', state_inc=''):
    print('Plotting: Entanglement entropy distribution vs Delta for: ', algo_inc, state_inc)
    if not db:
        db = load_database(h5_src, 'entanglement_entropy_midchain', algo_inc, state_inc)

    for Lkey in db['keys']['L']:
        figrows, figcols = get_optimal_subplot_num(1 + len(db['keys']['l']))
        fig, axes = plt.subplots(nrows=figrows, ncols=figcols, figsize=(4.5 * figcols, 4.5 * figrows), sharey='all',
                                 sharex='all')
        fig.tight_layout(pad=5, w_pad=1.0, h_pad=1.0)
        fig.subplots_adjust(wspace=0.2, hspace=0.2)
        axes_used = []
        for lidx, (lkey, ax) in enumerate(zip(db['keys']['l'], np.ravel(axes))):
            ed_palette = itertools.cycle(sns.color_palette("Set2"))
            current_palette = itertools.cycle(sns.color_palette())
            lamb = None
            size = None
            nlabel = {'line': [], 'text': []}
            for algoidx, algokey in enumerate(db['keys']['algo']):
                for stateidx, statekey in enumerate(db['keys']['state']):
                    if not contains(algokey, algo_inc) or not contains(statekey, state_inc):
                        continue
                    dsetkeys = [x for x in db['dsets'] if
                                Lkey in x and lkey in x and algokey in x and statekey in x]
                    if not dsetkeys:
                        continue

                    xdata, ydata, edata, ndata = [], [], [], []
                    lgnd, style = None, None
                    color = next(ed_palette) if "states" in statekey else next(current_palette)
                    for idx, dsetkey in enumerate(dsetkeys):
                        meta = db['dsets'][dsetkey]
                        midx = meta['midx']
                        ndata.append(meta['num'])
                        xdata.append(meta['d'])
                        ydata.append(get_data(meta['datanode']['avg'], 'L_', 'f8')[midx])
                        edata.append(get_data(meta['datanode']['ste'], 'L_', 'f8')[midx])

                        lamb = lamb if lamb else meta['l']
                        size = size if size else meta['L']
                        style = style if style else meta['style']
                        lgnd = lgnd if lgnd else (
                            "ED" if "states" in statekey else
                            re.sub(r'[\W_]', ' ', str(algokey + " " + statekey)))

                    sortIdx = np.argsort(np.asarray(xdata))
                    xdata = np.asarray(xdata)[sortIdx]
                    ydata = np.asarray(ydata)[sortIdx]
                    edata = np.asarray(edata)[sortIdx]
                    line = ax.errorbar(x=xdata, y=ydata, yerr=edata, label=lgnd, color=color, capsize=2,
                                       elinewidth=0.3, markeredgewidth=0.8, marker=style['mstyle'],
                                       linestyle=style['lstyle'],
                                       linewidth=style['lwidth'], alpha=style['lalpha'])
                    nlabel['line'].append(line)
                    nlabel['text'].append('{}'.format(ndata))
                    # for i, n in enumerate(ndata):
                    #     min_S = db['min']['mvg']
                    #     max_S = db['max']['mvg']
                    #     ytext = min_S - 0.06 * (max_S - min_S) * i
                    #     # ax.annotate(txt, (val['x'][i], val['y'][i]), textcoords='data',
                    #     #             xytext=[val['x'][i] - 0.2, ytext],
                    #     #             color=val['color'], fontsize='x-small')
                    #     ax.annotate(n, xy=(xdata[i], ydata[i]), xytext=(xdata[i], ytext), color=color,
                    #                 alpha=0.8, fontsize='x-small')
                    axes_used.append(lidx) if not lidx in axes_used else axes_used

            ax.set_title('$\lambda = {:.4f}$'.format(lamb))
            ax.set_xlabel('$\Delta = \log \\bar J - \log \\bar h$')
            ax.set_ylabel('$S_E(L/2)$')
            # ax.legend(nlabel['line'], nlabel['text'], title='Realizations',
            #           loc='lower right', framealpha=0.2, fontsize='small', labelspacing=0.25, ncol=2)
            fig.suptitle('$L = {}$'.format(size))

        prettify_plot(fig, axes, rows=figrows, cols=figcols, axes_used=axes_used, ymin=db['min']['mvg'], ymax=db['max']['mvg'], nlabel=nlabel)
        if plotdir != '':
            plt.savefig(plotdir + '/Smid_vs_Delta_' + Lkey + '.pdf', format='pdf')
            plt.savefig(plotdir + '/Smid_vs_Delta_' + Lkey + '.png', format='png')

    # for dsetkey, dsetmeta in db['dsets'].items():
    #     k = dsetmeta['keys']
    #     for win_idx, win in enumerate(variance_window_limits):
    #         statenode = h5_src[k['L']][k['l']][k['d']][k['algo']][k['state']]
    #         datanode = h5_src[dsetkey]
    #         idx = get_v_inced_index_list(statenode, win)
    #         if not idx:
    #             data = datanode['data']
    #         else:
    #             data = datanode['data'][idx]
    #         if np.any(np.isnan(data)):
    #             raise ValueError("Data contains nan's")
    #         num = datanode['num'][()]
    #         avg = datanode['avg'][()]
    #         datarange = [0, db['max']]
    #
    #         hist, edge = np.histogram(data, bins=bins, range=datarange, density=False)
    #         bincentres = [(edge[i] + edge[i + 1]) / 2. for i in range(len(edge) - 1)]
    #         width = np.diff(edge)
    #         norm = np.dot(hist, width)
    #         hist = hist / norm
    #         dsetmeta['data'] = data
    #         dsetmeta['hist'] = hist
    #         dsetmeta['edge'] = edge
    #         if "states" in statenode.name:
    #             nicename = "ED e=[" + statenode.attrs["efmt"] + "]"
    #         else:
    #             nicename = re.sub(r'[\W_]', ' ', str(k['algo'] + " " + k['state']))
    #
    #         dsetmeta['legend'] = nicename + ' (' + str(num) + ')'

    # Let's try plotting
    # num = 0
    # for Lkey in db['keys']['L']:
    #     figrows = len(db['keys']['l'])
    #     figcols = len(db['keys']['state']) * len(db['keys']['algo'])
    #     fig, axes = plt.subplots(nrows=figrows, ncols=figcols, figsize=(7 * figrows, 7 * figcols))
    #     fig.tight_layout(pad=5, w_pad=1.0, h_pad=1.0)
    #     fig.subplots_adjust(wspace=0.3, hspace=0.3)
    #     used_ax = 0
    #     for lidx, lkey in enumerate(db['keys']['l']):
    #         for algoidx, algokey in enumerate(db['keys']['algo']):
    #             for stateidx, statekey in enumerate(db['keys']['state']):
    #                 if num > 10:
    #                     continue
    #                 # Let's histograms into a matrix at fixed lambda varying delta
    #                 rows = bins
    #                 cols = len(db['keys']['d'])
    #                 data2d = []
    #                 hist2d = np.zeros(shape=(rows, cols))
    #                 deltas = np.zeros(shape=(cols))
    #                 dsetkeys = [x for x in db['dsets'] if
    #                             Lkey in x and lkey in x and algokey in x and statekey in x]
    #                 if not dsetkeys:
    #                     continue
    #                 for idx, dsetkey in enumerate(dsetkeys):
    #                     dsetmeta = db['dsets'][dsetkey]
    #                     deltas[idx] = dsetmeta['d']
    #                     hist2d[:, idx] = dsetmeta['hist']
    #                     data2d.append(dsetmeta['data'])
    #                 sortIdx = np.argsort(deltas)
    #                 deltas = deltas[sortIdx]
    #                 hist2d = hist2d[:, sortIdx]
    #                 extent = [np.min(deltas), np.max(deltas), 0, db['max']]
    #
    #                 ax = axes[lidx, stateidx * len(db['keys']['algo']) + algoidx]
    #                 im = ax.imshow(hist2d,
    #                                origin='lower',
    #                                aspect='auto',
    #                                extent=extent,
    #                                interpolation='nearest',
    #                                cmap=plt.get_cmap('viridis'),  # use nicer color map
    #                                )
    #                 plt.colorbar(im, ax=ax)
    #                 # fig.colorbar(im, orientation='vertical')
    #                 ax.set_xlabel('x')
    #                 ax.set_ylabel('y')
    #                 num = num + 1
    #                 used_ax = used_ax + 1
    #     for ax in np.ravel(axes)[used_ax:]:
    #         fig.delaxes(ax)
    # h5close(h5_src)
    # plt.show()
    # exit(0)

# print('Finding unique')
# path_L = h5py_unique_finder(h5_src, filter='L_', dep=1)
# path_l = h5py_unique_finder(h5_src, filter='l_', dep=2)
# path_d = h5py_unique_finder(h5_src, filter='d_', dep=3)
# # We make one subplot for each system size L
# # Collect 1d-histograms of S_E for each delta
# print('Collecting data')
# hists = {'path_L': path_L, 'path_l': path_l, 'path_d': path_d}
# max_S = 0
# for L in path_L:
#     for l in path_l:
#         for d in path_d:
#             path = L + '/' + l + '/' + d
#             if not path in h5_src:
#                 continue
#             basenode = h5_src[path]
#             hists[path] = {}
#             hists[path]['length'] = basenode.attrs['model_size']
#             hists[path]['d'] = basenode.attrs['d']
#             hists[path]['l'] = basenode.attrs['l']
#             for algokey, algopath, algonode in h5py_group_iterator(g=basenode, filter=algo_inc, dep=1):
#                 for statekey, statepath, statenode in h5py_group_iterator(g=algonode, filter=state_inc, dep=1):
#                     for win_idx, win in enumerate(variance_window_limits):
#                         idx = get_v_inced_index_list(statenode, win)
#                         hists[path][algopath + statepath] = {}
#                         for datakey, datapath, datanode in h5py_node_finder(g=statenode, filter='entanglement_entropy',
#                                                                             dep=8):
#                             fullpath = path + algopath + statepath + datapath
#                             print('Processing path', fullpath)
#                             if not idx:
#                                 data = datanode['data']
#                             else:
#                                 data = datanode['data'][idx]
#                             if np.any(np.isnan(data)):
#                                 raise ValueError("Data contains nan's")
#                             hists[path][algopath + '/' + statepath]['name'] = datanode.name
#                             hists[path][algopath + '/' + statepath]['path'] = datapath
#                             hists[path][algopath + '/' + statepath]['node'] = datanode
#                             max_S = np.max([max_S, datanode['max'][()]])

# Now we have a list full of data


# hists['']
# num = datanode['num'][()]
# avg = datanode['avg'][()]
# datarange = [np.min(data), np.max(data)]
# hist, edges = np.histogram(data, bins=bins, range=datarange, density=False)
# bincentres = [(edges[i] + edges[i + 1]) / 2. for i in range(len(edges) - 1)]
# widths = np.diff(edges)
# norm = np.dot(hist, widths)
# if "states" in statekey:
#     color = next(ed_palette)
#     nicename  = "ED e=[" + statenode.attrs["efmt"] + "]"
#     lwidth = 2.4
#     lalpha = 0.8
# else:
#     color = next(current_palette)
#     nicename = re.sub(r'[\W_]', ' ',str(algokey + " " + statekey))
#     lwidth = 1.4
#     lalpha = 0.9
#     max_S = np.max([max_S, np.max(data)])
# nicename = nicename + ' (' + str(num) + ')'
# ax.step(bincentres, hist / norm, where='mid', label=nicename, linewidth=lwidth,alpha=lalpha,color=color)
# ax.axvline(avg, linestyle='dashed', linewidth=lwidth,alpha=lalpha,color=color)


# # One figure per unique_l, unique_J and unique_h
# for l in path_l:
#     for d in path_d:
#         # In each figure we want one subplot per unique_L
#         rows, cols = get_optimal_subplot_num(len(path_L))
#         fig,axes = plt.subplots(nrows=rows, ncols=cols, figsize=(7 * cols, 7 * rows))
#         fig.tight_layout(pad=5, w_pad=1.0, h_pad=1.0)
#         fig.subplots_adjust(wspace=0.3, hspace=0.3)
#         used_ax = 0
#         delt = 0
#         lamb = 0
#         for ax, L in zip(np.ravel(axes), path_L):
#             if not L in h5_src or not l in h5_src[L] or not d in h5_src[L][l]:
#                 continue
#             basenode = h5_src[L][l][d]
#             ed_palette = itertools.cycle(sns.color_palette("Set2"))
#             current_palette = itertools.cycle(sns.color_palette("colorblind", 5))
#             chain_length = basenode.attrs['model_size']
#             delt = basenode.attrs['d']
#             lamb = basenode.attrs['l']
#             max_S = 0
#             for algokey,algopath,algonode in h5py_group_iterator(g=basenode,filter=algo_inc,dep=1):
#                 for statekey,statepath,statenode in h5py_group_iterator(g=algonode,filter=state_inc,dep=1):
#                     for win_idx, win in enumerate(variance_window_limits):
#                         idx = get_v_inced_index_list(statenode, win)
#                         for datakey,datapath,datanode in h5py_node_finder(g=statenode,filter='entanglement_entropy',dep=8):
#                             if not idx:
#                                 data = datanode['data']
#                             else:
#                                 data = datanode['data'][idx]
#                             if np.any(np.isnan(data)):
#                                 raise ValueError("Data contains nan's")
#                             num = datanode['num'][()]
#                             avg = datanode['avg'][()]
#                             datarange = [np.min(data), np.max(data)]
#                             hist, edges = np.histogram(data, bins=bins, range=datarange, density=False)
#                             bincentres = [(edges[i] + edges[i + 1]) / 2. for i in range(len(edges) - 1)]
#                             widths = np.diff(edges)
#                             norm = np.dot(hist, widths)
#                             if "states" in statekey:
#                                 color = next(ed_palette)
#                                 nicename  = "ED e=[" + statenode.attrs["efmt"] + "]"
#                                 lwidth = 2.4
#                                 lalpha = 0.8
#                             else:
#                                 color = next(current_palette)
#                                 nicename = re.sub(r'[\W_]', ' ',str(algokey + " " + statekey))
#                                 lwidth = 1.4
#                                 lalpha = 0.9
#                                 max_S = np.max([max_S, np.max(data)])
#                             nicename = nicename + ' (' + str(num) + ')'
#                             ax.step(bincentres, hist / norm, where='mid', label=nicename, linewidth=lwidth,alpha=lalpha,color=color)
#                             ax.axvline(avg, linestyle='dashed', linewidth=lwidth,alpha=lalpha,color=color)
#
#                 if max_S > 0:
#                     ax.set_xlim(0,max_S)
#                 ax.set_xlabel('$S_E$')
#                 ax.set_ylabel('$P(S_E)$')
#                 ax.set_title('$L = ' + str(chain_length) + '$')
#                 # ax.set_xlim(1e-21,100)
#             used_ax = used_ax + 1
#             ax.legend()
#         fig.suptitle(
#             'Distribution of mid-chain entanglement entropy @ $\Delta = $' + str(delt) + '$\lambda = $' + str(lamb))
#         for ax in np.ravel(axes)[used_ax:]:
#             fig.delaxes(ax)
#         if plotdir != '':
#             Jh = re.sub('/', '_', str(d))
#             plt.savefig(plotdir + '/S_distribution_' + l + '_' + d + '.pdf', format='pdf')
#             plt.savefig(plotdir + '/S_distribution_' + l + '_' + d + '.png', format='png')
