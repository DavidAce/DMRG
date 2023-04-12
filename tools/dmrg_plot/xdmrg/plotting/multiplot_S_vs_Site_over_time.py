from .tools import *
from src.io.h5ops import *
import matplotlib.pyplot as plt
from .filter import *


# import random

def multiplot_S_vs_Site_over_time(src, plotdir='', algo_filter='', state_filter='', time=[0, 10], time_num=0):
    print('Plotting: S vs Site up to time', time, 'for : ', algo_filter, state_filter)
    h5_src = h5open(src, 'r')
    path_L = h5py_unique_finder(h5_src, filter='L_', dep=1)
    path_l = h5py_unique_finder(h5_src, filter='l_', dep=2)
    path_J = h5py_unique_finder(h5_src, filter='J_', dep=3)
    path_h = h5py_unique_finder(h5_src, filter='h_', dep=4)
    path_d = []
    for J in path_J:
        for h in path_h:
            path_d.append(J + '/' + h)

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
                for ax, L in zip(np.ravel(axes), path_L):
                    ed_palette = itertools.cycle(sns.color_palette("Set2"))
                    current_palette = itertools.cycle(sns.color_palette())
                    max_S = 0
                    basenode = h5_src[L][l][d]
                    chain_length = basenode.attrs['model_size']
                    delt = basenode.attrs['delta']
                    lamb = basenode.attrs['lambda']
                    for algokey, algopath, algonode in h5py_group_iterator(g=basenode, filter=algo_filter, dep=1):
                        for statekey, statepath, statenode in h5py_group_iterator(g=algonode, filter=state_filter,
                                                                                  dep=1):
                            for datakey, datapath, datanode in h5py_node_finder(g=statenode,
                                                                                filter='entanglement_entropies',
                                                                                dep=8):
                                if "states" in statekey:
                                    ndata = datanode['num'][()]
                                    ydata = np.array(datanode['avg'])
                                    edata = np.array(datanode['ste'])
                                    xdata = range(len(ydata))
                                else:
                                    idx = get_time_filtered_index_list(statenode, time)
                                    if (len(idx[0]) == 0):
                                        continue
                                    if not idx:
                                        print(idx)
                                        data = datanode['data']
                                    else:
                                        data = datanode['data'][:, idx[0]]
                                    # data = datanode['data'][:,idx]
                                    axis = 1
                                    std = np.nanstd(data, axis=axis)
                                    num = np.shape(data)[axis]
                                    ndata = num
                                    ydata = np.nanmean(data, axis=axis)
                                    edata = std / np.sqrt(num)
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
                                    nicename = re.sub(r'[\W_]', ' ', str(algokey + " " + statekey)) + ' (random ' + str(
                                        ndata) + ')'
                                    ax.errorbar(x=xdata, y=ydata, yerr=edata, label=nicename, capsize=2,
                                                elinewidth=0.3, markeredgewidth=0.8, linewidth=0.3)
                                    max_S = np.max([max_S, np.max(ydata)])

                    ax.set_xlabel('Site $l$')
                    ax.set_ylabel('$S_E(l)$')
                    ax.set_title('$L = ' + str(chain_length) + '$')
                    used_ax = used_ax + 1
                    ax.legend()
                    if not np.isnan(max_S):
                        # ax.set_ylim(0, max_S * 1.2)
                        ax.set_ylim(0, 1.0)
                    fig.suptitle('Entanglement entropy vs site @ $\Delta = ' + str(delt) + '\quad \lambda = ' + str(
                        lamb) + '$\n'
                                'Simulation time filter: ' + str(time[0]) + ' to ' + str(time[1]) + ' minutes')
                for ax in np.ravel(axes)[used_ax:]:
                    fig.delaxes(ax)

                if plotdir != '':
                    Jh = re.sub('/', '_', str(d))
                    plt.savefig(plotdir + '/S_vs_Site_' + l + '_' + str(Jh) + '_time_' + str(time_num) + '.pdf', format='pdf')
                    # plt.savefig(plotdir + '/S_vs_Site_' + l + '_' + str(Jh) + '_time_'+ str(num) +'.png', format='png')

    h5close(h5_src)

    # # One figure per unique_L
    # # In each figure we want one subplot per unique_l and unique_h
    # for L in path_L:
    #     rows, cols = get_optimal_subplot_num(len(path_d)*len(path_l))
    #     fig, axes = plt.subplots(nrows=rows, ncols=cols, figsize=(3.5 * cols, 3.5 * rows))
    #     fig.tight_layout(pad=5, w_pad=1.0, h_pad=1.0)
    #     fig.subplots_adjust(wspace=0.3, hspace=0.3)
    #     used_ax = 0
    #     delt = 0
    #     lamb = 0
    #     length = 0
    #     for d in path_d:
    #         for l in path_l:
    #             max_S = 0
    #             ax = np.ravel(axes)[used_ax]
    #             for key in h5_src[L][l][d].keys():
    #                 node         = h5_src[L][l][d][key]['entanglement_entropies']
    #                 length       = node.attrs['chain_length']
    #                 lamb         = node.attrs['lambda']
    #                 delt         = node.attrs['delta']
    #                 ydata        = np.asarray(node['avg'])
    #                 edata        = np.asarray(node['ste'])
    #                 xdata        = range(len(ydata))
    #                 nicename = re.sub(r'[\W_]', ' ', str(key))
    #                 legend = nicename + ' (' + str(node['num'][()]) + ')'
    #                 ax.errorbar(x=xdata, y=ydata, yerr=edata, label=legend, capsize=2, elinewidth=0.3, markeredgewidth=0.8)
    #                 max_S = np.max([max_S, np.max(ydata)])
    #                 if not np.isnan(max_S):
    #                     ax.set_ylim(0, max_S * 1.2)
    #
    #             if(delt == 0 and lamb == 0 and length == 16):
    #                 S = [0, 0.374471, 0.472247, 0.523218, 0.552671, 0.573738, 0.587836, 0.594864, 0.596525, 0.594907, 0.586336, 0.572526, 0.554162, 0.523051, 0.471988, 0.373509, 0]
    #                 error = [0, 0.00102476, 0.00115219, 0.00123002, 0.00127545, 0.00130646, 0.00133261, 0.00134472, 0.00134827, 0.00133873, 0.00132296, 0.00130463, 0.00127609, 0.00123216, 0.00115692, 0.001028340, 0]
    #                 X = range(len(S))
    #                 ax.errorbar(x=X, y=S, yerr=error,label='ED', capsize=2, elinewidth=0.3,markeredgewidth=0.8)
    #             if(delt == 0 and lamb == 0 and length == 20):
    #                 S = [0,0.373, 0.471249, 0.522481, 0.55781, 0.580094, 0.596833, 0.607267, 0.613621, 0.618359, 0.618011, 0.617369, 0.612451, 0.605299, 0.593242, 0.57919, 0.557428, 0.523598, 0.475332, 0.375457,0]
    #                 X = range(len(S))
    #                 error = [0,0.00101869, 0.00114315, 0.00121947, 0.00126443, 0.00129272, 0.00132133, 0.00134265, 0.00134894, 0.00135087, 0.00135763, 0.00135267, 0.00135099, 0.00133986, 0.00132104, 0.00129401, 0.00126394, 0.0012185, 0.00114213, 0.00101948,0]
    #                 ax.errorbar(x=X, y=S, yerr=error,label='ED', capsize=2, elinewidth=0.3,markeredgewidth=0.8)
    #             if(delt == 0 and lamb == 0 and length == 24):
    #                 S = [0,0.37461736352091407, 0.4725763540937676, 0.5227031165354324, 0.5542200744364754, 0.5758221102784965, 0.5938080916214405, 0.6061505552727909, 0.615498525184133, 0.6226707105840117, 0.6260032405610142, 0.6288501829396799, 0.6299770044181504, 0.6311934512689692, 0.6291701530943633, 0.6219089008669284, 0.6143429792361674, 0.6050723264602603, 0.5946600879240719, 0.5798447124123601, 0.5566384691004214, 0.5245997914863004, 0.47141376336639035, 0.3727081401273881,0]
    #                 X = range(len(S))
    #                 error = [0,0.0010166, 0.00114359, 0.00121893, 0.0012695, 0.00129574, 0.0013202, 0.00133471, 0.00134827, 0.00135663, 0.00137044, 0.0013722, 0.00136483, 0.00136554, 0.00136365, 0.00136013, 0.00134999, 0.00133435, 0.00131441, 0.00128837, 0.00125913, 0.00121221, 0.00114432, 0.00101821,0]
    #                 ax.errorbar(x=X, y=S, yerr=error,label='ED', capsize=2, elinewidth=0.3,markeredgewidth=0.8)
    #
    #             if(delt == 0 and lamb == 0 and length == 28):
    #                 S = [0,0.373628, 0.471354, 0.524944, 0.557051, 0.578698, 0.594773, 0.608878, 0.618423,
    #                      0.625545, 0.632077, 0.637698, 0.641019, 0.642286, 0.643614, 0.642264, 0.638517,
    #                      0.637435, 0.632265, 0.626954, 0.619498, 0.607997, 0.595997, 0.579926, 0.556749,
    #                      0.524608, 0.473192, 0.373848,0]
    #                 X = range(len(S))
    #                 error = [0,0.00101788, 0.00114517, 0.00122205, 0.00126571, 0.00129543, 0.00131021,
    #                          0.00133502, 0.00135903, 0.00136457, 0.00136902, 0.00137835, 0.00137895,
    #                          0.00138376, 0.00138815, 0.00138321, 0.00137414, 0.00137163, 0.00137223,
    #                          0.00136042, 0.00134497, 0.00133506, 0.00131516, 0.00129404, 0.00126898,
    #                          0.00122117, 0.00114738, 0.00101925,0]
    #                 ax.errorbar(x=X, y=S, yerr=error,label='ED', capsize=2, elinewidth=0.3,markeredgewidth=0.8)
    #
    #             if(delt == 0 and lamb == 0 and length == 32):
    #                 S = [0,0.374089, 0.471729, 0.522618, 0.555267, 0.577773, 0.595709, 0.608997, 0.617723,
    #                      0.628575, 0.633758, 0.639552, 0.644951, 0.647703, 0.650074, 0.651359, 0.651945,
    #                      0.650796, 0.649983, 0.646755, 0.644332, 0.637778, 0.632488, 0.626712, 0.619418,
    #                      0.607465, 0.595365, 0.578531, 0.556161, 0.522405, 0.472208, 0.3755,0]
    #                 X = range(len(S))
    #                 error = [0,0.00101825, 0.00114809, 0.00121922, 0.00126688, 0.0012968, 0.00132185,
    #                          0.00133899, 0.00135762, 0.00136437, 0.00138122, 0.00138054, 0.00138744,
    #                          0.001388, 0.00139313, 0.001395, 0.00139175, 0.00140125, 0.00139372, 0.0013929,
    #                          0.00138127, 0.00137573, 0.00136636, 0.00136352, 0.00135417, 0.00134357,
    #                          0.001319, 0.00129467, 0.00125492, 0.00121383, 0.00113988, 0.00101578,0]
    #                 ax.errorbar(x=X, y=S, yerr=error,label='ED', capsize=2, elinewidth=0.3,markeredgewidth=0.8)
    #
    #             if(delt == 0 and lamb == 0 and length == 36):
    #                 S = [0,0.373552, 0.473465, 0.526165, 0.559651, 0.579821, 0.596138, 0.6084, 0.618548,
    #                      0.627826, 0.635767, 0.641627, 0.646387, 0.651642, 0.653668, 0.65635, 0.658692,
    #                      0.660192, 0.658136, 0.659071, 0.657466, 0.655574, 0.652659, 0.651128, 0.64598,
    #                      0.641066, 0.635263, 0.626342, 0.618127, 0.60799, 0.59512, 0.578946, 0.554724,
    #                      0.522394, 0.472139, 0.373299,0]
    #                 X = range(len(S))
    #                 error = [0,0.0010177, 0.00114565, 0.00121946, 0.00126327, 0.0012966, 0.00132245,
    #                          0.00134516, 0.00135709, 0.00136254, 0.00137484, 0.00138625, 0.00139435,
    #                          0.00139612, 0.00140128, 0.00140103, 0.00140282, 0.00140423, 0.00140601,
    #                          0.00140425, 0.00139926, 0.00139787, 0.00139663, 0.0013912, 0.00138992,
    #                          0.00138378, 0.00136701, 0.00135737, 0.00135367, 0.00133836, 0.00131677,
    #                          0.00129163, 0.00126626, 0.00121859, 0.00114422, 0.00101824,0]
    #                 ax.errorbar(x=X, y=S, yerr=error,label='ED', capsize=2, elinewidth=0.3,markeredgewidth=0.8)
    #
    #
    #
    #
    #             ax.set_xlabel('Site $l$')
    #             ax.set_ylabel('Entanglement entropy $S(l)$')
    #             ax.set_title('$L = $' + str(length) + '\quad $\Delta = $' + str(delt) + ' $\quad \lambda = $' + str(lamb))
    #             used_ax = used_ax + 1
    #             ax.legend()
    #     for ax in np.ravel(axes)[used_ax:]:
    #         fig.delaxes(ax)
    #     fig.suptitle('Entanglement entropy vs site')
    #
    #     if plotdir != '':
    #         plt.savefig(plotdir + '/S_vs_l_' + L + '.pdf', format='pdf')
    #

    # used_ax = 0
    # for i, (path_L, node_L) in enumerate(h5py_node_finder(g=h5_src, filter='L_')):
    #
    #     datasets = h5py_node_finder(node_L, filter='S_vs_Site')
    #     num_deltas   = []
    #     all_deltas   = []
    #     for dataset in [x[1] for x in datasets]:
    #         num_deltas.append(dataset.shape[1])
    #         all_deltas.append(dataset[0][:].T[0])
    #     max_deltas = max(num_deltas)
    #     loc        = np.argmax(num_deltas)
    #     val_deltas = all_deltas[loc]
    #
    #     rows, cols = get_optimal_subplot_num(max_deltas)
    #     fig, axes = plt.subplots(nrows=rows, ncols=cols, figsize=(3.5 * cols, 3.5 * rows))
    #     fig.tight_layout(pad=5, w_pad=1.0, h_pad=1.0)
    #     fig.subplots_adjust(wspace=0.3, hspace=0.3)
    #     fig.suptitle('Entanglement entropy vs Site -- ' + type + ' values')
    #     used_ax = 0
    #
    #     for ax, delta in zip(np.ravel(axes), val_deltas):
    #         max_S = 0
    #         for dataset in [x[1] for x in datasets]:
    #             available_deltas = dataset[0][:].T[0]
    #             idx = np.where(available_deltas == delta)[0]
    #             if len(idx) == 0:
    #                 continue
    #             else:
    #                 idx=idx[0]
    #                 if (type =='typical'):
    #                     ax.plot(dataset[5][idx], label='$\lambda = $' + str(dataset.attrs['lambda']))
    #                     max_S = np.max([max_S, np.max(dataset[5][idx])])
    #                     ax.set_ylim(0, max_S*1.2)
    #                 elif(type == 'average'):
    #                     ax.errorbar(x=range(len(dataset[2][idx])), y=dataset[2][idx], yerr=dataset[4][idx],
    #                                 label='$\lambda = $' + str(dataset.attrs['lambda']), capsize=2, elinewidth=0.3,
    #                                 markeredgewidth=0.8)
    #                     max_S = np.max([max_S, np.max(dataset[2][idx])])
    #                     ax.set_ylim(0, max_S*1.2)
    #
    #                     # ax.plot(dataset[2][idx], label='$\lambda = $' + str(dataset.attrs['lambda']))
    #
    #         ax.legend()
    #         ax.set_xlabel('Site')
    #         ax.set_ylabel(dataset.attrs['ylabel'])
    #         ax.set_title('$\Delta = ' + str(delta) + '$, $L =' + dataset.attrs['chain_length'] + '$')
    #         used_ax = used_ax + 1
    #
    #     for ax in np.ravel(axes)[used_ax:]:
    #         fig.delaxes(ax)
    #
    #     if plotdir != '':
    #         plt.savefig(plotdir + '/S_vs_Site_' + type + '_L_' + dataset.attrs['chain_length'] + '.pdf', format='pdf')
    # h5close(h5_src)


def multiplot_S_vs_l_foreach_lambda_old(src, plotdir='', type='typical'):
    print('Plotting:     S vs l at every delta, for each lambda -- ' + type)
    h5_src = h5open(src, 'r')
    used_ax = 0
    for i, (path_L, node_L) in enumerate(h5py_node_finder(g=h5_src, filter='L_')):

        datasets = h5py_node_finder(node_L, filter='S_vs_Site')
        num_deltas = []
        all_deltas = []
        for dataset in [x[1] for x in datasets]:
            num_deltas.append(dataset.shape[1])
            all_deltas.append(dataset[0][:].T[0])
        max_deltas = max(num_deltas)
        loc = np.argmax(num_deltas)
        val_deltas = all_deltas[loc]

        rows, cols = get_optimal_subplot_num(max_deltas)
        fig, axes = plt.subplots(nrows=rows, ncols=cols, figsize=(3.5 * cols, 3.5 * rows))
        fig.tight_layout(pad=5, w_pad=1.0, h_pad=1.0)
        fig.subplots_adjust(wspace=0.3, hspace=0.3)
        fig.suptitle('Entanglement entropy vs Site -- ' + type + ' values')
        used_ax = 0

        for ax, delta in zip(np.ravel(axes), val_deltas):
            max_S = 0
            for dataset in [x[1] for x in datasets]:
                available_deltas = dataset[0][:].T[0]
                idx = np.where(available_deltas == delta)[0]
                if len(idx) == 0:
                    continue
                else:
                    idx = idx[0]
                    if (type == 'typical'):
                        ax.plot(dataset[5][idx], label='$\lambda = $' + str(dataset.attrs['lambda']))
                        max_S = np.max([max_S, np.max(dataset[5][idx])])
                        ax.set_ylim(0, max_S * 1.2)
                    elif (type == 'average'):
                        ax.errorbar(x=range(len(dataset[2][idx])), y=dataset[2][idx], yerr=dataset[4][idx],
                                    label='$\lambda = $' + str(dataset.attrs['lambda']), capsize=2, elinewidth=0.3,
                                    markeredgewidth=0.8)
                        max_S = np.max([max_S, np.max(dataset[2][idx])])
                        ax.set_ylim(0, max_S * 1.2)

                        # ax.plot(dataset[2][idx], label='$\lambda = $' + str(dataset.attrs['lambda']))

            ax.legend()
            ax.set_xlabel('Site')
            ax.set_ylabel(dataset.attrs['ylabel'])
            ax.set_title('$\Delta = ' + str(delta) + '$, $L =' + dataset.attrs['chain_length'] + '$')
            used_ax = used_ax + 1

        for ax in np.ravel(axes)[used_ax:]:
            fig.delaxes(ax)

        if plotdir != '':
            plt.savefig(plotdir + '/S_vs_Site_' + type + '_L_' + dataset.attrs['chain_length'] + '.pdf', format='pdf')
    h5close(h5_src)
