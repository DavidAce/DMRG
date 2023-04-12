from .tools import *
from src.io.h5ops import *
import numpy as np
import matplotlib.pyplot as plt
from src.measurement.compute_statistics import *


def write_plot_time_vs_Delta(src, tgt, sim):
    print('Writing:     Plot Time vs Delta')
    h5_tgt = h5open(tgt, 'a')
    h5_src = h5open(src, 'r')

    for i, (path_l, node_l) in enumerate(h5py_group_iterator(h5_src, filter='l_')):
        num_h = len(node_l['J_0'].keys())
        SE_stats = np.zeros([7, num_h])
        lamb = 0
        delta = 0
        chain_length = h5py_dataset_finder(node_l['J_0'], filter='length', num=1)[0][0]
        for j, (path_h, node_h) in enumerate(h5py_group_iterator(node_l, filter='h_')):
            H0 = h5py_node_finder(node_h, filter='Hamiltonian', num=1)[0][1]
            J_log_mean = H0[0][3]
            h_log_mean = H0[0][4]
            lamb = H0[0][9]
            delta = (J_log_mean - h_log_mean)
            SE_stats[0, j] = delta
            SE_stats[1, j], \
            SE_stats[2, j], \
            SE_stats[3, j], \
            SE_stats[4, j], \
            SE_stats[5, j], \
            SE_stats[6, j] = mps_simulation_time_statistics(node_h, sim, compute_statistics=True)
        groupname = path_l + '/plots'
        dsetname = 'Time_vs_Delta'
        h5_tgt.require_group(groupname)
        h5_tgt[groupname].create_dataset(dsetname, data=SE_stats)
        h5_tgt[groupname][dsetname].attrs['lambda'] = lamb
        h5_tgt[groupname][dsetname].attrs['delta'] = delta
        h5_tgt[groupname][dsetname].attrs['chain_length'] = str(chain_length)
        h5_tgt[groupname][dsetname].attrs['xlabel'] = '$\Delta = \log J_\mu - \log h_\mu $'
        h5_tgt[groupname][dsetname].attrs['ylabel'] = 'Seconds'
        h5_tgt[groupname][dsetname].attrs['title'] = 'Simulation Time [s]'
        h5_tgt[groupname][dsetname].attrs['std'] = 'Standard Deviation'
        h5_tgt[groupname][dsetname].attrs['ste'] = 'Standard Error'
        h5_tgt[groupname][dsetname].attrs['col0'] = '$\Delta$'
        h5_tgt[groupname][dsetname].attrs['col1'] = '$N$'
        h5_tgt[groupname][dsetname].attrs['col2'] = '$\langle t \\rangle$'
        h5_tgt[groupname][dsetname].attrs['col3'] = '$\sigma_{t}$'
        h5_tgt[groupname][dsetname].attrs['col4'] = '$\sigma_{t} / \sqrt{N}$'
        h5_tgt[groupname][dsetname].attrs['col5'] = '$\exp{(\langle \log t \\rangle)}$'
        h5_tgt[groupname][dsetname].attrs['col6'] = 'Median'

    h5close(h5_tgt)
    h5close(h5_src)


def multiplot_time_vs_Delta(src, plotdir='', type='typical'):
    print('Plotting:     Time vs Delta -- ' + type)
    h5_src = h5open(src, 'r')
    num_L = len(h5py_node_finder(g=h5_src, filter='L_'))
    rows, cols = get_optimal_subplot_num(num_L)
    fig, axes = plt.subplots(nrows=rows, ncols=cols, figsize=(3.5 * cols, 3.5 * rows))
    fig.tight_layout(pad=5, w_pad=1.0, h_pad=1.0)
    fig.subplots_adjust(wspace=0.3, hspace=0.3)
    used_ax = 0
    for ax, (path_L, node_L) in zip(np.ravel(axes), h5py_node_finder(g=h5_src, filter='L_')):
        for i, (path_l, node_l) in enumerate(h5py_node_iterator(g=node_L, filter='l_')):
            dset = h5py_node_finder(node_l, filter='Time_vs_Delta')[0][1]
            stats = np.array(dset).transpose()

            if type == 'typical':
                ax.plot(np.array(stats[:, 0]), np.array(stats[:, 5]), label='$\lambda = $' + str(dset.attrs['lambda']))
                ax.set_title(dset.attrs['title'] + ' - Typical value ' + dset.attrs['col5'])
            elif type == 'average':
                ax.errorbar(x=stats[:, 0], y=stats[:, 2], yerr=stats[:, 4], label='$\lambda = $' + str(dset.attrs['lambda']), capsize=2, elinewidth=0.3,
                            markeredgewidth=0.8)
                ax.set_title(dset.attrs['title'] + ' - Average value ' + dset.attrs['col2'])

            # ax.errorbar(x=stats[:,0], y=stats[:,1], yerr=stats[:,2], label=dset.attrs['col1'])
            # ax.plot(np.array(stats[:,0]), np.array(stats[:,3]), label='$\lambda = $' + str(dset.attrs['lambda']))
            # ax.plot(np.array(stats[:,0]), np.array(stats[:,4]), label=dset.attrs['col4'])
            ax.set_ylim(1e-1, 5e4)
            ax.set_yscale('log', nonpositive='clip')
            # ax.set_title (dset.attrs['title'] + ' - Typical value '+ dset.attrs['col3'])
            ax.set_xlabel(dset.attrs['xlabel'])
            ax.set_ylabel(dset.attrs['ylabel'])
        used_ax = used_ax + 1
    ax.legend()
    for ax in np.ravel(axes)[used_ax:]:
        fig.delaxes(ax)

    if plotdir != '':
        plt.savefig(plotdir + '/Time_vs_Delta_' + type + '.pdf', format='pdf')
    h5close(h5_src)
