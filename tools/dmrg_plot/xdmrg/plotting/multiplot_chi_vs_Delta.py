from .tools import *
from src.io.h5ops import *
import numpy as np
import matplotlib.pyplot as plt
from src.measurement.compute_statistics import *


def write_plot_chi_vs_Delta(src, tgt, sims):
    print('Writing:     Plot Chi vs Delta', sims)
    if isinstance(sims, list) and len(sims) > 1:
        dsetname = 'Chi_vs_Delta_best'
        write_plot_chi_vs_Delta_impl(src=src, tgt=tgt, sim=sims, dsetname=dsetname, legend='best')
    for sim in sims:
        suffix = sim.split("/")[-1]
        dsetname = 'Chi_vs_Delta_' + suffix
        write_plot_chi_vs_Delta_impl(src=src, tgt=tgt, sim=sim, dsetname=dsetname, legend=suffix)


def write_plot_chi_vs_Delta_impl(src, tgt, sim, dsetname, legend=''):
    h5_tgt = h5open(tgt, 'a')
    h5_src = h5open(src, 'r')

    for i, (path_l, node_l) in enumerate(h5py_group_iterator(h5_src, filter='l_')):
        num_h = len(node_l['J_0'].keys())
        SE_stats = np.zeros([7, num_h])
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
            SE_stats[6, j] = mps_chi_statistics(node_h, sim, mode='middle', compute_statistics=True)
        groupname = path_l + '/plots'
        h5_tgt.require_group(groupname)
        h5_tgt[groupname].create_dataset(dsetname, data=SE_stats)
        h5_tgt[groupname][dsetname].attrs['lambda'] = lamb
        h5_tgt[groupname][dsetname].attrs['delta'] = delta
        h5_tgt[groupname][dsetname].attrs['chain_length'] = str(chain_length)
        h5_tgt[groupname][dsetname].attrs['xlabel'] = '$\Delta = \log J_\mu - \log h_\mu $'
        h5_tgt[groupname][dsetname].attrs['ylabel'] = '$\chi$'
        h5_tgt[groupname][dsetname].attrs['title'] = 'Mid-chain bond dimension'
        h5_tgt[groupname][dsetname].attrs['std'] = 'Standard Deviation'
        h5_tgt[groupname][dsetname].attrs['ste'] = 'Standard Error'
        h5_tgt[groupname][dsetname].attrs['col0'] = '$\Delta$'
        h5_tgt[groupname][dsetname].attrs['col1'] = '$N$'
        h5_tgt[groupname][dsetname].attrs['col2'] = '$\langle \chi \\rangle$'
        h5_tgt[groupname][dsetname].attrs['col3'] = '$\sigma_{\chi}$'
        h5_tgt[groupname][dsetname].attrs['col4'] = '$\sigma_{\chi} / \sqrt{N}$'
        h5_tgt[groupname][dsetname].attrs['col5'] = '$\exp{(\langle \log \chi \\rangle)}$'
        h5_tgt[groupname][dsetname].attrs['col6'] = 'Median'
        h5_tgt[groupname][dsetname].attrs['legend'] = str(legend)

    h5close(h5_tgt)
    h5close(h5_src)


def multiplot_chi_vs_Delta(src, plotdir='', type='typical'):
    print('Plotting:     Chi vs Delta -- ' + type)
    h5_src = h5open(src, 'r')
    num_L = len(h5py_node_finder(g=h5_src, filter='L_'))
    rows, cols = get_optimal_subplot_num(num_L)
    fig, axes = plt.subplots(nrows=rows, ncols=cols, figsize=(3.5 * cols, 3.5 * rows))
    fig.tight_layout(pad=5, w_pad=1.0, h_pad=1.0)
    fig.subplots_adjust(wspace=0.3, hspace=0.3)
    used_ax = 0
    for ax, (path_L, node_L) in zip(np.ravel(axes), h5py_node_finder(g=h5_src, filter='L_')):
        for i, (path_l, node_l) in enumerate(h5py_node_iterator(g=node_L, filter='l_')):
            dsets = h5py_node_finder(node_l, filter='Time_vs_Delta')
            for elem in dsets:
                dset = elem[1]
                stats = np.array(dset).transpose()
                if type == 'typical':
                    ax.plot(np.array(stats[:, 0]), np.array(stats[:, 5]), label='$\lambda = $' + str(dset.attrs['lambda']))
                    ax.set_title(dset.attrs['title'] + ' - Typical value ' + dset.attrs['col5'])
                elif type == 'average':
                    ax.errorbar(x=stats[:, 0], y=stats[:, 2], yerr=stats[:, 4], label='$\lambda = $' + str(dset.attrs['lambda']), ls='--', capsize=2,
                                elinewidth=0.3, markeredgewidth=0.8)
                    fig.suptitle(dset.attrs['title'] + ' - Average value ' + ' vs ' + dset.attrs['col0'])
                ax.set_xlabel(dset.attrs['xlabel'])
                ax.set_ylabel(dset.attrs['ylabel'])
                ax.set_title('$L=' + dset.attrs['chain_length'] + '$')
        used_ax = used_ax + 1
        ax.legend()
    for ax in np.ravel(axes)[used_ax:]:
        fig.delaxes(ax)

    if plotdir != '':
        plt.savefig(plotdir + '/Chi_vs_Delta_' + type + '.pdf', format='pdf')
    h5close(h5_src)
