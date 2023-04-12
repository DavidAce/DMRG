from .tools import *
from src.io.h5ops import *
import numpy as np
import matplotlib.pyplot as plt
from src.measurement.compute_statistics import *


def write_scatter_e_vs_delta(src, tgt, sim):
    print('Writing:     Scatter e vs Delta')
    h5_tgt = h5open(tgt, 'a')
    h5_src = h5open(src, 'r')

    for i, (path_l, node_l) in enumerate(h5py_group_iterator(h5_src, filter='l_')):
        num_h = len(node_l['J_0'].keys())
        chain_length = h5py_dataset_finder(node_l['J_0'], filter='length', num=1)[0][0]
        for j, (path_h, node_h) in enumerate(h5py_group_iterator(node_l, filter='h_')):
            H0 = h5py_node_finder(node_h, filter='Hamiltonian', num=1)[0][1]
            J_log_mean = H0[0][3]
            h_log_mean = H0[0][4]
            lamb = H0[0][9]
            delta = (J_log_mean - h_log_mean)
            E_data = mps_energy_statistics(node_h, sim, compute_statistics=False)

            groupname = path_l + '/delta_' + str(j)
            dsetname = 'energy_density'
            h5_tgt.require_group(groupname)
            h5_tgt[groupname].create_dataset(dsetname, data=E_data)
            h5_tgt[groupname][dsetname].attrs['lambda'] = str(lamb)
            h5_tgt[groupname][dsetname].attrs['delta'] = str(delta)
            h5_tgt[groupname][dsetname].attrs['chain_length'] = str(chain_length)
            h5_tgt[groupname][dsetname].attrs['ylabel'] = '$\epsilon$'
            h5_tgt[groupname][dsetname].attrs['title'] = 'Energy density'
            h5_tgt[groupname][dsetname].attrs['col0'] = 'Energy density'
    h5close(h5_tgt)
    h5close(h5_src)


def multiplot_e_vs_delta_and_lambda_scatter(src, plotdir=''):
    print('Plotting:     S vs l at every delta, for each lambda')
    h5_src = h5open(src, 'r')
    for i, (path_L, node_L) in enumerate(h5py_node_finder(g=h5_src, filter='L_')):
        lambdanodes = h5py_node_finder(node_L, filter='l_')
        num_lambdas = len(lambdanodes)
        rows, cols = get_optimal_subplot_num(num_lambdas)
        fig, axes = plt.subplots(nrows=rows, ncols=cols, figsize=(6 * cols, 6 * rows))
        fig.tight_layout(pad=5, w_pad=1.0, h_pad=1.0)
        fig.subplots_adjust(wspace=0.3, hspace=0.3)
        used_ax = 0
        chain_length = 0
        for ax, (path_l, node_l) in zip(np.ravel(axes), h5py_node_finder(g=node_L, filter='l_')):
            for (path_d, node_d) in h5py_node_finder(node_l, filter='delta_'):
                delta = node_d['energy_density'].attrs['delta']
                lamb = node_d['energy_density'].attrs['lambda']
                chain_length = node_d['energy_density'].attrs['chain_length']
                yaxis = node_d['energy_density']

                xaxis = np.full(len(yaxis), delta)
                ax.scatter(xaxis, yaxis, s=0.25)
                ax.set_ylim(0, 1)
                ax.set_xlabel('$\Delta$')
                ax.set_ylabel('$\epsilon$')
                ax.set_title('$\lambda =$' + str(lamb))
            used_ax = used_ax + 1
        fig.suptitle('Energy density vs Delta, $L=' + str(chain_length) + '$')
        for ax in np.ravel(axes)[used_ax:]:
            fig.delaxes(ax)

        if plotdir != '':
            plt.savefig(plotdir + '/E_vs_Delta_scatter' + '_L_' + str(chain_length) + '.pdf', format='pdf')
    h5close(h5_src)


def multiplot_e_vs_delta_and_lambda_probdist(src, plotdir=''):
    print('Plotting:     S vs l at every delta, for each lambda')
    h5_src = h5open(src, 'r')
    for i, (path_L, node_L) in enumerate(h5py_node_finder(g=h5_src, filter='L_')):
        lambdanodes = h5py_node_finder(node_L, filter='l_')
        num_lambdas = len(lambdanodes)
        rows, cols = get_optimal_subplot_num(num_lambdas)
        fig, axes = plt.subplots(nrows=rows, ncols=cols, figsize=(6 * cols, 6 * rows))
        fig.tight_layout(pad=5, w_pad=1.0, h_pad=1.0)
        fig.subplots_adjust(wspace=0.3, hspace=0.3)

        used_ax = 0
        chain_length = 0
        for ax, (path_l, node_l) in zip(np.ravel(axes), h5py_node_finder(g=node_L, filter='l_')):
            for (path_d, node_d) in h5py_node_finder(node_l, filter='delta_'):
                delta = node_d['energy_density'].attrs['delta']
                lamb = node_d['energy_density'].attrs['lambda']
                chain_length = node_d['energy_density'].attrs['chain_length']
                yaxis = node_d['energy_density']
                xaxis = np.full(len(yaxis), delta)

                # density = stats.gaussian_kde(yaxis)
                n, x, _ = ax.hist(yaxis, bins=np.linspace(0, 1, 40),
                                  histtype=u'step', density=False)
                # ax.plot(x, density(x))

                # ax.scatter(xaxis,yaxis,s=0.25)
                # ax.set_ylim(0,1)
                ax.set_xlabel('$\epsilon$')
                ax.set_ylabel('$P(\epsilon)$')
                ax.set_title('$\lambda =$' + str(lamb))
            used_ax = used_ax + 1
        fig.suptitle('Distribution of energy densities, $L=' + str(chain_length) + '$')
        for ax in np.ravel(axes)[used_ax:]:
            fig.delaxes(ax)

        if plotdir != '':
            plt.savefig(plotdir + '/E_vs_Delta_probdist' + '_L_' + str(chain_length) + '.pdf', format='pdf')
    h5close(h5_src)
