from .tools import *
import matplotlib.pyplot as plt
from src.measurement.compute_statistics import *


# from merge_settings import *


def process_par_distribution(src, tgt, sims):
    if isinstance(sims, list) and len(sims) > 1:
        process_par_distribution_impl(src=src, tgt=tgt, sim=sims, name='best var', legend='best var', best_by_variance=True)
        process_par_distribution_impl(src=src, tgt=tgt, sim=sims, name='best par', legend='best par', best_by_variance=False)
    for sim in sims:
        suffix = sim.split("/")[-1]
        process_par_distribution_impl(src=src, tgt=tgt, sim=sim, name=suffix, legend=suffix)


def process_par_distribution_impl(src, tgt, sim, name, legend='', best_by_variance=True):
    print('Writing:     par distribution: ', str(name))
    h5_tgt = h5open(tgt, 'a')
    h5_src = h5open(src, 'r')
    for i, (path_l, node_l) in enumerate(h5py_group_iterator(h5_src, filter='l_')):
        for j, (path_h, node_h) in enumerate(h5py_group_iterator(node_l, filter='h_')):
            groupname = path_l + path_h + '/par_distribution_' + name
            dsetname_hist = 'spin_components_hist'
            dsetname_cent = 'spin_components_cent'
            chi_max = h5py_node_finder(node_h, filter='measurements', num=1)[0][1].get('simulation_progress')['chi_max'][-1]
            H0 = h5py_node_finder(node_h, filter='Hamiltonian', num=1)[0][1]
            chain_length = h5py_dataset_finder(node_h, filter='length', num=1)[0][0]

            J_log_mean = H0[0][3]
            h_log_mean = H0[0][4]
            lamb = H0[0][9]
            delta = (J_log_mean - h_log_mean)
            data_x, data_y, data_z = np.array(mps_par_statistics(hdf5_set=node_h, data_path=sim, best_by_variance=best_by_variance))
            # edge = np.linspace(-1.01,1.01,6)
            hist_x, bins = np.histogram(data_x, bins=3, range=[-1.5, 1.5], density=False)
            hist_y, bins = np.histogram(data_y, bins=3, range=[-1.5, 1.5], density=False)
            hist_z, bins = np.histogram(data_z, bins=3, range=[-1.5, 1.5], density=False)
            cent = (bins[:-1] + bins[1:]) / 2

            h5_tgt.require_group(groupname)
            dset_exists = dsetname_hist in h5_tgt[groupname]
            h5_tgt[groupname].require_dataset(dsetname_hist, data=[hist_x, hist_y, hist_z], shape=np.shape([hist_x, hist_y, hist_z]), dtype=np.float64)
            h5_tgt[groupname].require_dataset(dsetname_cent, data=cent, shape=np.shape([cent]), dtype=np.float64)
            if not dset_exists:
                h5_tgt[groupname][dsetname_hist].attrs['lambda'] = str(lamb)
                h5_tgt[groupname][dsetname_hist].attrs['delta'] = str(delta)
                h5_tgt[groupname][dsetname_hist].attrs['chain_length'] = str(chain_length)
                h5_tgt[groupname][dsetname_hist].attrs['title'] = '$\Delta =' + str(delta) + '$'
                h5_tgt[groupname][dsetname_hist].attrs['label'] = '$\Delta =' + str(delta) + '$'
                h5_tgt[groupname][dsetname_hist].attrs['xlabel'] = 'Spin Projection'
                h5_tgt[groupname][dsetname_hist].attrs['ylabel'] = 'Histogram'
                h5_tgt[groupname][dsetname_hist].attrs['legend'] = str(legend)
                h5_tgt[groupname][dsetname_hist].attrs['chi_max'] = chi_max
                h5_tgt[groupname][dsetname_hist].attrs['var_thresh'] = variance_threshold_upper

    h5close(h5_tgt)
    h5close(h5_src)


def multiplot_par_distribution(src, plotdir=''):
    print('Plotting:     Histogram of spin components')
    h5_src = h5open(src, 'r')
    unique_L = h5py_unique_finder(g=h5_src, filter='L_')
    unique_J = h5py_unique_finder(g=h5_src, filter='J_')
    unique_h = h5py_unique_finder(g=h5_src, filter='h_')
    unique_l = h5py_unique_finder(g=h5_src, filter='l_')

    # One figure per unique_l, unique_J and unique_h
    for l in unique_l:
        for J in unique_J:
            for h in unique_h:
                # In each figure we want one subplot per unique_L
                rows, cols = get_optimal_subplot_num(len(unique_L))
                fig, axes = plt.subplots(nrows=rows, ncols=cols, figsize=(3.5 * cols, 3.5 * rows))
                fig.tight_layout(pad=5, w_pad=1.0, h_pad=1.0)
                fig.subplots_adjust(wspace=0.3, hspace=0.3)
                used_ax = 0
                for ax, (path_L, node_L) in zip(np.ravel(axes), h5py_node_finder(g=h5_src, filter='L_')):
                    node_projections = h5py_node_finder(node_L[l][J][h], filter='par_distribution_xDMRG')
                    for node in node_projections:
                        hist = h5py_node_finder(node[1], filter='spin_components_hist')[0][1]
                        cent = h5py_node_finder(node[1], filter='spin_components_cent')[0][1][0]
                        nicename = re.sub(r'[\W_]', ' ', str(hist.attrs['legend']))
                        legend = nicename
                        print(np.array(hist)[0, :])
                        print(np.array(cent))
                        print()
                        width = 0.1
                        centers = np.array([cent - width, cent, cent + width])
                        ax.bar(x=centers[0], height=np.array(hist)[0, :].T, linewidth=0.6, alpha=0.6, align='center', width=0.1, label=legend + ' sx')
                        ax.bar(x=centers[1], height=np.array(hist)[1, :].T, linewidth=0.6, alpha=0.6, align='center', width=0.1, label=legend + ' sy')
                        ax.bar(x=centers[2], height=np.array(hist)[2, :].T, linewidth=0.6, alpha=0.6, align='center', width=0.1, label=legend + ' sz')

                        ax.set_xlabel('Spin projection')
                        ax.set_ylabel('Histogram')
                        # ax.set_xscale('log')
                        ax.set_title('$L = $' + hist.attrs['chain_length'])
                    used_ax = used_ax + 1
                    ax.legend()
                delt = hist.attrs['delta']
                lamb = hist.attrs['lambda']
                fig.suptitle('Distribution of bond dimensions for each $L$ @ $\Delta = $' + str(delt) + '$\lambda = $' + str(lamb))

                for ax in np.ravel(axes)[used_ax:]:
                    fig.delaxes(ax)
                if plotdir != '':
                    plt.savefig(plotdir + '/par_distribution_' + l + '_' + J + '_' + h + '.pdf', format='pdf')

    h5close(h5_src)
