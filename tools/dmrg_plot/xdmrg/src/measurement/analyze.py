from src.plotting.multiplot import *


def write_boxplot_S_vs_site_foreach_Delta(src, tgt, sim):
    print('Writing:     Boxplot S vs site')
    h5_tgt = h5open(tgt, 'a')
    h5_src = h5open(src, 'r')
    for i, (path_l, node_l) in enumerate(h5py_group_iterator(h5_src, keypattern='l_')):
        for j, (path_h, node_h) in enumerate(h5py_group_iterator(node_l, keypattern='h_')):
            SE_data = np.array(mps_entanglement_entropy_statistics(node_h, sim, mode='all', compute_statistics=False))
            groupname = path_l + path_h
            dsetname = 'S_vs_L'
            H0 = h5py_node_finder(node_h, keypattern='H_0', num=1)
            J_log_mean = H0[0][1].attrs['J_log_mean']
            h_log_mean = H0[0][1].attrs['h_log_mean']
            lamb = H0[0][1].attrs['lambda']
            delta = J_log_mean - h_log_mean

            h5_tgt.require_group(groupname)
            h5_tgt[groupname].create_dataset(dsetname, data=SE_data)
            h5_tgt[groupname][dsetname].attrs['lambda'] = str(lamb)
            h5_tgt[groupname][dsetname].attrs['delta'] = str(delta)
            h5_tgt[groupname][dsetname].attrs['title'] = '$\Delta =' + str(delta) + '$'
            h5_tgt[groupname][dsetname].attrs['label'] = '$\Delta =' + str(delta) + '$'
            h5_tgt[groupname][dsetname].attrs['xlabel'] = 'Postion'
            h5_tgt[groupname][dsetname].attrs['xaxis'] = range(0, SE_data.shape[1], 1)
            h5_tgt[groupname][dsetname].attrs['ylabel'] = 'Entanglement Entropy'
    h5close(h5_tgt)
    h5close(h5_src)


def write_histogram_S_foreach_Delta(src, tgt, sim):
    print('Writing:     Histogram S for each Delta')
    h5_tgt = h5open(tgt, 'a')
    h5_src = h5open(src, 'r')
    for i, (path_l, node_l) in enumerate(h5py_group_iterator(h5_src, keypattern='l_')):
        chain_length = h5py_dataset_finder(node_l['J_0'], keypattern='length', num=1)[0][0]
        for j, (path_h, node_h) in enumerate(h5py_group_iterator(node_l, keypattern='h_')):
            groupname = path_l + path_h + '/S_histogram'
            dsetname1 = 'S_hist'
            dsetname2 = 'S_edges'
            H0 = h5py_node_finder(node_h, keypattern='Hamiltonian', num=1)[0][1]
            J_log_mean = H0[0][3]
            h_log_mean = H0[0][4]
            lamb = H0[0][9]
            delta = (J_log_mean - h_log_mean)
            SE_data = np.array(mps_entanglement_entropy_statistics(node_h, sim, mode='middle', compute_statistics=False))
            SE_hist, SE_edges = np.histogram(SE_data, bins=100, density=False)

            try:
                del h5_tgt[groupname][dsetname1]
            except:
                pass
            try:
                del h5_tgt[groupname][dsetname2]
            except:
                pass

            h5_tgt.require_group(groupname)
            h5_tgt[groupname].create_dataset(dsetname1, data=SE_data)
            h5_tgt[groupname].create_dataset(dsetname2, data=SE_edges)
            h5_tgt[groupname][dsetname1].attrs['lambda'] = str(lamb)
            h5_tgt[groupname][dsetname1].attrs['delta'] = str(delta)
            h5_tgt[groupname][dsetname1].attrs['chain_length'] = str(chain_length)
            h5_tgt[groupname][dsetname1].attrs['title'] = '$\Delta =' + str(delta) + '$'
            h5_tgt[groupname][dsetname1].attrs['label'] = '$\Delta =' + str(delta) + '$'
            h5_tgt[groupname][dsetname1].attrs['xlabel'] = 'Entanglement Entropy'
            h5_tgt[groupname][dsetname1].attrs['bin_edges'] = SE_edges
            h5_tgt[groupname][dsetname1].attrs['ylabel'] = 'Histogram'
    h5close(h5_tgt)
    h5close(h5_src)


def write_histogram_var_foreach_Delta(src, tgt, sim):
    print('Writing:     Histogram Var for each Delta')
    h5_tgt = h5open(tgt, 'a')
    h5_src = h5open(src, 'r')
    for i, (path_l, node_l) in enumerate(h5py_group_iterator(h5_src, keypattern='l_')):
        chain_length = h5py_dataset_finder(node_l['J_0'], keypattern='length', num=1)[0][0]
        for j, (path_h, node_h) in enumerate(h5py_group_iterator(node_l, keypattern='h_')):
            groupname = path_l + path_h + '/var_histogram'
            dsetname1 = 'var_hist'
            dsetname2 = 'var_edges'

            H0 = h5py_node_finder(node_h, keypattern='Hamiltonian', num=1)[0][1]
            J_log_mean = H0[0][3]
            h_log_mean = H0[0][4]
            lamb = H0[0][9]
            delta = (J_log_mean - h_log_mean)
            var_data = np.array(mps_variance_statistics(node_h, sim, compute_statistics=False))
            var_hist, var_edges = np.histogram(var_data, bins=100, density=False)
            try:
                del h5_tgt[groupname][dsetname1]
            except:
                pass
            try:
                del h5_tgt[groupname][dsetname2]
            except:
                pass

            h5_tgt.require_group(groupname)
            h5_tgt[groupname].create_dataset(dsetname1, data=var_data)
            h5_tgt[groupname].create_dataset(dsetname2, data=var_edges)
            h5_tgt[groupname][dsetname1].attrs['lambda'] = lamb
            h5_tgt[groupname][dsetname1].attrs['delta'] = delta
            h5_tgt[groupname][dsetname1].attrs['chain_length'] = str(chain_length)
            h5_tgt[groupname][dsetname1].attrs['title'] = '$\Delta =' + str(delta) + '$'
            h5_tgt[groupname][dsetname1].attrs['label'] = '$\Delta =' + str(delta) + '$'
            h5_tgt[groupname][dsetname1].attrs['xlabel'] = 'Variance $\sigma^2 $'
            h5_tgt[groupname][dsetname1].attrs['bin_edges'] = var_edges
            h5_tgt[groupname][dsetname1].attrs['ylabel'] = 'Variance Histogram'
    h5close(h5_tgt)
    h5close(h5_src)
