import os
import scipy.stats as sc
from statistics.write_statistics_table import *

import numba as nb
from numba.typed import List as TypedList
from numba.extending import overload, register_jitable

nb.set_num_threads(8)


def write_data_to_node(data, tgt_node, axis):
    std = np.nanstd(data, axis=axis)
    num = np.shape(data)[axis]

    tgt_node.create_dataset(name='avg', data=np.nanmean(data, axis=axis))
    tgt_node.create_dataset(name='std', data=std)
    tgt_node.create_dataset(name='ste', data=std / np.sqrt(num))
    tgt_node.create_dataset(name='med', data=np.nanmedian(data, axis=axis))
    tgt_node.create_dataset(name='max', data=np.nanmax(data, axis=axis))
    tgt_node.create_dataset(name='min', data=np.nanmin(data, axis=axis))
    tgt_node.create_dataset(name='q25', data=np.nanpercentile(data, 0.25, axis=axis))
    tgt_node.create_dataset(name='q75', data=np.nanpercentile(data, 0.75, axis=axis))
    tgt_node.create_dataset(name='num', data=num)
    tgt_node.create_dataset(name='data', data=data)
    dama = np.ma.masked_where(np.asarray(data) <= 0, np.asarray(data))
    tgt_node.create_dataset(name='typ',
                            data=sc.gmean(dama, axis=axis, dtype=dama.dtype))  # Geometric average


@nb.njit(nogil=False, parallel=True, cache=True)
def get_correlation_matrix_avg(corr):
    rows, cols, _ = np.shape(corr)
    corr_avg = np.empty((rows, cols))
    for i in nb.prange(rows):
        for j in nb.prange(cols):
            corr_avg[i, j] = np.mean(corr[i, j, :])
    return corr_avg


@nb.njit(nogil=False, parallel=True, cache=True)
def get_connected_correlation_matrix_avg(corr, expv):
    rows, cols, num = np.shape(corr)
    conc_data = np.empty((rows, cols, num))
    conc_avg = np.empty((rows, cols))
    for i in nb.prange(rows):
        for j in nb.prange(cols):
            conc_data[i, j, :] = corr[i, j, :] - expv[i, :] * expv[j, :]
            conc_avg[i, j] = np.mean(conc_data[i, j, :])
    return conc_avg, conc_data


def write_statistics_dset(meta, dsets, h5_tgt):
    # Props contains the names of the datasets
    dset_name, dset_path, dset_node = meta
    # dset_path = meta[1]
    if not dset_name in dsets:
        return
    data = dset_node[()]
    if "correlation" in dset_name:
        # Here we have N matrices of shape LxL.
        # We start by computing <<S_i S_j>>
        rows, cols, num = np.shape(data)
        corr_avg = get_correlation_matrix_avg(data)

        tgt_node = h5_tgt.require_group(dset_path)
        tgt_node.create_dataset(name='data', data=data)
        tgt_node.create_dataset(name='avg', data=corr_avg)
        tgt_node.create_dataset(name='num', data=num)

        # Next we compute the connected correlation
        #     <<S_i S_j>>_c =  <<S_i S_j> - <S_i><S_j>>
        # where <S_i> is the local expectation values.
        # For that we need to find the dataset "expectation_values_xx"
        dset_axis = dset_name.split("_")[-1]
        expv_name = "expectation_values_{}".format(dset_axis)
        if expv_name in dset_node.parent:
            expv = dset_node.parent[expv_name][()]
            conc_avg, conc_data = get_connected_correlation_matrix_avg(data, expv)
            conc_path = "{}/connected_matrix_{}".format(dset_node.parent.name, dset_axis)
            tgt_node = h5_tgt.require_group(conc_path)
            tgt_node.create_dataset(name='data', data=conc_data)
            tgt_node.create_dataset(name='avg', data=conc_avg)
            tgt_node.create_dataset(name='num', data=num)
        return
    elif (dset_name == "schmidt_midchain"):
        data = np.array(data.view(dtype=np.complex128).real)
    tgt_node = h5_tgt.require_group(dset_path)
    write_data_to_node(data=data, tgt_node=tgt_node, axis=1)


def write_statistics_ed(h5_ed_src, tgt_node, L):
    measurements_node = tgt_node.require_group(tgt_node.name + "/measurements")
    parity_candidates = ["/", "/0/", "/1/"]
    if "energies" in h5_ed_src:
        data = None
        for par in parity_candidates:
            path = "energies" + par + str(L)
            if path in h5_ed_src:
                if np.shape(h5_ed_src[path]) is None:
                    print("Shape", path, "is None")
                    continue
                if data is None:
                    data = np.array(h5_ed_src[path])
                else:
                    data = np.concatenate([data, np.array(h5_ed_src[path])], axis=0)
        energy_node = measurements_node.require_group("energy")
        write_data_to_node(data=data, tgt_node=energy_node, axis=0)
        data = data / L
        energy_per_site_node = measurements_node.require_group("energy_per_site")
        write_data_to_node(data=data, tgt_node=energy_per_site_node, axis=0)

    if "relativeEnergies" in h5_ed_src:
        data = None
        for par in parity_candidates:
            path = "relativeEnergies" + par + str(L)
            if path in h5_ed_src:
                if np.shape(h5_ed_src[path]) is None:
                    print("Shape", path, "is None")
                    continue
                if data is None:
                    data = np.array(h5_ed_src[path])
                else:
                    data = np.concatenate([data, np.array(h5_ed_src[path])], axis=0)
        energy_dens_node = measurements_node.require_group("energy_dens")
        write_data_to_node(data=data, tgt_node=energy_dens_node, axis=0)

    if "entropies" in h5_ed_src:
        data = None
        for par in parity_candidates:
            path = "entropies" + par + str(L)
            if path in h5_ed_src:
                if np.shape(h5_ed_src[path]) is None:
                    print("Shape", path, "is None")
                    continue
                if data is None:
                    data = np.array(h5_ed_src[path])
                else:
                    data = np.concatenate([data, np.array(h5_ed_src[path])], axis=1)

        # Pad the data
        zero_row = np.zeros([1, data.shape[1]])
        data = np.append(zero_row, data, axis=0)
        data = np.append(data, zero_row, axis=0)
        entropy_node = measurements_node.require_group("entanglement_entropies")
        write_data_to_node(data=data, tgt_node=entropy_node, axis=1)
        width = np.shape(data)[0]
        middle = int(width / 2)
        data = data[middle, :]
        entropy_node = measurements_node.require_group("entanglement_entropy_midchain")
        write_data_to_node(data=data, tgt_node=entropy_node, axis=0)


def write_statistics(src, tgt, reqs):
    h5_src = h5open(src, 'r')
    print('Averaging fes')
    for tablename, tablepath, tablenode in h5py_node_iterator(node=h5_src, keypattern=reqs['fes'], dep=20,
                                                              excludeKeys=['.db', 'cronos', 'dsets', 'tables', 'bondpoint'], nodeType=h5py.Dataset):
        write_statistics_table2((tablename, tablepath, tablenode), reqs['fes'], tgt)

    with h5open(tgt, 'a') as h5_tgt:
        print('Averaging dsets')
        for dsetname, dsetpath, dsetnode in h5py_node_iterator(node=h5_src, keypattern=reqs['dsets'], dep=20,
                                                               excludeKeys=['.db', 'cronos', 'tables', 'iter_', 'fes', 'bondpoint'], nodeType=h5py.Dataset):
            write_statistics_dset((dsetname, dsetpath, dsetnode), reqs['dsets'], h5_tgt)

    print('Averaging tables')
    for tablename, tablepath, tablenode in h5py_node_iterator(node=h5_src, keypattern=reqs['tables'], dep=20,
                                                              excludeKeys=['.db', 'cronos', 'dsets', 'iter_', 'fes', 'bondpoint'], nodeType=h5py.Dataset):
        write_statistics_table2((tablename, tablepath, tablenode), reqs['tables'], tgt)

    print('Averaging bondpoint')
    for tablename, tablepath, tablenode in h5py_node_iterator(node=h5_src, keypattern=reqs['bondpoint'], dep=20,
                                                              excludeKeys=['.db', 'cronos', 'dsets', 'tables', 'fes'], nodeType=h5py.Dataset):
        write_statistics_table2((tablename, tablepath, tablenode), reqs['bondpoint'], tgt)

    with h5open(tgt, 'a') as h5_tgt:
        # Find ED data
        for dirName, subdirList, fileList in os.walk("ed_data"):
            subdirList.sort()
            fileList.sort()
            print("Found fileList:", fileList)
            for src_filename in fileList:
                print("Checking src_filename:", src_filename)
                try:
                    h5_ed_src = h5py.File(dirName + '/' + src_filename, 'r', swmr=True)
                except:
                    continue

                if "ed-l_0.0000-d_[-0.5...0.5]-e[0.00-1.00].h5" in src_filename:
                    print("Found ED file:", src_filename)
                    for key, path, node in h5py_node_finder(h5_tgt, keypattern="d_", dep=3, num=0):
                        # node on the target file refers to a point in the phase diagram, encoded in its path
                        # to find the same phase diagram point in ED, just compare the path!
                        if path in h5_ed_src:
                            print("Found path in h5_ed_src:", path)
                            for algo_key, algo_node in h5_ed_src[path].items():
                                node.copy(source=algo_node, dest=node)
                if "ed-l_0.1000" in src_filename:
                    print("Found ED file:", src_filename)
                    for key, path, node in h5py_node_finder(h5_tgt, keypattern="d_", dep=3, num=0):
                        # node on the target file refers to a point in the phase diagram, encoded in its path
                        # to find the same phase diagram point in ED, just compare the path!
                        print("Searching for path", path, "in", h5_ed_src)
                        if path in h5_ed_src:
                            for algo_key, algo_node in h5_ed_src[path].items():
                                print("Found path in h5_ed_src:", path)
                                node.copy(source=algo_node, dest=node)

    h5close(h5_src)
