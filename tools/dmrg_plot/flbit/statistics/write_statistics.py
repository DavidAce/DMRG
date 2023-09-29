import numpy as np
import scipy.stats as sc
from plotting.tools import get_timepoints
from .write_statistics_table import *


def write_stats_to_node(data, tgt_node, axis):
    print(f'writing stat to node: {tgt_node.name}')
    std = np.nanstd(data, axis=axis)
    num = np.shape(data)[axis]
    dama = np.ma.masked_invalid(np.ma.masked_equal(np.abs(data), 0))
    # dabs = np.abs(data)
    # dama = np.manp.ma.masked_invalid(dabs)
    dtyp = np.exp(np.nanmean(np.log(np.abs(data)), axis=axis))

    if '_typ' in tgt_node.name:
        if np.any(np.isnan(dama)):
            print('dama has nans: \n{}', dama)
        tgt_node.create_dataset(name='avg', data=dtyp)
        tgt_node.create_dataset(name='typ', data=dtyp)
    else:
        tgt_node.create_dataset(name='avg', data=np.nanmean(data, axis=axis))
        tgt_node.create_dataset(name='typ', data=dtyp)

    tgt_node.create_dataset(name='std', data=std)
    tgt_node.create_dataset(name='ste', data=std / np.sqrt(num))
    tgt_node.create_dataset(name='med', data=np.nanmedian(data, axis=axis))
    tgt_node.create_dataset(name='max', data=np.nanmax(data, axis=axis))
    tgt_node.create_dataset(name='min', data=np.nanmin(data, axis=axis))
    tgt_node.create_dataset(name='q25', data=np.nanpercentile(data, 0.25, axis=axis))
    tgt_node.create_dataset(name='q75', data=np.nanpercentile(data, 0.75, axis=axis))
    tgt_node.create_dataset(name='num', data=num)
    tgt_node.create_dataset(name='data', data=data, compression="gzip", compression_opts=9 )


def get_renyi(numprobs, alpha, pn_cutoff=0.0):
    print(f'get_renyi({alpha=}, {pn_cutoff=})')
    # Assume indices are n, site, time
    data = np.copy(numprobs)
    if pn_cutoff > 0:
        data[data < pn_cutoff] = 0.0
        # Renormalize
        norm = np.linalg.norm(data, axis=0, keepdims=True)
        data = data/norm

    with np.errstate(divide='ignore'):
        r = 1.0 / (1 - alpha) * np.log(np.nansum(data ** alpha, axis=0))
    r[np.isinf(r)] = 0.0
    print(np.shape(data), '-->', np.shape(r))
    return r

def get_tdata(crononode):
    tdata = []
    for key, node in crononode.items():
        tdata.append(
            node['measurements']['physical_time'][0].astype(float)
        )
    tdata = np.sort(np.asarray(tdata))
    return tdata

def get_pinfty(numprobs, dsetnode):
    print('get_pinfty()')
    modelnode = dsetnode.parent.parent['model']
    hamiltonian = modelnode['hamiltonian']
    db = {'vals': {}}
    db['vals']['L'] = modelnode['model_size'][()]
    db['vals']['r'] = hamiltonian['J2_span'][()]
    db['vals']['x'] = modelnode['xi_Jcls'][()]
    db['vals']['w'] = (modelnode['J1_wdth'][()], modelnode['J2_wdth'][()], modelnode['J3_wdth'][()])
    tdata = get_tdata(dsetnode.parent['cronos'])
    # idx_num, idx_ent = find_saturation_idx3(tdata, db)
    t = get_timepoints(tdata,db)
    idx_num = t.idx_num_saturated
    ptavg = np.mean(numprobs[:, idx_num:, :], axis=1)
    with np.errstate(divide='ignore'):
        return -np.nansum(ptavg * np.log(ptavg), axis=0)


def get_matching_prop(props, dsetpath):
    for prop in props.keys():
        if dsetpath.endswith(prop):
            return prop
    return None

def write_statistics_dset(meta, props, h5_tgt):
    print(f'writing dset stats: {meta[0]}')
    # Props contains the names of the datasets
    dsetname = meta[0]
    dsetpath = meta[1]
    dsetnode = meta[2]
    dsetprop = get_matching_prop(props, dsetpath)
    dsetaxis = props.get(dsetprop).get('axis')
    dsetcopy = props.get(dsetprop).get('copy')
    if not dsetaxis:
        dsetaxis = 0
    if (dsetname == 'schmidt_midchain'):
        dsetdata = np.array(dsetnode.view(dtype=np.complex128).real)
    if (dsetname == 'number_probabilities'):
        Lhalf = np.shape(dsetnode)[1]//2
        print('reading number probabilities ...')
        numprobs = dsetnode[:,Lhalf,:,:]
        if props.get(dsetprop).get('hartley'):
            hartley_number_entropy_data = get_renyi(numprobs, alpha=1e-3, pn_cutoff=1e-8)
            hartley_number_entropy_path = dsetnode.parent.name + '/hartley_number_entropies'
            tgt_node = h5_tgt.require_group(hartley_number_entropy_path)
            print(f'writing dset hartley_number_entropies along axis {dsetaxis}')
            write_stats_to_node(data=hartley_number_entropy_data, tgt_node=tgt_node, axis=1)
        if props.get(dsetprop).get('renyi2'):
            print(f'writing dset renyi2_number_entropies along axis {dsetaxis}')
            renyi2_number_entropy_data = get_renyi(numprobs, alpha=2, pn_cutoff=0)
            renyi2_number_entropy_path = dsetnode.parent.name + '/renyi2_number_entropies'
            tgt_node = h5_tgt.require_group(renyi2_number_entropy_path)
            write_stats_to_node(data=renyi2_number_entropy_data, tgt_node=tgt_node, axis=1)
        if props.get(dsetprop).get('pinfty'):
            print(f'writing dset pinfty_number_entropies along axis {dsetaxis}')
            pinfty_number_entropy_data = get_pinfty(numprobs, dsetnode)
            pinfty_number_entropy_path = dsetnode.parent.name + '/pinfty_number_entropies'
            tgt_node = h5_tgt.require_group(pinfty_number_entropy_path)
            write_stats_to_node(data=pinfty_number_entropy_data, tgt_node=tgt_node, axis=0)

        if dsetcopy:
            print(f'deep copying dset: {dsetname} {np.shape(numprobs)}')
            tgt_node = h5_tgt.require_group(dsetnode.parent.name)
            zlvl = 2
            zlbl = "gzip"
            tgt_node.create_dataset(name=dsetname, data=numprobs, compression=zlbl, compression_opts=zlvl)
        else:
            tgt_node = h5_tgt.require_group(dsetpath)
            write_stats_to_node(data=numprobs, tgt_node=tgt_node, axis=dsetaxis)
    else:
        if dsetcopy:
            print(f'deep copying dset: {dsetname} {np.shape(dsetnode)}')
            tgt_node = h5_tgt.require_group(dsetnode.parent.name)
            zlvl = None if np.shape(dsetnode) == () else (1 if 'corrmat' in dsetname else 6)
            zlbl = None if np.shape(dsetnode) == () else "gzip"
            tgt_node.create_dataset(name=dsetname, data=dsetnode, compression=zlbl, compression_opts=zlvl)
        else:
            tgt_node = h5_tgt.require_group(dsetpath)
            # print('writing dset "{}" along axis {}'.format(dsetname, dsetaxis))
            write_stats_to_node(data=dsetnode[()], tgt_node=tgt_node, axis=dsetaxis)

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
        write_stats_to_node(data=data, tgt_node=energy_node, axis=0)
        data = data / L
        energy_per_site_node = measurements_node.require_group("energy_per_site")
        write_stats_to_node(data=data, tgt_node=energy_per_site_node, axis=0)

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
        write_stats_to_node(data=data, tgt_node=energy_dens_node, axis=0)

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
        write_stats_to_node(data=data, tgt_node=entropy_node, axis=1)
        width = np.shape(data)[0]
        middle = int(width / 2)
        data = data[middle, :]
        entropy_node = measurements_node.require_group("entanglement_entropy")
        write_stats_to_node(data=data, tgt_node=entropy_node, axis=0)


def write_statistics(src, tgt, reqs):
    h5_src = h5open(src, 'r')
    write_statistics_crono4.number_probabilities_path = None
    write_statistics_crono4.hartley_number_entropy_data = None


    with h5py.File(tgt, 'w') as h5_tgt:
        print('Averaging dsets')
        for dsetname, dsetpath, dsetnode in h5py_node_iterator(node=h5_src, keypattern=reqs['dsets'], dep=20, excludeKeys=['.db', 'cronos', 'iter_'],
                                                               nodeType=h5py.Dataset):
            print('Found dset: {}'.format(dsetpath))
            write_statistics_dset((dsetname, dsetpath, dsetnode), reqs['dsets'], h5_tgt)

    print('Averaging tables')
    for tablename, tablepath, tablenode in h5py_node_iterator(node=h5_src, keypattern=reqs['tables'], dep=20, excludeKeys=['.db', 'cronos', 'dsets', 'iter_'],
                                                              nodeType=h5py.Dataset):
        write_statistics_table2((tablename, tablepath, tablenode), reqs['tables'], tgt)


    with tb.File(tgt, 'a') as h5f:
        print('Averaging cronos v4')
        node_cache = {}
        done_crono = {}
        for crononame, cronopath, crononode in h5py_node_iterator(node=h5_src, keypattern='cronos', dep=20, excludeKeys=['.db', 'model', 'tables', 'dsets'],
                                                                  nodeType=h5py.Group, godeeper=False):
            print('found cronos:', cronopath)
            if done := done_crono.get(cronopath):
                continue
            else:
                for iternode in h5py_node_iterator(node=crononode, keypattern='iter_', dep=1, nodeType=h5py.Group, godeeper=False):
                    done_crono[cronopath] = write_statistics_crono4(iternode, reqs['cronos'], h5f, node_cache)
                    # print('found iter:', iternode[1])
                    if done := done_crono.get(cronopath):
                        print('{} is done'.format(cronopath))
                        break


    with h5py.File(tgt, 'a') as h5_tgt:
        # Find ED data
        print('Finding ED data')
        for dirName, subdirList, fileList in os.walk("ed_data"):
            subdirList.sort()
            fileList.sort()
            for src_filename in fileList:
                try:
                    h5_ed_src = h5py.File(dirName + '/' + src_filename, 'r', swmr=True)
                except:
                    continue

                if "ed-l_0.0000-d_[-0.5...0.5]-e[0.00-1.00].h5" in src_filename:
                    for key, path, node in h5py_node_finder(h5_tgt, keypattern="d_", dep=3, num=0):
                        # node on the target file refers to a point in the phase diagram, encoded in its path
                        # to find the same phase diagram point in ED, just compare the path!
                        if path in h5_ed_src:
                            for algo_key, algo_node in h5_ed_src[path].items():
                                node.copy(source=algo_node, dest=node)
                    break

                # We now have some ED data but we must infer what its parameters are
                # The numbers in paths such as L_16/l_0/J_0/h_0 simply enumerate
                # the values in batch simulations. The numbers are NOT the values
                # of these parameters, except in the case of L

                # Start by finding out the values in the current ED data
                # ed_lambda = h5_ed_src["lambda"][()]
                # ed_delta  = h5_ed_src["delta"][()]
                # ed_J_mean = h5_ed_src["J_mean"][()]
                # ed_h_mean = h5_ed_src["h_mean"][()]
                # ed_efmt   = h5_ed_src["efmt"][()]
                #
                # for key, path, node in h5py_node_finder(h5_tgt, filter="d_", dep=3, num=0):
                #     tgt_lambda = node.attrs["lambda"]
                #     tgt_delta  = node.attrs["delta"]
                #     tgt_J_mean = node.attrs["J_mean"]
                #     tgt_h_mean = node.attrs["h_mean"]
                #     tgt_L      = node.attrs["model_size"]
                #     if(ed_lambda != tgt_lambda):
                #         continue
                #     if(ed_delta != tgt_delta):
                #         continue
                #     if(ed_J_mean != tgt_J_mean):
                #         continue
                #     if(ed_h_mean != tgt_h_mean):
                #         continue
                #     # Now we know he ED and target data refer to the same point on the
                #     # phase diagram. Now we just have to match L data
                #     for ed_L in h5_ed_src["Ls"]:
                #         if ed_L == tgt_L:
                #             tgt_node = h5_tgt.require_group(path + "/ed-e" + ed_efmt + "/states")
                #             write_statistics_ed(h5_ed_src=h5_ed_src, tgt_node=tgt_node,L=ed_L)
                #             tgt_node.attrs["efmt"] = h5_ed_src["efmt"]
                #             tgt_node.attrs["emin"] = h5_ed_src["emin"]
                #             tgt_node.attrs["emax"] = h5_ed_src["emax"]

    h5close(h5_src)
