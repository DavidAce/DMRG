from src.plotting.tools import *
from src.io.h5ops import *
import numpy as np
import matplotlib.pyplot as plt
from src.measurement.compute_stats_hists import *
from tqdm import tqdm
import os
from pathlib import PurePath


def create_dataset(h5_node, data, data_attrs):
    for root_name, root in data.items():  # Typically orig or best
        for group_name, group in root.items():  # Typically entanglement_entropies, schmidt, etc
            for dset_name, dset in group.items():  # Typically data, avg, med, etc
                h5_node.require_group(root_name + '/' + group_name)
                h5_node[root_name][group_name].create_dataset(dset_name, data=dset)
            for name, val in data_attrs.items():
                if root_name in h5_node and group_name in h5_node[root_name]:
                    h5_node[root_name][group_name].attrs[name] = val
        for name, val in data_attrs.items():
            if root_name in h5_node:
                h5_node[root_name].attrs[name] = val


def create_dataset2(h5_node, dset_name, data, stats, num_points, data_attrs, stats_attrs):
    h5_node.create_dataset(dset_name + '_data', data=data)
    h5_node.create_dataset(dset_name + '_stats', data=stats)

    for name, val in data_attrs.items():
        h5_node[dset_name + '_data'].attrs[name] = val
        h5_node[dset_name + '_stats'].attrs[name] = val

    for name, val in stats_attrs.items():
        h5_node[dset_name + '_stats'].attrs[name] = val
    h5_node[dset_name + '_stats'].attrs['num_points'] = num_points


def compute_all_averages(src, tgt, roots, props, ed, ed_props, label=''):
    print('Writing all averages and distributions')
    h5_tgt = h5open(tgt, 'a')
    h5_src = h5open(src, 'r')  # Default ~1MB, set to 1024MB instead?
    h5_ed_tgt = h5open(ed, 'a')
    num_path_L = len(h5py_unique_finder(h5_src, keypattern='L_', dep=1))
    num_path_l = len(h5py_unique_finder(h5_src, keypattern='l_', dep=2))
    num_path_J = len(h5py_unique_finder(h5_src, keypattern='J_', dep=3))
    num_path_h = len(h5py_unique_finder(h5_src, keypattern='h_', dep=4))
    total_iter = num_path_L * num_path_l * num_path_J * num_path_h
    with tqdm(total=total_iter) as pbar:
        for L_id, (path_L, node_L) in enumerate(h5py_group_iterator(h5_src, keypattern='L_')):
            for l_id, (path_l, node_l) in enumerate(h5py_group_iterator(node_L, keypattern='l_')):
                for J_id, (path_J, node_J) in enumerate(h5py_group_iterator(node_l, keypattern='J_')):
                    for h_id, (path_h, node_h) in enumerate(h5py_group_iterator(node_J, keypattern='h_')):
                        firstnode = h5py_node_finder(node_h, keypattern='mbl', num=1, dep=1)[0][1]
                        chi_max = firstnode['xDMRG/results'].get('sim_status')['chi_max'][-1]
                        chain_length = firstnode['xDMRG/results'].get('measurements')['length'][-1]
                        J_mean = firstnode['xDMRG/model'].get('Hamiltonian')['J_mean'][-1]
                        h_mean = firstnode['xDMRG/model'].get('Hamiltonian')['h_mean'][-1]
                        lamb = firstnode['xDMRG/model'].get('Hamiltonian')['lambda'][-1]
                        delta = firstnode['xDMRG/model'].get('Hamiltonian')['delta'][-1]
                        data = mps_statistics(node_h, roots, props)
                        groupname = path_L + path_l + path_J + path_h
                        h5_tgt.require_group(groupname)
                        data_attrs = {
                            'chain_length': chain_length,
                            'delta': delta,
                            'lambda': lamb,
                            'J_mean': J_mean,
                            'h_mean': h_mean,
                            'chi_max': chi_max
                        }
                        create_dataset(h5_tgt[groupname], data, data_attrs)

                        # Try to find corresponding ED data in the ed_data folder
                        ed_data = {}
                        if not groupname in h5_ed_tgt:
                            h5_ed_tgt.require_group(groupname)
                        for dirName, subdirList, fileList in os.walk("ed_data"):
                            subdirList.sort()
                            fileList.sort()
                            for src_filename in fileList:
                                try:
                                    h5_ed_src = h5py.File(dirName + '/' + src_filename, 'r', swmr=True)
                                except:
                                    continue
                                src_filestem = PurePath(src_filename).stem
                                if src_filestem not in ed_data:
                                    ed_data[src_filestem] = {}
                                for prop_name, prop in ed_props.items():
                                    if prop_name not in ed_data:
                                        ed_data[src_filestem][prop_name] = {}
                                    for key in h5_ed_src[prop['path']].keys():
                                        if str(chain_length) == key:
                                            subgroupname = src_filestem + '/' + prop_name
                                            if not subgroupname in h5_ed_tgt[groupname]:
                                                h5_ed_tgt[groupname].require_group(subgroupname)
                                            rawdata = np.atleast_2d(h5_ed_src[prop['path']][key])
                                            ed_data[src_filestem][prop_name] = get_processed_data(prop_name, rawdata)
                        create_dataset(h5_ed_tgt[groupname], ed_data, data_attrs)
                        create_dataset(h5_tgt[groupname], ed_data, data_attrs)
                        pbar.update(1)

    h5close(h5_tgt)
    h5close(h5_src)
    h5close(h5_ed_src)
    h5close(h5_ed_tgt)

    #
    #
    #
    #     for i, (path_l, node_l) in enumerate(h5py_group_iterator(h5_src, filter='l_')):
    #         num_h = len(node_l['J_0'].keys())
    #         # num_sites = node_l['J_0/h_0/mbl_'].get('xDMRG')['chain_length'][0] + 1
    #         chain_length = int(h5py_dataset_finder(node_l['J_0'], filter='length', num=1)[0][0])
    #         # chain_length = h5py_node_finder(node_l['J_0'], filter='mbl_', num=1)[0][1][sim]['measurements']['2site']['length'][0]
    #         # print(chain_length)
    #         SE_stats = np.zeros([6,chain_length+1])
    #         for j, (path_h, node_h) in enumerate(h5py_group_iterator(node_l, filter='h_')):
    #             H0 = h5py_node_finder(node_h,filter='Hamiltonian', num=1)[0][1]
    #             J_log_mean = H0[0][3]
    #             h_log_mean = H0[0][4]
    #             lamb       = H0[0][9]
    #             delta      = (J_log_mean - h_log_mean)
    #             SE_stats[0,:] = delta
    #             SE_stats[1:, :] = mps_entanglement_entropy_statistics(node_h, sim, mode="all", compute_statistics=True)
    #
    #             groupname = path_l + '/' + path_h + '/plots'
    #             h5_tgt.require_group(groupname)
    #             dset_exists = dsetname in h5_tgt[groupname]
    #             h5_tgt[groupname].require_dataset(dsetname, data=SE_stats, shape=np.shape(SE_stats), dtype=float)
    #             pbar.update(1)
    #             if not dset_exists:
    #                 h5_tgt[groupname][dsetname].attrs['lambda'] = str(lamb)
    #                 h5_tgt[groupname][dsetname].attrs['delta']  = str(delta)
    #                 h5_tgt[groupname][dsetname].attrs['lambda_str'] = '$\lambda = $' + str(lamb)
    #                 h5_tgt[groupname][dsetname].attrs['delta_str'] = '$\Delta = $' + str(delta)
    #                 h5_tgt[groupname][dsetname].attrs['chain_length'] = str(chain_length)
    #                 h5_tgt[groupname][dsetname].attrs['xlabel'] = '$\Delta = \log J_\mu - \log h_\mu $'
    #                 h5_tgt[groupname][dsetname].attrs['ylabel'] = '$S_E$'
    #                 h5_tgt[groupname][dsetname].attrs['title'] = 'Mid-Chain Entanglement Entropy'
    #                 h5_tgt[groupname][dsetname].attrs['std']  = 'Standard Deviation'
    #                 h5_tgt[groupname][dsetname].attrs['ste']  = 'Standard Error'
    #                 h5_tgt[groupname][dsetname].attrs['col0'] = '$\Delta$'
    #                 h5_tgt[groupname][dsetname].attrs['col1'] = '$N$'
    #                 h5_tgt[groupname][dsetname].attrs['col2'] = '$\langle S_E \\rangle$'
    #                 h5_tgt[groupname][dsetname].attrs['col3'] = '$\sigma_{S_E}$'
    #                 h5_tgt[groupname][dsetname].attrs['col4'] = '$\sigma_{S_E} / \sqrt{N}$'
    #                 h5_tgt[groupname][dsetname].attrs['col5'] = 'Median'
    #                 # h5_tgt[groupname][dsetname].attrs['col6'] = '$\exp{(\langle \log S_E \\rangle)}$'
    #                 h5_tgt[groupname][dsetname].attrs['legend'] = str(legend)
    #                 h5_tgt[groupname][dsetname].attrs['var_thresh'] = variance_threshold
    # h5close(h5_tgt)
    # h5close(h5_src)
