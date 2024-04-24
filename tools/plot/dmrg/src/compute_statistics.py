from src.compute_values import *
from itertools import compress
# from src.settings import *
from tqdm import tqdm


def discard_state(hdf5_sim):
    variance = hdf5_sim['xDMRG/measurements/full/energy_variance_per_site_mpo'][0]
    energy = hdf5_sim['xDMRG/measurements/full/energy_per_site_mpo'][0]
    energy_min = hdf5_sim['xDMRG/measurements'].get('simulation_progress')['energy_min'][-1]
    energy_max = hdf5_sim['xDMRG/measurements'].get('simulation_progress')['energy_max'][-1]
    energy_dens = (energy - energy_min) / (energy_max - energy_min)
    return np.abs(energy_dens - 0.5) > energy_window or variance > variance_threshold_upper


def best_state_selector(node, paths):
    # This assumes that the paths in projections are full
    # starting from the simulation, i.e.
    #   paths = xDMRG/, xDMRG/projections/sx_up , xDMRG/projections/sx_dn ...
    # It also assumes that the first element is the "original" simulation and the rest are projections
    spin_comp_ok = []
    variances = []
    if isinstance(paths, list) and len(paths) > 1:
        for path in paths:
            spin_component = node[path + '/measurements/full/spin_component_sx'][0]
            spin_comp_ok.append(np.abs(spin_component) > 0.999)
            variances.append(node[path + '/measurements/full/energy_variance_per_site_mpo'][0])

        # The "best" state is one with low variance, and whose spin is OK, i.e. close to either + or - 1
        # At this stage there may be different situations:
        #   1) All spins are OK, and one of them is lowest. <---- most common case
        #   2) No spin component is OK.
        #   3) The state with lowest variance has a spin component not OK.

        # Filter out non-ok spin projections
        variances_filtered = list(compress(variances, spin_comp_ok))
        if len(variances_filtered) == 0:
            variances_filtered = variances

        best_var_idx = variances_filtered.index(min(variances_filtered))
        return paths[best_var_idx]


def best_parity_selector(node, original, projections):
    # This assumes that the paths in projections are full
    # starting from the simulation, i.e.
    #   original = xDMRG/
    #   projections = xDMRG/projections/sx_up , xDMRG/projections/sx_dn ...

    if isinstance(projections, list) and len(projections) > 1:

        spin_component_original = node[original[0] + '/measurements/full/spin_component_sx'][0]
        spin_component_A = node[projections[0] + '/measurements/full/spin_component_sx'][0]
        # spin_component_B        = node[projections[1] + '/measurements/full/spin_component_sx'][0]
        if spin_component_original * spin_component_A > 0:
            return projections[0]
        else:
            return projections[1]
    else:
        return original


def best_variance_selector(node, paths):
    # This assumes that the paths in projections are full
    # starting from the simulation, i.e.
    #   xDMRG/projections/sx_up
    #   xDMRG/projections/sx_dn
    # etc.

    if isinstance(paths, list) and len(paths) > 1:
        variances = []
        for path in paths:
            variances.append(node[path + '/measurements/full/energy_variance_per_site_mpo'][0])
        idx = variances.index(min(variances))
        return paths[idx]
    else:
        return paths


def mps_entanglement_entropy_statistics(hdf5_set, data_path, mode='all', compute_statistics=True):
    chain_entropy = []
    key_sorted = sorted(hdf5_set.keys(), key=natural_keys)
    # projection_paths = [x for x in data_path if any(y in x for y in  ['sx_up','sx_dn'])]
    # original_path    = [x for x in data_path if not any(y in x for y in  ['sx_up','sx_dn'])]
    total_iter = len(key_sorted)
    with tqdm(total=min(max_realizations, total_iter)) as pbar:
        if isinstance(data_path, list) and len(data_path) > 1:
            for i, sim in enumerate(key_sorted):
                pbar.update(1)
                if i > max_realizations:
                    break
                if discard_state(hdf5_set[sim]):
                    chain_entropy.append(np.nan * mps_entanglement_entropy(hdf5_set[sim][data_path[0]], mode=mode))
                    continue
                best_path = best_state_selector(hdf5_set[sim], data_path)
                chain_entropy.append(mps_entanglement_entropy(hdf5_set[sim][best_path], mode=mode))
        else:
            for i, sim in enumerate(key_sorted):
                pbar.update(1)
                if i > max_realizations:
                    break
                if discard_state(hdf5_set[sim]):
                    chain_entropy.append(np.nan * mps_entanglement_entropy(hdf5_set[sim][data_path], mode=mode))
                    continue

                chain_entropy.append(mps_entanglement_entropy(hdf5_set[sim][data_path], mode=mode))

    if all(isinstance(elem, np.ndarray) for elem in chain_entropy):
        results = []
        length = len(chain_entropy[0])
        chain_matrix = np.array(chain_entropy)
        non_nan_iters = np.count_nonzero(~np.isnan(chain_matrix[:, 0]))
        for i in range(length):
            results.append([length,
                            np.nanmean(chain_matrix[:, i]),
                            np.nanstd(chain_matrix[:, i]),
                            np.nanstd(chain_matrix[:, i]) / np.sqrt(non_nan_iters),
                            np.nanmedian(chain_matrix[:, i])])
            # np.exp(np.nanmean(np.log(chain_matrix[:,i]))),
        return np.array(results).T

    else:
        # Filter out weird zeros if they are present, otherwise they ruin the typical value calculation
        if np.any(np.array(chain_entropy) > 1e-10):
            chain_entropy = list(filter(lambda a: a > 1e-10, chain_entropy))
        # chain_entropy = chain_entropy[np.where(np.array(chain_entropy) > 1e-8)[0]]
        non_nan_iters = np.count_nonzero(~np.isnan(chain_entropy))
        if compute_statistics:
            return len(chain_entropy), \
                   np.nanmean(chain_entropy, axis=0), \
                   np.nanstd(chain_entropy, axis=0), \
                   np.nanstd(chain_entropy, axis=0) / np.sqrt(non_nan_iters), \
                   np.nanmedian(chain_entropy, axis=0)
            # np.exp(np.nanmean(np.log(chain_entropy),axis=0)), \
        else:
            return np.array(chain_entropy).T.tolist()


def mps_simulation_time_statistics(hdf5_set, data_path, compute_statistics=True):
    simulation_time = []
    key_sorted = sorted(hdf5_set.keys(), key=natural_keys)
    for sim in key_sorted:
        simulation_time.append(hdf5_set[sim][data_path]['measurements'].get('simulation_progress')['wall_time'][-1])
    if compute_statistics:
        return len(simulation_time), \
               np.mean(simulation_time, axis=0), \
               np.std(simulation_time, axis=0), \
               np.std(simulation_time, axis=0) / np.sqrt(len(simulation_time)), \
               np.exp(np.mean(np.log(simulation_time), axis=0)), \
               np.median(simulation_time)
    else:
        return np.array(simulation_time).T.tolist()


def mps_chi_statistics(hdf5_set, data_path, mode='all', compute_statistics=True):
    chi = []
    key_sorted = sorted(hdf5_set.keys(), key=natural_keys)
    for sim in key_sorted:
        if isinstance(data_path, list) and len(data_path) > 1:
            best_path = best_variance_selector(hdf5_set[sim], data_path)
            chi.append(mps_bond_dimension(hdf5_set[sim][best_path], mode=mode))
        else:
            chi.append(mps_bond_dimension(hdf5_set[sim][data_path], mode=mode))
    if compute_statistics:
        return len(chi), \
               np.mean(chi, axis=0), \
               np.std(chi, axis=0), \
               np.std(chi, axis=0) / np.sqrt(len(chi)), \
               np.exp(np.mean(np.log(chi), axis=0)), \
               np.median(chi, axis=0)
    else:
        return np.array(chi).T.tolist()


def mps_par_statistics(hdf5_set, data_path, best_by_variance=True):
    par_sx = []
    par_sy = []
    par_sz = []
    projection_paths = [x for x in data_path if any(y in x for y in ['sx_up', 'sx_dn'])]
    key_sorted = sorted(hdf5_set.keys(), key=natural_keys)
    if isinstance(data_path, list) and len(data_path) > 1:
        if best_by_variance:
            for sim in key_sorted:
                best_path = best_variance_selector(hdf5_set[sim], projection_paths)
                par_sx.append(np.asarray(hdf5_set[sim][best_path + '/measurements/full/spin_component_sx']))
                par_sy.append(np.asarray(hdf5_set[sim][best_path + '/measurements/full/spin_component_sy']))
                par_sz.append(np.asarray(hdf5_set[sim][best_path + '/measurements/full/spin_component_sz']))
        else:
            for sim in key_sorted:
                best_path = best_parity_selector(hdf5_set[sim], data_path, projection_paths)
                par_sx.append(np.asarray(hdf5_set[sim][best_path + '/measurements/full/spin_component_sx']))
                par_sy.append(np.asarray(hdf5_set[sim][best_path + '/measurements/full/spin_component_sy']))
                par_sz.append(np.asarray(hdf5_set[sim][best_path + '/measurements/full/spin_component_sz']))
    else:
        for sim in key_sorted:
            par_sx.append(np.asarray(hdf5_set[sim][data_path + '/measurements/full/spin_component_sx']))
            par_sy.append(np.asarray(hdf5_set[sim][data_path + '/measurements/full/spin_component_sy']))
            par_sz.append(np.asarray(hdf5_set[sim][data_path + '/measurements/full/spin_component_sz']))
    return par_sx, par_sy, par_sz


# def mps_entanglement_entropy_statistics(hdf5_set,data_path, mode='all', compute_statistics=True, compute_typical=True,best_by_variance=True):
#     chain_entropy = []
#     key_sorted = sorted(hdf5_set.keys(), key=natural_keys)
#     projection_paths = [x for x in data_path if any(y in x for y in  ['sx_up','sx_dn'])]
#     # original_path    = [x for x in data_path if not any(y in x for y in  ['sx_up','sx_dn'])]
#     if isinstance(data_path, list) and len(data_path) > 1:
#         if (best_by_variance):
#             for sim in key_sorted:
#                 best_path = best_variance_selector(hdf5_set[sim], projection_paths)
#                 chain_entropy.append(mps_entanglement_entropy(hdf5_set[sim][best_path], mode=mode))
#         else:
#             for sim in key_sorted:
#                 spin_component = hdf5_set[sim][data_path[0] + '/measurements/full/spin_component_sx'][0]
#                 if spin_component > 0:
#                     chain_entropy.append(mps_entanglement_entropy(hdf5_set[sim][data_path[1]], mode=mode))
#                 else:
#                     chain_entropy.append(mps_entanglement_entropy(hdf5_set[sim][data_path[2]], mode=mode))


def mps_energy_statistics(hdf5_set, data_path, compute_statistics=True):
    eden = []
    key_sorted = sorted(hdf5_set.keys(), key=natural_keys)
    for sim in key_sorted:
        emin = hdf5_set[sim][data_path]['sim_state']['energy_min'][0]
        emax = hdf5_set[sim][data_path]['sim_state']['energy_max'][0]
        efin = hdf5_set[sim][data_path]['sim_state']['energy_now'][0]
        eden.append((efin - emin) / (emax - emin))
    if compute_statistics:
        return len(eden), \
               np.mean(eden, axis=0), \
               np.std(eden, axis=0), \
               np.std(eden, axis=0) / np.sqrt(len(eden)), \
               np.exp(np.mean(np.log(eden), axis=0)), \
               np.median(eden)
    else:
        return np.array(eden).T.tolist()


def mps_variance_statistics(hdf5_set, data_path, compute_statistics=True):
    var = []
    key_sorted = sorted(hdf5_set.keys(), key=natural_keys)
    for sim in key_sorted:
        if isinstance(data_path, list) and len(data_path) > 1:
            best_path = best_variance_selector(hdf5_set[sim], data_path)
            var.append(hdf5_set[sim][best_path]['measurements']['2site']['energy_variance_per_site_mpo'][0])
        else:
            var.append(hdf5_set[sim][data_path]['measurements']['2site']['energy_variance_per_site_mpo'][0])
    if compute_statistics:
        return len(var), \
               np.mean(var, axis=0), \
               np.std(var, axis=0), \
               np.std(var, axis=0) / np.sqrt(len(var)), \
               np.exp(np.mean(np.log(var), axis=0)), \
               np.median(var)
    else:
        return np.array(var).T.tolist()
