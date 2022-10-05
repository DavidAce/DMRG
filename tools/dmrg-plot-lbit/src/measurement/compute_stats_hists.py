from tqdm import tqdm
from stats_settings import *
from src.general.natural_sort import *
from src.io.h5ops import *
import numpy as np


# def discard_state(hdf5_sim,path):
#     variance = hdf5_sim[path]['measurements/full/energy_variance_per_site_mpo'][0]
#     energy   = hdf5_sim[path]['measurements/full/energy_per_site_mpo'][0]
#     energy_min = hdf5_sim['xDMRG/measurements'].get('simulation_progress')['energy_min'][-1]
#     energy_max = hdf5_sim['xDMRG/measurements'].get('simulation_progress')['energy_max'][-1]
#     energy_dens = (energy - energy_min)/(energy_max - energy_min)
#     return np.abs(energy_dens - 0.5) > energy_window or variance > variance_threshold

# @profile
def best_state_selector(node, data_roots):
    # Data roots is a dict with keys 'orig' and 'proj', where 'proj' is a list
    # of paths to projections we want to consider.
    # This assumes that the paths in projections are full
    # starting from the simulation, i.e.
    #   paths = xDMRG/, xDMRG/projections/x , xDMRG/projections/+x ...

    # Check that the paths exist
    existing_data_roots = {'orig': "", 'proj': [], 'best': ""}
    existing_data_roots['orig'] = data_roots['orig'] if data_roots['orig'] in node else ""
    existing_data_roots['proj'] = [path for path in data_roots['proj'] if path in node]

    if len(existing_data_roots['proj']) <= 0:
        return existing_data_roots

    if len(existing_data_roots['proj']) == 1:
        existing_data_roots['best'] = existing_data_roots['proj'][0]
    else:
        proj_variances = [node[path + '/measurements/energy_variance_per_site'][0] for path in existing_data_roots['proj']]
        best_proj_idx = proj_variances.index(min(proj_variances))
        existing_data_roots['best'] = existing_data_roots['proj'][best_proj_idx]

    return existing_data_roots

    # existing_proj_paths = [path for path in proj_paths if path in node]
    # if not existing_proj_paths:
    #     return orig_path, False, proj_paths, True
    #
    #
    # if isinstance(existing_proj_paths,list) and len(existing_proj_paths) > 1:
    #     proj_variances          = [node[path + '/measurements/energy_variance_per_site'][0] for path in existing_proj_paths]
    #     best_proj_idx           = proj_variances.index(min(proj_variances))
    #     best_proj_path          = existing_proj_paths[best_proj_idx]
    #     return orig_path, False, best_proj_path, False
    # else:
    #     return orig_path, False, existing_proj_paths[0], False


# @profile
def get_data_multi(hdf5_sim, data_roots, data_props, processed_data):
    existing_data_roots = best_state_selector(hdf5_sim, data_roots)
    curr_root = {'orig': existing_data_roots['orig'], 'proj': existing_data_roots['best']}
    # all_data  = {'orig':{}, 'proj':{}}
    all_data = {'orig': {}, 'proj': {}}
    for prop_name, prop in data_props.items():
        for root in prop['roots']:  # Iterate through the different roots (i.e. original and best projection)
            try:
                if not curr_root[root]:  # The root does not exist in the node.. no point looking for it
                    continue
                if prop['chi_scaling']:
                    for path_chi, node_chi in h5py_group_iterator(hdf5_sim[curr_root[root]], keypattern='chi_'):
                        try:
                            chi_root = path_chi.rsplit('/', 1)[-1]
                            if not chi_root in all_data:
                                all_data[chi_root] = {}
                            if not chi_root in processed_data:
                                processed_data[chi_root] = {}
                            all_data[chi_root][prop_name] = np.atleast_1d(node_chi[prop['name']])
                        except Exception as err:
                            raise LookupError("all_data[" + chi_root + "][" + prop_name + "]: " + str(err))

                if 'dset' in prop and prop['dset']:  # Check if we are dealing with a dataset or table
                    try:
                        all_data[root][prop_name] = np.atleast_1d(hdf5_sim[curr_root[root]][prop['path']])
                        # if root == 'orig' and prop_name == 'schmidt' and np.shape(all_data[root][prop_name])[0] > 512:
                        #     print(np.shape(all_data[root][prop_name])[0])
                        #     raise ValueError("Too many schmidt values!")
                    except Exception as err:
                        print(hdf5_sim.keys())
                        print(hdf5_sim[curr_root[root]].keys())
                        raise LookupError("all_data[" + root + "][" + prop_name + "]: " + str(err))

                    # if prop['chi_scaling']:
                    #     for path_chi, node_chi in h5py_group_iterator(hdf5_sim[curr_root[root]], filter='chi_'):
                    #         chi_root = path_chi.rsplit('/', 1)[-1]
                    #         all_data[chi_root][prop_name] = np.atleast_1d(node_chi[prop['name']])

                elif 'attr' in prop and prop['attr']:  # This is an attribute
                    all_data[root][prop_name] = np.atleast_1d(hdf5_sim[curr_root[root]][prop['path']].attrs[prop['attr_name']])
                elif prop['name'] == 'renyi':
                    renyi = []
                    for path, item in h5py_node_iterator(hdf5_sim[curr_root[root]][prop['path']], keypattern='L_', dep=1):
                        if ('L_C' in path):
                            continue
                        Sq = []
                        lambdas = item[()].view(dtype=np.complex128).real
                        for q in prop['q']:
                            if q == 1:
                                Sq.append(-np.sum(lambdas ** 2 * np.log(lambdas ** 2)))
                            else:
                                Sq.append(1.0 / (1.0 - q) * np.log(np.sum(lambdas ** (2 * q))))
                        renyi.append(Sq)
                    mat = np.atleast_1d(renyi)
                    all_data[root][prop_name] = mat
                    continue
                elif prop['name'] == 'renyi2':
                    renyi = []
                    for path, item in h5py_node_iterator(hdf5_sim[curr_root[root]][prop['path']], keypattern='L_', dep=1):
                        if ('L_C' in path):
                            continue
                        Sq = []
                        lambdas = item[()].view(dtype=np.complex128).real
                        for q in prop['q']:
                            if q == 1:
                                Sq.append(-np.sum(lambdas ** 2 * np.log2(lambdas ** 2)))
                            else:
                                Sq.append(1.0 / (1.0 - q) * np.log2(np.sum(lambdas ** (2 * q))))
                        renyi.append(Sq)
                    mat = np.atleast_1d(renyi)
                    all_data[root][prop_name] = mat
                    continue

                else:  # This is a table
                    all_data[root][prop_name] = np.atleast_1d(hdf5_sim[curr_root[root]].get(prop['path'])[prop['col_name']][prop['col_idx']])

                if np.issubdtype(all_data[root][prop_name].dtype, np.integer):
                    all_data[root][prop_name] = np.array(all_data[root][prop_name], dtype=np.float64)
                # if (orig_discard and root == 'orig')  or  ( best_discard and root == 'proj'):
                #     all_data[root][prop_name][:] = np.nan
                if prop['dtype'] == np.complex128 or np.iscomplexobj(all_data[root][prop_name]):
                    all_data[root][prop_name] = all_data[root][prop['name']].view(dtype=np.complex128).real
            except Exception as er:
                print("Error in item [", hdf5_sim, "]: Could not read prop: [", prop_name, "] in root: [", root, "]. Reason: ", er)
                exit(1)
                return all_data
    return all_data


# @profile
def gather_all_realization_data(hdf5_set, data_roots, data_props, processed_data):
    key_sorted = sorted(hdf5_set.keys(), key=natural_keys)
    total_iter = len(key_sorted)

    with tqdm(total=min(total_iter, max_realizations)) as pbar:
        all_data = []
        for i, sim in enumerate(key_sorted[:max_realizations]):
            all_data.append(get_data_multi(hdf5_set[sim], data_roots, data_props, processed_data))
            pbar.update(1)
    return all_data


def sublists_same_length(inputlist):
    it = iter(inputlist)
    the_len = len(next(it))
    if not all(len(l) == the_len for l in it):
        return False
    else:
        return True
        # raise ValueError('not all lists have same length!')


def padarray(A, size, constant_values=0):
    t = size - len(A)
    return np.pad(A, pad_width=(0, t), mode='constant', constant_values=constant_values)


def get_processed_data(name, data, data_props=''):
    if name == 'schmidt':
        # print("shape: ", np.shape(data))
        # print("data:\n ", data)
        max_len = len(max(data, key=len))
        # print("schmidt max len : ", max_len)
        # print("shape before: ", np.shape(data))
        data = np.asarray([padarray(x, max_len, np.nan) for x in data])
        # print("shape after : ", np.shape(data))
    elif np.isnan(data).any():
        raise ValueError("NAN in data")
    num_elems = len(data[:, 0])
    tot_iters = len(data)
    avg = np.nanmean(data, axis=0)
    std = np.nanstd(data, axis=0)
    ste = np.nanstd(data, axis=0) / np.sqrt(num_elems)
    med = np.nanmedian(data, axis=0)
    processed_data = {'data': data, 'avg': avg, 'std': std, 'ste': ste, 'med': med, 'num': num_elems, 'tot': tot_iters}
    if name == 'renyi' and isinstance(data_props, dict) and name in data_props:
        processed_data['q'] = data_props[name]['q']
        processed_data['typ'] = np.exp(np.nanmean(np.log(data), axis=0))
    if name == 'renyi2' and isinstance(data_props, dict) and name in data_props:
        processed_data['q'] = data_props[name]['q']
        processed_data['typ'] = np.exp(np.nanmean(np.log(data), axis=0))
    if len(data[0]) > 0 and ~np.isnan(data).any():
        middle = int((len(data[0]) - 1) / 2)
        hist, edges = np.histogram(data[:, middle], bins=40, density=False)
        processed_data['hist'] = hist
        processed_data['edges'] = edges
    return processed_data


# @profile
def mps_statistics(hdf5_set, data_roots, data_props):
    processed_data = {'orig': {}, 'proj': {}}
    # In gather_all_data we modify processed_data to contain the possible roots (in case there is any  chi_scaling analysis for instance)
    all_data = gather_all_realization_data(hdf5_set, data_roots, data_props, processed_data)

    chi_roots = [str(key) for key in processed_data.keys() if "chi_" in key]
    for name in data_props:
        for root in processed_data:
            try:
                # Extract a type of data with name "name" from all_data
                if data_props[name]['chi_scaling'] and "chi_" in root:
                    # This is for chi_scaling data.
                    # The idea is that we copy chi_# into processed_data[chi_#], but if it's not found, we replace it
                    # with the final data in orig to avoid biasing the data.
                    data = []
                    for i in range(len(all_data)):
                        try:
                            if root in all_data[i]:
                                data.append(np.array(all_data[i][root][name]))
                            else:
                                # No data for this chi size present, so we use "mock" data
                                data.append(np.array(all_data[i]['orig'][name]))
                        except Exception as err:
                            raise LookupError("all_data[", str(i), "][ ", root, " ] [ ", name, " ]: ", err)
                    data = np.array(data)
                else:
                    data = np.array([np.array(all_data[i][root][name]) for i in range(len(all_data)) if root in all_data[i] and name in all_data[i][root]])
                if np.size(data) == 0:
                    continue
                processed_data[root][name] = get_processed_data(name, data, data_props)

            except KeyError as k:
                print(" KeyError when getting data from all_data: no key ", k, " in: all_data[i][ ", root, " ] [ ", name, " ]")
                pass
            except Exception as err:
                print("Error when getting data: ", err)

    return processed_data
