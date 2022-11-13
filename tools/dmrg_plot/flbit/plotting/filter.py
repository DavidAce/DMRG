from src.io.h5ops import *


def get_time_maximum(node):
    if 'states' in node.name:
        return None
    result = h5py_node_finder(node=node, keypattern="algorithm_time", num=1, dep=10)
    if len(result) == 0:
        return None
    if (len(result) > 1):
        raise LookupError("Multiple algorithm time nodes found. Try being more specific:\n{}".format(result))
    time_key, time_path, time_node = result[0]
    if 'data' in time_node:
        return int(np.ceil(np.max(time_node['data']) / 60))
    else:
        raise LookupError("Node 'algorithm_time' not found in node: " + node.name)


def get_time_filtered_index_list(node, time_limits):
    if 'states' in node.name:
        return np.empty(shape=[0])
    result = h5py_node_finder(node=node, keypattern="algorithm_time", num=1, dep=10)
    if len(result) == 0:
        return np.empty(shape=[0])
    if (len(result) > 1):
        raise LookupError("Multiple algorithm time nodes found. Try being more specific:\n{}".format(result))
    time_key, time_path, time_node = result[0]
    if 'data' in time_node:
        boollist = np.logical_and(
            np.asarray(time_node['data']) >= time_limits[0] * 60,
            np.asarray(time_node['data']) <= time_limits[1] * 60
        )
        return np.where(boollist)[:]
    else:
        raise LookupError("Node 'algorithm_time' not found in node: " + node.name)


def get_v_filtered_index_list(node, filter_limits):
    if 'states' in node.name:
        return np.empty(shape=[0])
    result = h5py_node_finder(node=node, keypattern="energy_variance", num=1, dep=10)
    if len(result) == 0:
        return np.empty(shape=[0])
    if (len(result) > 1):
        raise LookupError("Multiple variance nodes found. Try being more specific:\n{}".format(result))
    var_key, var_path, var_node = result[0]
    if 'data' in var_node:
        boollist = np.logical_and(
            np.asarray(var_node['data']) >= filter_limits[0],
            np.asarray(var_node['data']) <= filter_limits[1]
        )
        return np.where(boollist)[:]
    else:
        raise LookupError("Node 'energy_variance' not found in node: " + node.name)


def get_v_filtered_edata(ydata, xdata, edata, ndata, meta, filter_limits, axis=1):
    if not filter_limits:
        return ydata, xdata, edata, ndata
    idx = get_v_filtered_index_list(meta['statenode'], filter_limits)
    if not idx:
        return ydata, xdata, edata, ndata
    idx = idx[0]
    if len(idx) == meta['datanode']['num'][()]:
        print('Did not filter')
        return ydata, xdata, edata, ndata
    data = np.array(meta['datanode']['data'])[:, idx]
    ndata = len(idx)
    ydata = np.nanmean(data, axis=axis)
    edata = np.nanstd(data, axis=axis) / np.sqrt(ndata)
    xdata = range(len(ydata))
    return ydata, xdata, edata, ndata


def get_e_filtered_index_list(node, filter_limits):
    if not 'energy_dens' in node:
        return np.empty(shape=[0])
    boollist = np.logical_and(
        np.asarray(node['energy_dens/data']) >= filter_limits[0],
        np.asarray(node['energy_dens/data']) <= filter_limits[1]
    )
    idxlist, vallist = np.where(boollist)

    return idxlist


def get_v_e_filtered_index_list(node, v_filter_limits, e_filter_limits):
    if 'states' in node.name:
        return np.empty(shape=[0])
    v_result = h5py_node_finder(node=node, keypattern="energy_variance_per_site", num=1, dep=10)
    e_result = h5py_node_finder(node=node, keypattern="energy_dens", num=1, dep=10)
    if len(v_result) == 0 and len(e_result):
        return np.empty(shape=[0])
    if (len(v_result) > 1):
        raise LookupError("Multiple variance nodes found. Try being more specific:\n{}".format(v_result))
    if (len(e_result) > 1):
        raise LookupError("Multiple energy nodes found. Try being more specific:\n{}".format(e_result))

    v_key, v_path, v_node = v_result[0]
    e_key, e_path, e_node = e_result[0]

    if 'data' in e_node and 'data' in v_node:
        e_boollist = np.logical_and(
            np.asarray(e_node['data']) >= e_filter_limits[0],
            np.asarray(e_node['data']) <= e_filter_limits[1],
        )
        v_boollist = np.logical_and(
            np.asarray(v_node['data']) >= v_filter_limits[0],
            np.asarray(v_node['data']) <= v_filter_limits[1]
        )
        v_e_boollist = np.logical_and(e_boollist, v_boollist)
        return np.where(v_e_boollist)
    elif 'data' in v_node:
        return get_v_filtered_index_list(node, v_filter_limits)
    elif 'data' in e_node:
        return get_e_filtered_index_list(node, e_filter_limits)
    else:
        raise LookupError("Nodes 'energy_dens' and 'energy_variance_per_site' not found in node: " + node.name)

