from src.io.h5ops import *
from tqdm import tqdm
import asyncio
import cProfile


def background(f):
    def wrapped(*args, **kwargs):
        return asyncio.get_event_loop().run_in_executor(None, f, *args, **kwargs)

    return wrapped


# @background
def mbl_collector(tgt_node, mbl_node, mbl_id, mbl_num, data_props):
    for algo_path, algo_node in [(x, mbl_node[x]) for x in data_props['algo_list'] if x in mbl_node]:
        for state_path, state_node in h5py_group_iterator(algo_node, keypattern=data_props['state_list']):
            groupname = algo_path + '/' + state_path + '/'
            tgt_node.require_group(groupname)
            dest_node = tgt_node[groupname]
            for table_path, table_node in h5py_dataset_iterator(state_node, keypattern=data_props['table_list']):
                table_name = groupname + table_path
                if (table_name in dest_node):
                    # testid = mbl_id
                    # print(np.shape(table_node[-1]))
                    dest_node[table_name][mbl_id] = table_node[-1]
                else:
                    table_shape = np.shape([table_node[-1]])
                    table_max_shape = (mbl_num * table_shape[0],)
                    table_chk_shape = (int(mbl_num / 10) * table_shape[0],)
                    table = dest_node.create_dataset(name=table_name, data=[table_node[-1]], shape=table_shape,
                                                     maxshape=table_max_shape, chunks=table_chk_shape,
                                                     dtype=table_node[-1].dtype)
                    table.resize(size=table_max_shape)
            for point_path, point_node in h5py_node_iterator(state_node, keypattern=data_props['point_list']):
                for dset_path, dset_node in h5py_dataset_iterator(point_node, keypattern=data_props['dset_list']):
                    dset_name = groupname + point_path + dset_path
                    if (dset_name in dest_node):
                        dest_node[dset_name][:, mbl_id] = dset_node
                    else:
                        dset = dest_node.create_dataset(name=dset_name, data=dset_node, shape=(dset_node.shape[0], 1),
                                                        maxshape=(dset_node.shape[0], mbl_num),
                                                        chunks=(dset_node.shape[0], int(mbl_num / 10)),
                                                        dtype=dset_node.dtype)
                        dset.resize(size=(dset_node.shape[0], mbl_num))


def collect_data(src, tgt, data_props):
    print('Collecting data')
    h5_tgt = h5open(tgt, 'a', driver='core')
    h5_src = h5open(src, 'r')  # Default ~1MB, set to 1024MB instead?
    pr = cProfile.Profile()
    for path_L, node_L in h5_src.items():
        for path_l, node_l in node_L.items():
            for path_J, node_J in node_l.items():
                for path_h, node_h in node_J.items():
                    mbl_num = len(node_h.items())
                    tgt_name = path_L + '/' + path_l + '/' + path_J + '/' + path_h + '/'
                    tgt_node = h5_tgt.require_group(tgt_name)
                    print("Collecting from", tgt_name)
                    with tqdm(total=mbl_num) as pbar:
                        # Collect the data from all realizations
                        # Enter a state in the node and take what you need
                        for mbl_id, (mbl_path, mbl_node) in enumerate(node_h.items()):
                            pr.enable()
                            mbl_collector(tgt_node, mbl_node, mbl_id, mbl_num, data_props)
                            pr.disable()
                            pbar.update(1)
                    pbar.close()

    h5close(h5_tgt)
    h5close(h5_src)
    pr.print_stats(sort="tottime")
