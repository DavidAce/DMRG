import numpy as np
import h5py
from ..general.natural_sort import *


def h5open(h5filenames, permission='r', chunk_cache_mem_size=1000 * 1024 ** 2, driver=None, swmr=False):
    if len(h5filenames) == 1 or isinstance(h5filenames, str):
        try:
            return h5py.File(h5filenames, permission, swmr=swmr, rdcc_nbytes=chunk_cache_mem_size, rdcc_nslots=chunk_cache_mem_size, driver=driver)
        except Exception as err:
            print('Could not open file [{}]: {}'.format(h5filenames, err))
            exit(1)

    else:
        h5files = []
        for name in h5filenames:
            h5files.append(h5open(name, permission))
        return h5files


def h5close(h5files):
    if isinstance(h5files, h5py.File):  # Just HDF5 files
        try:
            h5files.flush()
            h5files.close()
        except:
            pass  # Was already closed

    elif isinstance(h5files, list):
        for file in h5files:
            h5close(file)
    else:
        print("Trying to close h5 files that aren't in a list.")
        exit(1)


def h5py_dataset_iterator(node, keypattern=None, dep=1):
    if dep == 0:
        return
    for key in sorted(node.keys(), key=natural_keys):
        # node = node[key]
        item = node[key]
        path = item.name
        if isinstance(item, h5py.Dataset):
            if keypattern:
                if isinstance(keypattern, str) and keypattern in path:
                    yield (key, path, item)
                elif isinstance(keypattern, list) or isinstance(keypattern, dict):
                    if any(f in path for f in keypattern):
                        yield (key, path, item)
                elif keypattern and keypattern in path:
                    yield (key, path, item)
            else:
                yield (key, path, item)
        elif isinstance(item, h5py.Group) and dep > 0:  # test for group (go down)
            yield from h5py_dataset_iterator(node=item, keypattern=keypattern, dep=dep - 1)


def h5py_group_iterator(node, keypattern=None, dep=1, excludeKeys=None):
    if dep == 0:
        return
    for key in sorted(node.keys(), key=natural_keys):
        item = node[key]
        path = item.name
        excluded = False
        if excludeKeys:
            if isinstance(excludeKeys, list):
                excluded = any(e in key for e in excludeKeys)
            elif isinstance(excludeKeys, str):
                excluded = excludeKeys in key
        if isinstance(item, h5py.Group) and not excluded:
            if keypattern:
                if isinstance(keypattern, str) and keypattern in path:
                    yield (key, path, item)
                elif isinstance(keypattern, list) or isinstance(keypattern, dict):
                    if any(f in path for f in keypattern):
                        yield (key, path, item)
                elif keypattern and keypattern in path:
                    yield (key, path, item)
            else:
                yield (key, path, item)
            if dep > 0:  # (go to subgroup)
                yield from h5py_group_iterator(node=item, keypattern=keypattern, dep=dep - 1)


def h5py_node_iterator(node, keypattern=None, dep=1, excludeKeys=None, nodeType=None, godeeper=True):
    if dep == 0:
        return
    # print(' going into {} dep {}'.format(node.name, dep))
    for key in sorted(node.keys(), key=natural_keys):
        # print(' checking', key)
        item = node[key]
        path = item.name
        excluded = False
        if excludeKeys:
            if isinstance(excludeKeys, list):
                excluded = any(e in key for e in excludeKeys)
            elif isinstance(excludeKeys, str):
                excluded = excludeKeys in key

        typeisok = isinstance(item, nodeType) if nodeType else True
        if typeisok and not excluded:
            if keypattern:
                if isinstance(keypattern, str) and keypattern in key:
                    if not godeeper:
                        dep = 0
                    # print('yield 0 dep', dep)
                    yield (key, path, item)
                elif isinstance(keypattern, list):
                    if any((f in key or path.endswith(f)) for f in keypattern):
                        if not godeeper:
                            dep = 0
                        # print('yield 1')
                        yield (key, path, item)
                elif isinstance(keypattern, dict):
                    if any((f in key or path.endswith(f)) for f in keypattern.keys()):
                        if not godeeper:
                            dep = 0
                        # print('yield 1')
                        yield (key, path, item)
                elif keypattern and keypattern in key:
                    if not godeeper:
                        dep = 0
                    # print('yield 2')
                    yield (key, path, item)
            else:
                if not godeeper:
                    dep = 0
                # print('yield 3 dep', dep)
                yield (key, path, item)
        if isinstance(item, h5py.Group) and dep > 0 and not excluded:  # test for group (go down)
            # print('yield from dep', dep)
            yield from h5py_node_iterator(node=item, keypattern=keypattern, dep=dep - 1, excludeKeys=excludeKeys, nodeType=nodeType, godeeper=godeeper)


def h5py_node_finder(node, keypattern=None, num=None, dep=1, includePath=True, excludeKeys=None, nodeType=None):
    matches = []
    for (key, path, node) in h5py_node_iterator(node=node, keypattern=keypattern, dep=dep, excludeKeys=excludeKeys, nodeType=nodeType):
        if nodeType and not isinstance(node, nodeType):
            continue
        if excludeKeys:
            if isinstance(excludeKeys, str) and excludeKeys in path:
                continue
            if isinstance(excludeKeys, list) and any(key in path for key in excludeKeys):
                continue
            if isinstance(excludeKeys, dict) and any(key in path for key in excludeKeys):
                continue
        if isinstance(num, int) and num > 0 and num <= len(matches):
            break
        matches.append((key, path, node))
    if includePath:
        return matches
    else:
        return [x[2] for x in matches]


def h5py_unique_finder(node, keypattern=None, dep=1):
    matches = h5py_node_finder(node=node, keypattern=keypattern, dep=dep)
    matches = [x[1].split("/")[-1] for x in matches]
    return list(sorted(set(matches)))


def h5py_dataset_finder(node, keypattern=None, num=None, dep=1, includePath=True, excludeKeys=None):
    return h5py_node_finder(node=node, keypattern=keypattern, num=num, dep=dep, includePath=includePath, excludeKeys=excludeKeys, nodeType=h5py.Dataset)


def h5py_group_finder(node, keypattern=None, num=None, dep=1, includePath=True, excludeKeys=None):
    return h5py_node_finder(node=node, keypattern=keypattern, num=num, dep=dep, includePath=includePath, excludeKeys=excludeKeys, nodeType=h5py.Group)


def load_component_cplx(hdf5_obj, path, keypattern=None, type=np.complex128):
    key_sorted = sorted(hdf5_obj[path].keys(), key=natural_keys)
    ret_list = []
    if not keypattern:
        for key in key_sorted:
            ret_list.append(np.asarray(hdf5_obj[path][key].value.view(dtype=np.complex128)))
    else:
        for key in filter(lambda list: keypattern in list, key_sorted):
            ret_list.append(np.asarray(hdf5_obj[path][key].value.view(dtype=np.complex128)))
    return ret_list

    # if len(ret_list) == 1:
    #     return ret_list[0]
    # else:
    #
