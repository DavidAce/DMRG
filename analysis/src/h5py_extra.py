
import h5py
import numpy as np

def h5py_dataset_iterator(g, prefix=''):
    for key in g.keys():
        item = g[key]
        path = '{}/{}'.format(prefix, key)
        if isinstance(item, h5py.Dataset): # test for dataset
            yield (np.asarray(item), dict(g.attrs.items()), key, path)
        elif isinstance(item, h5py.Group): # test for group (go down)
            yield from h5py_dataset_iterator(item, path)


def h5py_get_datasets(f, dset_name=''):
    temp =[]
    for (dset, attr, key, path) in h5py_dataset_iterator(f):
        if(key == dset_name):
            temp.append([dset, attr, key, path ])
    return temp

def h5py_get_datasets(f, group_name='/', dset_name=''):
    temp =[]
    if group_name in f:
        g = f[group_name]
    else:
        return []
    for (dset, attr, key, path) in h5py_dataset_iterator(g):
        if(key == dset_name):
            temp.append([dset, attr, key, path ])
    return temp
