import numpy as np
import matplotlib.pyplot as plt
from src.io.h5ops import *
import scipy.stats as sc

tgt_file = "ed_data/ed-l_0.1000-d_0.0-e[0.00-1.00].h5"
parity_select = '0-'
parities = ['0', '1']
h5_tgt = h5open(tgt_file, 'w')

src_items = {
    'ED_data-0-1-Interacting-0.1.h5':
        {
            'path': "new_ed_data/ED_data-0-1-Interacting-0.1.h5",
            'lamb': 0.1,
            'efmt': "0.00-1.00",
            'emin': 0.00,
            'emax': 1.00,
        },
}
src_dsets = {
    'energies': {
        'propername': 'energy',
        'persitename': 'energy_per_site'
    },
    'entropies': {
        'propername': 'entanglement_entropies',
        'midchainname': 'entanglement_entropy_midchain'
    },
    'relativeEnergies': {
        'propername': 'energy_dens',
    },
    'normalizedEnergies': {
        'propername': 'energy_dens',
    }
}


def standard_path(Lval, lval, dval):
    return "L_{}/l_{:.4f}/d_{:+.4f}".format(int(Lval), lval, dval)


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
    tgt_node.create_dataset(name='data', data=data, compression="gzip")
    dama = np.ma.masked_where(np.asarray(data) <= 0, np.asarray(data))
    tgt_node.create_dataset(name='typ',
                            data=sc.gmean(dama, axis=axis, dtype=dama.dtype))  # Geometric average


for src_key, src_val in src_items.items():
    h5_src = h5open(src_val['path'], 'r')
    for src_dsetkey, src_dsetval in src_dsets.items():
        for dsetkey, dsetpath, dsetnode in h5py_node_iterator(node=h5_src, keypattern=src_dsetkey, dep=1):
            # Now there are two possibilities.
            # Either the data is subdivided into deltas [-0.5,-0.45...], or into parities [0,1]
            # We can tell by checking if the group name is a single digit or not
            has_delta = False
            has_parity = False
            if len(dsetnode.keys()) > 0:
                has_delta = all(not x.isdigit() for x in dsetnode.keys())
            for dkey, dpath, dnode in h5py_group_iterator(node=dsetnode, dep=1):
                # We are now inside the next level, either some delta or parity
                if dkey.isdigit():
                    dval = 0
                else:
                    dval = float(dkey)
                pnode = dnode
                # We can now check if the convention parity is at this level
                if parity_select in dnode:
                    pnode = dnode[parity_select]
                print("pnode", pnode.name)
                for Lkey, Lpath, Lnode in h5py_node_iterator(node=pnode, dep=1):
                    print("Lnode", Lnode)
                    if isinstance(Lnode, h5py.Group):
                        Lnode = Lnode["1"]
                    if not np.shape(Lnode):  # Some data is missing
                        continue
                    Lval, lval = int(Lkey), float(src_val['lamb'])

                    dlt_path = standard_path(Lval, lval, dval)
                    dlt_node = h5_tgt.require_group(dlt_path)
                    dlt_node.attrs.create(name="delta", data=dval)
                    dlt_node.attrs.create(name="lambda", data=lval)
                    dlt_node.attrs.create(name="model_size", data=Lval)
                    stt_node = h5_tgt.require_group(dlt_path + "/ed-e" + src_val["efmt"] + "/states")
                    stt_node.attrs.create(name="efmt", data=src_val["efmt"])
                    stt_node.attrs.create(name="emin", data=src_val["emin"])
                    stt_node.attrs.create(name="emax", data=src_val["emax"])
                    msm_node = h5_tgt.require_group(dlt_path + "/ed-e" + src_val["efmt"] + "/states/measurements")
                    if "propername" in src_dsetval:
                        tgt_node = h5_tgt.require_group(msm_node.name + "/" + src_dsetval["propername"])
                        if "entropies" in src_dsetval["propername"]:
                            # Pad the data
                            print("before data", Lnode)
                            data = np.asarray(Lnode)
                            zero_row = np.zeros([1, data.shape[1]])
                            data = np.append(zero_row, data, axis=0)
                            data = np.append(data, zero_row, axis=0)
                            write_data_to_node(data=data, tgt_node=tgt_node, axis=1)
                        else:
                            write_data_to_node(data=np.asarray(Lnode), tgt_node=tgt_node, axis=0)
                    if "persitename" in src_dsetval:
                        tgt_node = h5_tgt.require_group(msm_node.name + "/" + src_dsetval["persitename"])
                        write_data_to_node(data=np.asarray(Lnode) / Lval, tgt_node=tgt_node, axis=0)
                    if "midchainname" in src_dsetval:
                        tgt_node = h5_tgt.require_group(msm_node.name + "/" + src_dsetval["midchainname"])
                        middle = int(np.shape(Lnode)[0] / 2)
                        write_data_to_node(data=np.asarray(Lnode)[middle, :], tgt_node=tgt_node, axis=0)
