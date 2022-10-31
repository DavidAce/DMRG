from src.io.h5ops import *
import tables as tb
import tables.parameters
import warnings
from tables import NaturalNameWarning
from timeit import default_timer as timer
import numba as nb

warnings.filterwarnings('ignore', category=NaturalNameWarning)

tb.parameters.COND_CACHE_SLOTS = 1024  # Maximum number of conditions for table queries to be kept in memory.
tb.parameters.METADATA_CACHE_SIZE = 10485760  # Size (in bytes) of the HDF5 metadata cache.
#    NODE_CACHE_SLOTS: Maximum number of nodes to be kept in the metadata cache.
#    It is the number of nodes to be kept in the metadata cache. Least recently used nodes are unloaded from memory when this number of loaded nodes is reached. To load a node again, simply access it as usual. Nodes referenced by user variables and, in general, all nodes that are still open are registered in the node manager and can be quickly accessed even if they are not in the cache.
#    Negative value means that all the touched nodes will be kept in an internal dictionary. This is the faster way to load/retrieve nodes. However, and in order to avoid a large memory comsumption, the user will be warned when the number of loaded nodes will reach the -NODE_CACHE_SLOTS value.
#    Finally, a value of zero means that any cache mechanism is disabled.
tables.parameters.NODE_CACHE_SLOTS = 4096
nb.set_num_threads(8)

stattables = [
    'avg',
    'std',
    'ste',
    'med',
    'max',
    'min',
]


# constant_cols = [
#     'iter', 'step', 'position', 'num','num_iters','max_iters'
#     'chi_lim_init', 'chi_lim_max','chi_lim',
#     'bond_dimension_current', 'bond_dimension_max',
#     'phys_time', 'algo_type', 'algo_stop']

@nb.njit(nogil=False, parallel=True, cache=True)
def get_stat(data, colname, statkey):
    result = 0.0
    if statkey == 'avg':
        result = np.mean(data)
    elif statkey == 'std':
        result = np.std(data)
    elif statkey == 'ste':
        num = np.shape(data)[0]
        result = np.std(data) / np.sqrt(num)
    elif statkey == 'med':
        result = np.median(data)
    elif statkey == 'max':
        result = np.max(data)
    elif statkey == 'min':
        result = np.min(data)
    return result


@nb.njit(nogil=False, parallel=True, cache=True)
def get_stats(data):
    num = np.shape(data)[0]
    # if not num:
    #     return 0, 0, 0, 0, 0, 0
    std = np.std(data)
    avg = np.mean(data)
    med = np.median(data)
    min = np.min(data)
    max = np.max(data)
    return avg, std, std / np.sqrt(num), med, max, min


def get_dtype(tablenode, req_columns, num=True, cast=True):
    names = []
    formats = []
    if num:
        names = ['num']
        formats = ['i8']
    for col, (dtype, offset) in tablenode.dtype.fields.items():
        if req_columns == 'ALL' or col in req_columns:
            names.append(col)
            if dtype.base.kind == 'V' and dtype.alignment == 1 and dtype.name == 'void128':  # Check if complex. These show up as "V1"
                formats.append('complex128')
            elif np.issubdtype(dtype, np.uint8) or not cast:  # np.uint8 is a bool!
                formats.append(('{}={}{}'.format(dtype.shape, dtype.base.kind, dtype.alignment)))
            else:
                formats.append(('{}=f{}'.format(dtype.shape, dtype.alignment)))  # This is for averages. Averaging ints, for instance, gives floats
    return np.dtype({"names": names, "formats": formats})


def write_statistics_table(meta, props, h5_tgt):
    # Props contains the names of the tables
    # We should generate the statistics for each column in the table
    table_name = meta[0]
    table_path = meta[1]
    table_dset = meta[2]
    for table_col_name in props[table_name]:
        data = table_dset[table_col_name]
        num = np.shape(data)[0]
        std = np.nanstd(data, axis=0)
        tgt_path = table_path + "/" + table_col_name
        tgt_node = h5_tgt.require_group(tgt_path)
        tgt_node.create_dataset(name='avg', data=np.nanmean(data, axis=0))
        tgt_node.create_dataset(name='std', data=std)
        tgt_node.create_dataset(name='ste', data=std / np.sqrt(np.shape(data)[0]))
        tgt_node.create_dataset(name='med', data=np.nanmedian(data, axis=0))
        tgt_node.create_dataset(name='max', data=np.nanmax(data, axis=0))
        tgt_node.create_dataset(name='min', data=np.nanmin(data, axis=0))
        tgt_node.create_dataset(name='q25', data=np.nanpercentile(data, 0.25, axis=0))
        tgt_node.create_dataset(name='q75', data=np.nanpercentile(data, 0.75, axis=0))
        tgt_node.create_dataset(name='num', data=num)
        tgt_node.create_dataset(name='data', data=data)
        # Generate some nice attributes for later plotting
    model_path = table_path[0: table_path.index("DMRG")] + "DMRG/model"
    if model_path in table_dset.file:
        model_node = table_dset.file[model_path]
        attrs_path = model_node.parent.parent.name
        attrs_node = h5_tgt.require_group(attrs_path)
        for key, dset in model_node.items():
            if '.db' in key:
                continue
            # print(key,dset,model_node)
            attrs_node.attrs.create(name=key, data=dset)
    else:
        raise LookupError("Could not find [" + model_path + "] in ", table_dset.file)


def write_statistics_table2(nodemeta, tablereqs, tgt):
    # Props contains the names of the tables as keys, and each value contains either "ALL" or a list of desired columns
    # We should generate the statistics for each column in the table, averaging each time point over all realizations
    t_tot = timer()
    t_stat = 0
    t_crt = 0
    t_pre = 0
    t_itr = 0
    t_app = 0
    tablename = nodemeta[0]  # The name of the table
    tablepath = nodemeta[1]  # The path of the table e.g. '...u_3/r_8/fLBIT/state_real/tables/iter_1'
    tablenode = nodemeta[2]  # The node of the table

    if not isinstance(tablenode, h5py.Dataset):
        raise TypeError("write_statistics_table2: meta should point to a h5py.Dataset. Got: ", nodemeta)
    point_node = tablenode.parent
    point_path = point_node.name
    fes = '/bond_' in point_path
    if fes:
        state_size = int(tablepath.split('L_')[1].split('/')[0])
        bond_fes = int(point_path.split('bond_')[1].split('/')[0])
        bond_max = 2 ** (state_size / 2)
        if (bond_fes > bond_max):
            return
        point_path = '{}'.format(point_node.parent.name)

    # open h5 file for appending
    with tb.File(tgt, 'a') as h5f:
        t_pre_start = timer()
        # We now have a source table. Let's get the final dtype
        dt_new = get_dtype(tablenode=tablenode, req_columns=tablereqs[tablename], num=True)
        dt_old = get_dtype(tablenode=tablenode, req_columns=tablereqs[tablename], num=False, cast=False)
        t_pre = t_pre + timer() - t_pre_start
        num = len(tablenode)
        numarray = np.array(num).reshape(1, 1)
        statgroup = '{}/{}'.format(point_path, tablename)

        # Copying the data itself
        if not fes and not '{}/data'.format(statgroup) in h5f:
            a = h5f.create_table(statgroup, name='data', title='{} data'.format(tablename),
                                 description=dt_old.newbyteorder('<'),
                                 createparents=True,
                                 expectedrows=num,
                                 filters=tb.Filters(3),
                                 track_times=False)
            data = tablenode.fields(dt_old.fields.keys())[()]
            a.append(data)

        # Create a dataset with just the number of realizations
        if not '{}/num'.format(statgroup) in h5f:
            h5f.create_earray(where=statgroup, name='num', title='Number of realizations',
                              obj=numarray,
                              createparents=True)
        else:
            a = h5f.get_node(statgroup, name='num', classname="Array")
            a.append(numarray)
            a.flush()

        statrows = {}
        for statkey in stattables:
            t_crt_start = timer()
            statpath = '{}/{}/{}'.format(point_path, tablename, statkey)
            stattitle = '{} {}'.format(tablename, statkey)
            dt = dt_old if statkey == 'data' else dt_new
            if not statpath in h5f:
                a = h5f.create_table(statgroup, name=statkey, title=stattitle,
                                     description=dt.newbyteorder('<'),
                                     createparents=True,
                                     expectedrows=num,
                                     filters=tb.Filters(3),
                                     track_times=False)
            else:
                a = h5f.get_node(statgroup, name=statkey, classname="Table")
            t_crt = t_crt + (timer() - t_crt_start)
            t_itr_start = timer()
            # create dataset row iterator
            a_row = a.row
            for r in a_row:
                r.append()
            statrows[statkey] = {'it': a_row, 'tb': a}

            t_itr = t_itr + timer() - t_itr_start

        t_stat_start = timer()

        for col, (dtype, offset) in dt_new.fields.items():
            # print("getting stats for table {} | col {} | dtype {}".format(tablepath, col,dtype))
            if col == 'num':
                stats = np.full(6, num)
            elif col == 'delta_t' or np.issubdtype(np.complex128, dtype):  # For constant complex values
                stats = np.full(6, tablenode.fields(col)[0].view(dtype=np.complex128))
            else:
                stats = get_stats(data=tablenode.fields(col)[()])

            for stat_tgt, stat_src in zip(statrows.values(), stats):
                stat_tgt['it'][col] = stat_src

        t_stat = t_stat + (timer() - t_stat_start)
        t_app_start = timer()
        for key, stat in statrows.items():
            stat['it'].append()
            stat['tb'].flush()
        t_app = t_app + timer() - t_app_start
    h5_tgt = h5open(tgt, 'a')
    model_path = tablepath[0: tablepath.index("DMRG")] + "DMRG/model"
    if model_path in tablenode.file:
        model_node = tablenode.file[model_path]
        attrs_path = model_node.parent.parent.name
        attrs_node = h5_tgt.require_group(attrs_path)
        for key, dset in model_node.items():
            if '.db' in key:
                continue
            if key in attrs_node.attrs:
                break
            # print(key,dset,model_node)
            attrs_node.attrs.create(name=key, data=dset)
    else:
        raise LookupError("Could not find [" + model_path + "] in ", tablenode.file)
    h5close(h5_tgt)

    t_tot = timer() - t_tot
    print('table2: {} | '
          'tot {:8.3e} | '
          'pre {:8.3e} {:.1f}% | '
          'crt {:8.3e} {:.1f}% | '
          'itr {:8.3e} {:.1f}% | '
          'stat {:8.3e} {:.1f}% | '
          'app {:8.3e} {:.1f}% | '
          .format(tablepath,
                  t_tot,
                  t_pre,
                  t_pre / t_tot * 100,
                  t_crt,
                  t_crt / t_tot * 100,
                  t_itr,
                  t_itr / t_tot * 100,
                  t_stat,
                  t_stat / t_tot * 100,
                  t_app,
                  t_app / t_tot * 100
                  ))
