import h5py

from .h5ops import *
from numbers import Integral
from line_profiler import profile

def get_bond_info(statenode, datanode):
    bavg_key = None
    bmax_key = None
    bavg = None
    bmax = None

    if maxnode := statenode.get('tables').get('status').get('max'):
        if 'bond_lim' in maxnode.dtype.fields.keys():
            bmax_key = 'bond_lim'
            bmax = np.max(maxnode[bmax_key][()])

    if avgnode := datanode.get('avg'):
        if 'bond_mid' in avgnode.dtype.fields.keys():
            bavg_key = 'bond_mid'
            bavg = np.max(avgnode[bavg_key][()])
    if bavg_key is None:
        if mmntnode := datanode.parent.get('measurements'):
            return get_bond_info(statenode, mmntnode)
        if crononode := statenode.get('cronos'):
            if mmntnode := crononode.get('measurements'):
                return get_bond_info(statenode, mmntnode)


    if bavg is None:
        print(f'WARNING: Could not find bavg in\n{statenode=}\nor\n{datanode=}')
        # raise LookupError(f'Could not find bavg in\n{statenode=}\nor\n{datanode=}')
    if bmax is None:
        print(f'WARNING: Could not find bmax in\n{statenode=}\nor\n{datanode=}')
        # raise LookupError(f'Could not find bmax in\n{statenode=}\nor\n{datanode=}')
    return bavg_key, bavg, bmax_key, bmax

def get_num_info(statenode, datanode):
    if avgnode := statenode.get('tables').get('status').get('avg'):
        return np.min(avgnode['num'][()])
    if isinstance(datanode, h5py.Dataset):
        if 'num' in datanode.dtype.fields:
            num = datanode['num']
            return np.min(num[()])
        if 'num' in datanode.parent:
            return get_num_info(statenode, datanode.parent)

    elif isinstance(datanode, h5py.Group):
        if 'num' in datanode and isinstance(datanode['num'], h5py.Dataset):
            num = datanode['num'][()]
            if np.isscalar(num):
                return num
            else:
                return np.min(num[()])
        if avgnode := datanode.get('avg'):
            if 'num' in avgnode.dtype.fields:
                return np.min(avgnode['num'][()])
        if maxnode := datanode.get('max'):
            if 'num' in maxnode.dtype.fields:
                return np.min(maxnode['num'][()])
        if data := datanode.get('data'):
            if len(np.shape(data)) == 1:
                return np.shape(data)[0]
            else:
                return np.shape(data)[1]
    raise LookupError(f'Failed to find num in {datanode=}')

def get_algorithm_time(statenode, datanode):
    if avgnode := statenode.get('tables').get('status').get('avg'):
        return np.mean(avgnode['algo_time'][()])/ 60
    if isinstance(datanode, h5py.Dataset):
        print(f'Looking for algorithm time in dataset: {datanode.name}')
        if 'algo_time' in datanode.dtype.fields:
            return np.max(datanode['algo_time'][()]) / 60
        elif 'algorithm_time' in datanode.dtype.fields:
            return np.max(datanode['algorithm_time'][()]) / 60
        elif 'status' in datanode.parent:
            return get_algorithm_time(statenode,datanode.parent['status'])
        elif 'status' in datanode.parent.parent:
            return get_algorithm_time(statenode,datanode.parent.parent['status'])
        elif 'measurements' in datanode.parent:
            return get_algorithm_time(statenode,datanode.parent['measurements'])
        elif 'measurements' in datanode.parent.parent:
            return get_algorithm_time(statenode,datanode.parent.parent['measurements'])
        elif 'avg' in datanode.parent:
            return get_algorithm_time(statenode,datanode.parent['avg'])

    elif isinstance(datanode, h5py.Group):
        # Check some typical datasets inside this group
        print(f'Looking for algorithm time in group: {datanode.name}')
        for link in ['avg', 'measurements', 'data', 'status']:
            if link in datanode and isinstance(datanode[link], h5py.Dataset):
                if 'algo_time' in datanode[link].dtype.fields:
                    return get_algorithm_time(statenode,datanode[link])
                elif 'algorithm_time' in datanode[link].dtype.fields:
                    return get_algorithm_time(statenode,datanode[link])
            elif link in datanode and isinstance(datanode[link], h5py.Group):
                print(f' Looking for algorithm time in group ({link}): {datanode[link].name}')
                if 'avg' in datanode[link] and isinstance(datanode[link]['avg'], h5py.Dataset):
                    if 'algo_time' in datanode[link]['avg'].dtype.fields:
                        return get_algorithm_time(statenode,datanode['avg'])
                    if 'algorithm_time' in datanode[link]['avg'].dtype.fields:
                        return get_algorithm_time(statenode,datanode[link]['avg'])
            elif link in datanode.parent:
                return get_algorithm_time(statenode,datanode.parent[link])

    raise LookupError(f'Failed to find algorithm time in {datanode.name=}')




def match_path(path, match):
    return path[0:path.index('/', path.index(match))]


def sort_db_vals(db, key):
    if len(db['vals'][key]) <= 1:
        return np.asarray(db['vals'][key])
    sort = np.argsort(np.asarray(db['vals'][key]))
    print(sort, np.asarray(db['vals'][key]))
    return np.asarray(db['vals'][key])[sort]

def find_between(s, first, last):
    try:
        start = s.index(first) + len(first)
        end = s.index(last, start)
        return s[start:end]
    except ValueError:
        return ""

def get_value(key_candidates, node_candidates, optional = True, default=None):
    if isinstance(key_candidates, str):
        keys = [key_candidates]
    elif isinstance(key_candidates, list):
        keys = key_candidates
    else:
        raise TypeError(f'key_candidates must be str or list: got {type(key_candidates)}: {key_candidates}')

    for key in keys:
        for node in node_candidates:
            if node is None:
                continue
            if key in node.attrs:
                return node.attrs[key][()]
            if isinstance(node, h5py.Dataset):
                if key in node.dtype.fields:
                    return node[key][0] # Just the first element
            if isinstance(node, h5py.Group):
                if key in node:
                    val = node[key][()]
                    if isinstance(val, np.ndarray):
                        return val[0]
                    return val
    if not optional:
        raise LookupError(f'could not find keys {keys} among node candidates: [{node_candidates}]')
    return default

def get_enum(key, node_candidates, choices, optional = True, default=None):
    val = get_value(key, node_candidates=node_candidates, optional=optional, default=default)
    if isinstance(val, Integral):
        return choices[val]
    if isinstance(val, str):
        return choices[val]

    if not optional:
        raise LookupError(f'could not match key [{key}:{val}] to enum choices: [{choices}]')
    return val

@profile
def load_isingmajorana_database(h5_src, meta, algo_filter=None, model_filter=None, state_filter=None, debug=False):
    # Gather data

    db = {
        'version': 1,
        'max': {'all': None, 'avg': None, 'mid': None, 'mvg': None},
        'min': {'all': None, 'avg': None, 'mid': None, 'mvg': None},
        'num': {'min': None, 'max': None},

        'tex': {
            'L': '$L$',
            'g': '$g$',
            'd': '$\Delta$',
            'algo': 'algo',
            'state': 'state',
            'data': 'data',
            'num': '$n$',
            'bavg': '$\\bar \chi$',
            'bmax': '$\lceil \chi \\rceil$',
            'tsim': '$\\bar t_\mathrm{sim}$',
        },
        'vals': {
            'L': set(),
            'g': set(),
            'd': set(),
            'ø': set()
        },
        'nodes': {},
        'filename': h5_src.filename,
        'plotdir': meta['common']['plotdir'],
        'cachedir': meta['common']['cachedir']
    }

    if debug:
        print('Looking for a dmrg algorithm group (depth:20) in {}'.format(h5_src.filename))

    # Collect the values that were asked for in meta
    include = {
        'L': set(),
        'g': set(),
        'd': set(),
        'ø': set(),
    }
    for mkey, mval in meta.items():
        if incl := mval.get('include'):
            if L := incl['L']: include['L'].update(L)
            if g := incl['g']: include['g'].update(g)
            if d := incl['d']: include['d'].update(d)
            if 'ø' in incl:
                if ø := incl['ø']: include['ø'].update(ø)

    print(f"{include=}")
    for algokey, algopath, algonode in h5py_node_iterator(node=h5_src, keypattern='DMRG', dep=20,
                                                          excludeKeys=['.db', 'TEBD', 'LBIT'],
                                                          nodeType=h5py.Group, godeeper=False):
        if debug:
            print(' Found algo {}'.format(algopath))
        modelnode = algonode['model']
        modelpath = modelnode.name
        hamiltonian = modelnode['hamiltonian']
        L = get_value('model_size', [modelnode, hamiltonian], optional=False)
        g = get_value(['g'], [hamiltonian], optional=False)
        d = get_value(['delta'], [hamiltonian], optional=False)
        ø = None
        # Skip if this point in the phase diagram was not asked for
        if include['L'] and not L in include['L']: continue
        if include['g'] and not g in include['g']: continue
        if include['d'] and not d in include['d']: continue

        db['vals']['L'].add(L)
        db['vals']['g'].add(g)
        db['vals']['d'].add(d)
        db['vals']['ø'].add(None)


        if debug:
            print(' Looking for datasets in meta: {}'.format(meta.keys()))
        for nodekey, nodepath, node in h5py_node_iterator(node=algonode,
                                                              excludeKeys=list(db['nodes'].keys()) + ['.db'],
                                                              dep=3):
            # Figure out if we have asked for this node
            match = False
            for mval in meta.values():
                bname = mval.get('groupbase')
                gname = mval.get('groupname')
                dname = mval.get('dsetname')
                if bname is not None and not bname in nodepath: # Asked for a specific base but it is not present
                    continue
                if gname is not None and dname is not None:
                    # This meta requires both
                    # Make sure we have gname/dname in the nodepath
                    if not f'{gname}/{dname}' in nodepath:
                        continue
                    if gname in nodepath and dname == nodekey and isinstance(node, h5py.Dataset):
                        match = True
                        break
                elif gname is None and dname is not None:
                    if dname == nodekey and isinstance(node, h5py.Dataset):
                        match = True
                        break
                elif gname is not None and dname is None:
                    if gname == nodekey and isinstance(node, h5py.Group):
                        match = True
                        break

            # Keep this node!
            if not match:
                continue
            print(f'Loading into database: {nodepath=}')
            db['nodes'][nodepath] = {}
            db['nodes'][nodepath]['version'] = db['version']
            db['nodes'][nodepath]['type'] = 'dataset' if isinstance(node, h5py.Dataset) else 'group'
            db['nodes'][nodepath]['vals'] = {}
            db['nodes'][nodepath]['vals']['L'] = L
            db['nodes'][nodepath]['vals']['g'] = g
            db['nodes'][nodepath]['vals']['d'] = d
            db['nodes'][nodepath]['vals']['ø'] = ø
            db['nodes'][nodepath]['vals']['plotdir'] = meta['common']['plotdir']
            db['nodes'][nodepath]['vals']['cachedir'] = meta['common']['cachedir']
            db['nodes'][nodepath]['vals']['filename'] = h5_src.filename
            db['nodes'][nodepath]['algonode'] = algonode
            db['nodes'][nodepath]['modelnode'] = modelnode
            db['nodes'][nodepath]['node'] = node
            db['nodes'][nodepath]['tex'] = {}
            db['nodes'][nodepath]['tex'] = {
                'keys': db['tex'],
                'vals': {
                    'L': f'{L}',
                    'g': f'{g}',
                    'd': f'{d}',
                    'ø': f'{ø}',
                    'algo': algokey,
                },
                'eqs': {
                    'L': '${}{}{}$'.format(db['tex']['L'].strip('$'), '{:}', L),
                    'g': '$g:{:.3f}$'.format(g),
                    'd': '$d:{:.2f}$'.format(d),
                    'ø': 'None',
                },
            }




    # Sort the keys so that we can iterate through them in order
    db['vals']['L'] = sorted(db['vals']['L'])
    db['vals']['g'] = sorted(db['vals']['g'])
    db['vals']['d'] = sorted(db['vals']['d'])
    return db

def load_database(h5_src, dsetname='', algo_filter='', state_filter='', debug=False):
    # Gather data
    print("Loading database:", dsetname)
    db = {
        'version': 1,
        'max': {'all': None, 'avg': None, 'mid': None, 'mvg': None},
        'min': {'all': None, 'avg': None, 'mid': None, 'mvg': None},
        'keys': {
            'L': [],
            'lambda': [],
            'delta': [],
            'algo': [],
            'state': [],
            'data': [],
        },
        'vals': {
            'L': [],
            'lambda': [],
            'delta': [],
        },
        'dsets': {},
    }
    for sizekey, sizenode in h5_src.items():
        for lambdakey, lambdanode in sizenode.items():
            for deltakey, deltanode in lambdanode.items():
                for algokey, algopath, algonode in h5py_group_iterator(node=deltanode, keypattern=algo_filter, dep=1):
                    for statekey, statepath, statenode in h5py_group_iterator(node=algonode, keypattern=state_filter, dep=1):
                        if debug:
                            print("Starting gather")
                        for datakey, datapath, datanode in h5py_node_finder(g=statenode,
                                                                            keypattern=dsetname,
                                                                            dep=3):
                            if datanode.name in db['dsets']:
                                continue
                            if debug:
                                print(sizekey, lambdakey, deltakey, algokey, statekey, datakey)
                            size = deltanode.attrs['model_size']
                            lamb = deltanode.attrs['lambda']
                            delt = deltanode.attrs['delta']
                            if debug:
                                print("Adding keys")
                            db['keys']['L'].append(sizekey) if not sizekey in db['keys']['L'] else db['keys']['L']
                            db['keys']['lambda'].append(lambdakey) if not lambdakey in db['keys']['lambda'] else db['keys']['lambda']
                            db['keys']['delta'].append(deltakey) if not deltakey in db['keys']['delta'] else db['keys']['delta']
                            db['keys']['algo'].append(algokey) if not algokey in db['keys']['algo'] else db['keys']['algo']
                            db['keys']['state'].append(statekey) if not statekey in db['keys']['state'] else db['keys']['state']
                            db['keys']['data'].append(datakey) if not datakey in db['keys']['data'] else db['keys']['data']
                            if debug:
                                print("Adding vals")
                            db['vals']['L'].append(size) if not size in db['vals']['L'] else db['vals']['L']
                            db['vals']['lambda'].append(lamb) if not lamb in db['vals']['lambda'] else db['vals']['lambda']
                            db['vals']['delta'].append(delt) if not delt in db['vals']['delta'] else db['vals']['delta']

                            if debug:
                                print("Adding extrema")
                            avgshape = np.shape(datanode['avg'][()])
                            datashape = np.shape(datanode['data'][()])
                            midx = int(avgshape[0] / 2) if avgshape or len(avgshape) > 1 else ()
                            db['max']['all'] = np.max([db['max']['all'], np.max(datanode['max'][()])]) if db['max']['all'] else np.max(datanode['max'][()])
                            db['max']['avg'] = np.max([db['max']['avg'], np.max(datanode['avg'][()])]) if db['max']['avg'] else np.max(datanode['avg'][()])
                            db['min']['all'] = np.min([db['min']['all'], np.min(datanode['min'][()])]) if db['min']['all'] else np.min(datanode['min'][()])
                            db['min']['avg'] = np.min([db['min']['avg'], np.min(datanode['avg'][()])]) if db['min']['avg'] else np.min(datanode['avg'][()])
                            db['max']['mid'] = np.max([db['max']['all'], np.max(datanode['max'][midx])]) if db['max']['all'] else np.max(datanode['max'][midx])
                            db['max']['mvg'] = np.max([db['max']['mvg'], np.max(datanode['avg'][midx])]) if db['max']['mvg'] else np.max(datanode['avg'][midx])
                            db['min']['mid'] = np.min([db['min']['all'], np.min(datanode['min'][midx])]) if db['min']['all'] else np.min(datanode['min'][midx])
                            db['min']['mvg'] = np.min([db['min']['mvg'], np.min(datanode['avg'][midx])]) if db['min']['mvg'] else np.min(datanode['avg'][midx])

                            if debug:
                                print("Adding node specific data")
                            db['dsets'][datanode.name] = {}
                            db['dsets'][datanode.name]['keys'] = {'L': sizekey, 'lambda': lambdakey,
                                                                  'delta': deltakey, 'algo': algokey,
                                                                  'state': statekey, 'data': datakey}

                            db['dsets'][datanode.name]['num'] = datanode['num'][()]
                            db['dsets'][datanode.name]['L'] = size
                            db['dsets'][datanode.name]['lambda'] = lamb
                            db['dsets'][datanode.name]['delta'] = delt
                            db['dsets'][datanode.name]['sizenode'] = sizenode
                            db['dsets'][datanode.name]['lambdanode'] = lambdanode
                            db['dsets'][datanode.name]['deltanode'] = deltanode
                            db['dsets'][datanode.name]['statenode'] = statenode
                            db['dsets'][datanode.name]['datanode'] = datanode
                            db['dsets'][datanode.name]['avgshape'] = avgshape
                            db['dsets'][datanode.name]['datashape'] = datashape
                            db['dsets'][datanode.name]['midx'] = midx
                            if debug:
                                print("Adding styles")

                            db['dsets'][datanode.name]['style'] = {}
                            if "states" in statekey:
                                db['dsets'][datanode.name]['style']['lwidth'] = 3.0
                                db['dsets'][datanode.name]['style']['lalpha'] = 1.0
                                db['dsets'][datanode.name]['style']['mstyle'] = None
                                db['dsets'][datanode.name]['style']['lstyle'] = 'solid'
                            else:
                                db['dsets'][datanode.name]['style']['lwidth'] = 1.4
                                db['dsets'][datanode.name]['style']['lalpha'] = 1.0
                                db['dsets'][datanode.name]['style']['mstyle'] = '.'
                                db['dsets'][datanode.name]['style']['lstyle'] = 'dotted'
                            if debug:
                                print("Finished gather")


    return db
