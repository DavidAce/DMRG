from dmrg_plot.common.io.h5ops import *
from plotting.filter import *
import itertools


def get_bond_info(statenode, datanode):
    bavg_key = None
    bmax_key = None
    bavg = None
    bmax = None
    try:
        bavg_key = next(k for k in ['bond_mid'] if k in datanode['avg'].dtype.fields.keys())
        bavg = np.max(datanode['avg'][bavg_key][()])
    except:
        pass
    try:
        bmax_key = next(k for k in ['bond_lim'] if k in statenode['tables']['status']['max'].dtype.fields.keys())
        bmax = np.max(statenode['tables']['status']['max'][bmax_key][()])
    except:
        pass
    return bavg_key, bavg, bmax_key, bmax


def load_time_database(h5_src, dsetname, algo_filter=None, state_filter=None, point_filter=None, debug=False):
    # Gather data
    print("Loading time database:", dsetname)
    db = {
        'max': {'all': None, 'avg': None, 'mid': None, 'mvg': None},
        'min': {'all': None, 'avg': None, 'mid': None, 'mvg': None},
        'num': {'min': None, 'max': None},
        'keys': {
            'L': [], 'J': [], 'w': [], 'b': [], 'f': [], 'u': [], 'r': [],
            'algo': [],
            'state': [],
            'point': [],
            'data': [],
        },
        'vals': {
            'L': [], 'J': [], 'w': [], 'b': [], 'f': [], 'u': [], 'r': []
        },
        'dsets': {},
    }
    for Lkey, Lnode in sorted(h5_src.items(), key=natural_items):
        for Jkey, Jnode in sorted(Lnode.items(), key=natural_items):
            for wkey, wnode in sorted(Jnode.items(), key=natural_items):
                for bkey, bnode in sorted(wnode.items(), key=natural_items):
                    for fkey, fnode in sorted(bnode.items(), key=natural_items):
                        for ukey, unode in sorted(fnode.items(), key=natural_items):
                            for rkey, rnode in sorted(unode.items(), key=natural_items):
                                for algokey, algopath, algonode in h5py_group_iterator(node=rnode, keypattern=algo_filter, dep=1):
                                    for statekey, statepath, statenode in h5py_group_iterator(node=algonode, keypattern=state_filter, dep=1):
                                        for pointkey, pointpath, pointnode in h5py_group_iterator(node=statenode,
                                                                                                  keypattern=point_filter,
                                                                                                  dep=1):
                                            # if not Jkey == "J[+0.00_+0.00_+0.00]":
                                            #     continue
                                            # if not wkey == "w[+1.00_+1.00_+0.10]":
                                            #     continue
                                            # if not fkey == "f_0.30":
                                            #     continue
                                            # if not bkey == "b_5.00":
                                            #     continue
                                            # if not ukey == "u_3":
                                            #     continue
                                            # if not rkey in ["r_8"]:
                                            #     continue
                                            # if not Lkey in ["L_16", "L_18", "L_20", "L_22", "L_24"]:
                                            #     continue

                                            if debug:
                                                print("Starting gather: Finding ", dsetname)
                                            for datakey, datapath, datanode in h5py_node_iterator(node=pointnode,
                                                                                                  keypattern=dsetname,
                                                                                                  dep=6,
                                                                                                  nodeType=h5py.Group):
                                                if datanode.name in db['dsets']:
                                                    continue
                                                if debug:
                                                    print(Lkey, Jkey, wkey, bkey, fkey, ukey, algokey, statekey,
                                                          pointkey, datakey)

                                                L = rnode.attrs['model_size']
                                                J = [rnode.attrs['J1_mean'], rnode.attrs['J2_mean'],
                                                     rnode.attrs['J3_mean']]
                                                w = [rnode.attrs['J1_wdth'], rnode.attrs['J2_wdth'],
                                                     rnode.attrs['J3_wdth']]
                                                b = rnode.attrs['J2_base']
                                                f = rnode.attrs['u_fmix']
                                                u = rnode.attrs['u_depth']
                                                r = rnode.attrs['J2_span']

                                                if debug:
                                                    print("Adding keys")
                                                db['keys']['L'].append(Lkey) if not Lkey in db['keys']['L'] else \
                                                db['keys']['L']
                                                db['keys']['J'].append(Jkey) if not Jkey in db['keys']['J'] else \
                                                db['keys']['J']
                                                db['keys']['w'].append(wkey) if not wkey in db['keys']['w'] else \
                                                db['keys']['w']
                                                db['keys']['b'].append(bkey) if not bkey in db['keys']['b'] else \
                                                db['keys']['b']
                                                db['keys']['f'].append(fkey) if not fkey in db['keys']['f'] else \
                                                db['keys']['f']
                                                db['keys']['u'].append(ukey) if not ukey in db['keys']['u'] else \
                                                db['keys']['u']
                                                db['keys']['r'].append(rkey) if not rkey in db['keys']['r'] else db['keys']['r']
                                                db['keys']['algo'].append(algokey) if not algokey in db['keys']['algo'] else db['keys']['algo']
                                                db['keys']['state'].append(statekey) if not statekey in db['keys']['state'] else db['keys']['state']
                                                db['keys']['point'].append(pointkey) if not pointkey in db['keys']['point'] else db['keys']['point']
                                                db['keys']['data'].append(datakey) if not datakey in db['keys']['data'] else db['keys']['data']
                                                if debug:
                                                    print("Adding vals")
                                                db['vals']['L'].append(L) if not L in db['vals']['L'] else db['vals']['L']
                                                db['vals']['J'].append(J) if not J in db['vals']['J'] else db['vals']['J']
                                                db['vals']['w'].append(w) if not w in db['vals']['w'] else db['vals']['w']
                                                db['vals']['b'].append(b) if not b in db['vals']['b'] else db['vals']['b']
                                                db['vals']['f'].append(f) if not f in db['vals']['f'] else db['vals']['f']
                                                db['vals']['u'].append(u) if not u in db['vals']['u'] else db['vals']['u']
                                                db['vals']['r'].append(r) if not r in db['vals']['r'] else db['vals']['r']

                                                if debug:
                                                    print("Adding extrema")
                                                # Iterate through all the time points to find extrema
                                                timesteps = len(datanode.keys())
                                                for tkey, tnode in itertools.islice(datanode.items(), 0, None, int(np.max([50, timesteps]) / 20)):
                                                    max = tnode['max'][()]
                                                    min = tnode['min'][()]
                                                    avg = tnode['avg'][()]
                                                    num = tnode['num'][()]
                                                    avgshape = np.shape(avg)
                                                    datashape = np.shape(tnode['data'][()])
                                                    midx = int(avgshape[0] / 2) if avgshape or len(avgshape) > 1 else ()
                                                    # db['max']['all'] = np.max([db['max']['all'], np.max(max)])   if db['max']['all'] else np.max(max)
                                                    # db['max']['avg'] = np.max([db['max']['avg'], np.max(avg)])   if db['max']['avg'] else np.max(avg)
                                                    # db['min']['all'] = np.min([db['min']['all'], np.min(min)])   if db['min']['all'] else np.min(min)
                                                    # db['min']['avg'] = np.min([db['min']['avg'], np.min(avg)])   if db['min']['avg'] else np.min(avg)
                                                    # db['max']['mid'] = np.max([db['max']['mid'], max[midx]])     if db['max']['mid'] else max[midx]
                                                    # db['max']['mvg'] = np.max([db['max']['mvg'], avg[midx]])     if db['max']['mvg'] else avg[midx]
                                                    # db['min']['mid'] = np.min([db['min']['all'], min[midx]])     if db['min']['all'] else min[midx]
                                                    # db['min']['mvg'] = np.min([db['min']['mvg'], avg[midx]])     if db['min']['mvg'] else avg[midx]
                                                    db['num']['min'] = np.min([db['num']['min'], num]) if db['num']['min'] else np.min(num)
                                                    db['num']['max'] = np.max([db['num']['max'], num]) if db['num']['max'] else np.min(num)

                                                if debug:
                                                    print("Adding node specific data")
                                                db['dsets'][datanode.name] = {}
                                                db['dsets'][datanode.name]['keys'] = {'L': Lkey, 'J': Jkey,
                                                                                      'w': wkey, 'b': bkey, 'f': fkey, 'u': ukey,
                                                                                      'algo': algokey, 'state': statekey, 'point': pointkey,
                                                                                      'data': datakey}

                                                db['dsets'][datanode.name]['L'] = L
                                                db['dsets'][datanode.name]['J'] = J
                                                db['dsets'][datanode.name]['w'] = w
                                                db['dsets'][datanode.name]['b'] = b
                                                db['dsets'][datanode.name]['f'] = f
                                                db['dsets'][datanode.name]['u'] = u
                                                db['dsets'][datanode.name]['r'] = r
                                                db['dsets'][datanode.name]['Lnode'] = Lnode
                                                db['dsets'][datanode.name]['Jnode'] = Jnode
                                                db['dsets'][datanode.name]['wnode'] = wnode
                                                db['dsets'][datanode.name]['bnode'] = bnode
                                                db['dsets'][datanode.name]['fnode'] = fnode
                                                db['dsets'][datanode.name]['unode'] = unode
                                                db['dsets'][datanode.name]['statenode'] = statenode
                                                db['dsets'][datanode.name]['pointnode'] = pointnode
                                                # db['dsets'][datanode.name]['mmntnode'] = pointnode['measurements']
                                                db['dsets'][datanode.name]['datanode'] = datanode
                                                db['dsets'][datanode.name]['avgshape'] = avgshape
                                                db['dsets'][datanode.name]['datashape'] = datashape
                                                db['dsets'][datanode.name]['midx'] = midx

                                                if debug:
                                                    print("Finished gather")

    # Sort the keys so that we can iterate through them in order
    # sortL = np.argsort(np.asarray(db['vals']['L']))
    # sortb = np.argsort(np.asarray(db['vals']['b']))
    # sortf = np.argsort(np.asarray(db['vals']['f']))
    # sortu = np.argsort(np.asarray(db['vals']['u']))
    #
    # db['vals']['L']   = np.asarray(db['vals']['L'])[sortL]
    # db['keys']['L']   = np.asarray(db['keys']['L'])[sortL]
    # db['vals']['b']  = np.asarray(db['vals']['b'])[sortb]
    # db['keys']['b']  = np.asarray(db['keys']['b'])[sortb]
    # db['vals']['f']  = np.asarray(db['vals']['f'])[sortf]
    # db['keys']['f']  = np.asarray(db['keys']['f'])[sortf]
    # db['vals']['u']  = np.asarray(db['vals']['u'])[sortu]
    # db['keys']['u']  = np.asarray(db['keys']['u'])[sortu]
    return db


def match_path(path, match):
    return path[0:path.index('/', path.index(match))]


def sort_db_vals(db, key):
    if len(db['vals'][key]) <= 1:
        return np.asarray(db['vals'][key])
    sort = np.argsort(np.asarray(db['vals'][key]))
    print(sort, np.asarray(db['vals'][key]))
    return np.asarray(db['vals'][key])[sort]


def load_time_database3(h5_src, meta, algo_filter=None, model_filter=None, state_filter=None, debug=False):
    # Gather data

    db = {
        'version': 3,
        'max': {'all': None, 'avg': None, 'mid': None, 'mvg': None},
        'min': {'all': None, 'avg': None, 'mid': None, 'mvg': None},
        'num': {'min': None, 'max': None},
        'keys': {
            'L': set(), 'J': set(), 'w': set(), 'x': set(), 'r': set(),
            'u': set(), 'f': set(), 'tstd': set(), 'cstd': set(), 'tgw8': set(), 'cgw8': set(),
            'algo': set(),
            'state': set(),
            'crono': set(),
            'data': set(),
        },
        'tex': {
            'L': '$L$',
            'J': '$J$', 'J1': '$J_1$', 'J2': '$J_2$', 'J3': '$J_3$',
            'w': '$w$', 'w1': '$w_1$', 'w2': '$w_2$', 'w3': '$w_3$',
            'x': '$\\xi_J$', 'r': '$r$',
            'u': '$d_u$', 'f': '$f$', 'tstd': '$\sigma_\\theta$', 'cstd': '$\sigma_c$', 'tgw8': '$w_\\theta$',
            'cgw8': '$w_c$',
            't': '$t$',
            'algo': 'algo',
            'state': 'state',
            'crono': 'crono',
            'data': 'data',
            'num': '$n$',
            'bavg': '$\\bar \chi$',
            'bmax': '$\lceil \chi \\rceil$',
            'time': '$\\bar t_\mathrm{sim}$',
            'tsim': '$\\bar t_\mathrm{sim}$',
        },
        'vals': {
            'L': set(),
            'J': set(), 'J1': set(), 'J2': set(), 'J3': set(),
            'w': set(), 'w1': set(), 'w2': set(), 'w3': set(),
            'x': set(), 'r': set(),
            'u': set(), 'f': set(), 'tstd': set(), 'cstd': set(), 'tgw8': set(), 'cgw8': set()
        },
        'dsets': {},
    }
    for algokey, algopath, algonode in h5py_node_iterator(node=h5_src, keypattern='fLBIT', dep=20,
                                                          excludeKeys=['.db', 'fDMRG', 'xDMRG', 'iDMRG', 'iTEBD'],
                                                          nodeType=h5py.Group, godeeper=False):
        if debug:
            print(algopath)
        modelnode = algonode['model']
        hamiltonian = modelnode['hamiltonian']
        L = modelnode['model_size'][()]
        J = hamiltonian['J1_mean', 'J2_mean', 'J3_mean'][0]
        w = hamiltonian['J1_wdth', 'J2_wdth', 'J3_wdth'][0]
        r = hamiltonian['J2_span'][0]
        x = hamiltonian['xi_Jcls'][0] if 'xi_Jcls' in hamiltonian.dtype.fields else hamiltonian['J2_xcls'][0]
        u = hamiltonian['u_depth'][0]
        f = hamiltonian['u_fmix'][0]
        tstd = hamiltonian['u_tstd'][0]
        cstd = hamiltonian['u_cstd'][0]
        tgw8 = 'ID' if hamiltonian['u_tgw8'][0] == 0 else 'EX'
        cgw8 = 'ID' if hamiltonian['u_cgw8'][0] == 0 else 'EX'
        # L if r == np.iinfo(np.uint64).max else r
        # Skip if not asked for
        if incl := meta.get('common').get('include_vals'):
            val_requested = True
            for tag, val in zip(['L', 'J', 'w', 'x', 'r', 'u', 'f', 'tstd', 'cstd', 'tgw8', 'cgw8'],
                                [L, J, w, x, r, u, f, tstd, cstd, tgw8, cgw8]):
                if tag in incl and not val in incl.get(tag):
                    val_requested = False
                    break
            if not val_requested:
                continue
            else:
                print([L, J, w, x, r, u, f, tstd, cstd, tgw8, cgw8])
        if debug:
            print("Adding keys")

        db['vals']['L'].add(L)
        db['vals']['J'].add(tuple(J))
        db['vals']['J1'].add(J[0])
        db['vals']['J2'].add(J[1])
        db['vals']['J3'].add(J[2])
        db['vals']['w'].add(tuple(w))
        db['vals']['w1'].add(w[0])
        db['vals']['w2'].add(w[1])
        db['vals']['w3'].add(w[2])
        db['vals']['x'].add(x)
        db['vals']['r'].add(r)
        db['vals']['u'].add(u)
        db['vals']['f'].add(f)
        db['vals']['tstd'].add(tstd)
        db['vals']['cstd'].add(cstd)
        db['vals']['tgw8'].add(tgw8)
        db['vals']['cgw8'].add(cgw8)

        db['keys']['algo'].add(algokey)
        for metakey, descr in meta.items():
            if 'include' in metakey or 'common' in metakey:
                continue
            if not 'dsetname' in descr:
                continue
            if debug:
                print('Looking for dataset {} in {}'.format(descr['dsetname'], modelnode.name))
            for datakey, datapath, datanode in h5py_node_iterator(node=modelnode,
                                                                  keypattern=descr['dsetname'],
                                                                  dep=4,
                                                                  nodeType=h5py.Group):
                print("Loading dset database version 3: {}".format(datapath))

                num = np.max(datanode['num'][()])
                for dname in [datanode.name]:
                    db['dsets'][dname] = {}
                    db['dsets'][dname]['version'] = db['version']
                    db['dsets'][dname]['keys'] = db['keys']
                    db['dsets'][dname]['vals'] = {}
                    db['dsets'][dname]['vals']['L'] = L
                    db['dsets'][dname]['vals']['J'] = J
                    db['dsets'][dname]['vals']['J1'] = J[0]
                    db['dsets'][dname]['vals']['J2'] = J[1]
                    db['dsets'][dname]['vals']['J3'] = J[2]
                    db['dsets'][dname]['vals']['w'] = w
                    db['dsets'][dname]['vals']['w1'] = w[0]
                    db['dsets'][dname]['vals']['w2'] = w[1]
                    db['dsets'][dname]['vals']['w3'] = w[2]
                    db['dsets'][dname]['vals']['x'] = x
                    db['dsets'][dname]['vals']['r'] = r
                    db['dsets'][dname]['vals']['u'] = u
                    db['dsets'][dname]['vals']['f'] = f
                    db['dsets'][dname]['vals']['tstd'] = tstd
                    db['dsets'][dname]['vals']['cstd'] = cstd
                    db['dsets'][dname]['vals']['tgw8'] = tgw8
                    db['dsets'][dname]['vals']['cgw8'] = cgw8
                    db['dsets'][dname]['vals']['num'] = num

                    db['dsets'][dname]['node'] = {}
                    db['dsets'][dname]['node']['L'] = h5_src[match_path(modelnode.name, 'L')]
                    db['dsets'][dname]['node']['J'] = h5_src[match_path(modelnode.name, 'J[')]
                    db['dsets'][dname]['node']['x'] = h5_src[match_path(modelnode.name, 'x')]
                    db['dsets'][dname]['node']['r'] = h5_src[match_path(modelnode.name, 'r')]
                    db['dsets'][dname]['node']['f'] = h5_src[match_path(modelnode.name, 'u[')]
                    db['dsets'][dname]['node']['u'] = h5_src[match_path(modelnode.name, 'u[')]
                    db['dsets'][dname]['node']['tstd'] = h5_src[match_path(modelnode.name, 'u[')]
                    db['dsets'][dname]['node']['cstd'] = h5_src[match_path(modelnode.name, 'u[')]
                    db['dsets'][dname]['node']['tgw8'] = h5_src[match_path(modelnode.name, 'u[')]
                    db['dsets'][dname]['node']['cgw8'] = h5_src[match_path(modelnode.name, 'u[')]
                    db['dsets'][dname]['node']['model'] = modelnode
                    db['dsets'][dname]['node']['data'] = datanode
                    db['dsets'][dname]['node']['avg'] = datanode['avg']

                    db['dsets'][dname]['tex'] = {
                        'keys': db['tex'],
                        'vals': {
                            'L': '{}'.format(L),
                            'J': '{}'.format(J),
                            'J1': '{}'.format(J[0]),
                            'J2': '{}'.format(J[1]),
                            'J3': '{}'.format(J[2]),
                            'w': '{}'.format(w),
                            'w1': '{}'.format(w[0]),
                            'w2': '{}'.format(w[1]),
                            'w3': '{}'.format(w[2]),
                            'x': '{:.2f}'.format(x),
                            'r': '$L$' if r == np.iinfo(np.uint64).max else '{}'.format(r),
                            'u': '{}'.format(u),
                            'f': '{:.2f}'.format(f),
                            'tstd': '{:.2f}'.format(tstd),
                            'cstd': '{:.2f}'.format(cstd),
                            'tgw8': '{}'.format(tgw8),
                            'cgw8': '{}'.format(cgw8),
                            'num': '{}'.format(num),
                            'algo': algokey,
                        },
                        'eqs': {
                            'L': '${}{}{}$'.format(db['tex']['L'].strip('$'), '{:}', L),
                            'J': '$J:{}$'.format(J),
                            'J1': '$J_1:{}$'.format(J[0]),
                            'J2': '$J_2:{}$'.format(J[1]),
                            'J3': '$J_3:{}$'.format(J[2]),
                            'w': '$w:{}$'.format(w),
                            'w1': '$w_1:{}$'.format(w[0]),
                            'w2': '$w_2:{}$'.format(w[1]),
                            'w3': '$w_3:{}$'.format(w[2]),
                            'x': '${}{}{:>.2f}$'.format(db['tex']['x'].strip('$'), '{:}', x),
                            'r': '$r = L$' if r == np.iinfo(np.uint64).max else '$r = {}$'.format(r),
                            'u': '${}{}{:>.0f}$'.format(db['tex']['u'].strip('$'), '{:}', u),
                            'f': '${}{}{:>.2f}$'.format(db['tex']['f'].strip('$'), '{:}', f),
                        },
                    }
        for statekey, statepath, statenode in h5py_group_iterator(node=algonode, keypattern=state_filter, dep=1):
            db['keys']['state'].add(statekey)
            for cronokey, cronopath, crononode in h5py_group_iterator(node=statenode, keypattern='cronos', dep=1):
                db['keys']['crono'].add(cronokey)

                for metakey, descr in meta.items():
                    if metakey == 'include' or metakey == 'common':
                        continue
                    objname = None
                    if 'colname' in descr:  # Consider only tables
                        objname = descr['colname']
                    elif 'dsetname' in descr:  # ... or datasets
                        objname = descr['dsetname']
                    else:
                        continue
                    print("Loading time database version 3: {}".format(cronopath))

                    for datakey, datapath, datanode in h5py_node_iterator(node=crononode,
                                                                          keypattern=descr['groupname'],
                                                                          dep=6,
                                                                          nodeType=h5py.Group):
                        if datanode.name in db['dsets']:
                            continue
                        num = 0
                        tsim = 0.0
                        bavg_key = 0.0
                        if 'colname' in descr:  # We have a table!
                            num = np.max(datanode['max']['num'][()])
                            tsim = np.max(datanode['avg']['algorithm_time'][()]) / 60
                            bavg_key, bavg, bmax_key, bmax = get_bond_info(statenode, datanode)
                        else:
                            mmntnode = datanode.parent['measurements']
                            num = np.shape(datanode['data'])[1]
                            tsim = np.max(mmntnode['avg']['algorithm_time'][()]) / 60
                            bavg_key, bavg, bmax_key, bmax = get_bond_info(statenode, mmntnode)

                        if debug:
                            print("Adding node data")
                        for dname in [datanode.name]:
                            db['dsets'][dname] = {}
                            db['dsets'][dname]['version'] = db['version']
                            db['dsets'][dname]['keys'] = db['keys']
                            db['dsets'][dname]['vals'] = {}
                            db['dsets'][dname]['vals']['L'] = L
                            db['dsets'][dname]['vals']['J'] = J
                            db['dsets'][dname]['vals']['J1'] = J[0]
                            db['dsets'][dname]['vals']['J2'] = J[1]
                            db['dsets'][dname]['vals']['J3'] = J[2]
                            db['dsets'][dname]['vals']['w'] = w
                            db['dsets'][dname]['vals']['w1'] = w[0]
                            db['dsets'][dname]['vals']['w2'] = w[1]
                            db['dsets'][dname]['vals']['w3'] = w[2]
                            db['dsets'][dname]['vals']['x'] = x
                            db['dsets'][dname]['vals']['r'] = r
                            db['dsets'][dname]['vals']['num'] = num
                            db['dsets'][dname]['vals']['u'] = u
                            db['dsets'][dname]['vals']['f'] = f
                            db['dsets'][dname]['vals']['tstd'] = tstd
                            db['dsets'][dname]['vals']['cstd'] = cstd
                            db['dsets'][dname]['vals']['tgw8'] = tgw8
                            db['dsets'][dname]['vals']['cgw8'] = cgw8
                            db['dsets'][dname]['vals']['tsim'] = tsim
                            db['dsets'][dname]['vals']['bavg'] = bavg
                            db['dsets'][dname]['vals']['bmax'] = bmax

                            db['dsets'][dname]['node'] = {}
                            db['dsets'][dname]['node']['L'] = h5_src[match_path(modelnode.name, 'L')]
                            db['dsets'][dname]['node']['J'] = h5_src[match_path(modelnode.name, 'J[')]
                            db['dsets'][dname]['node']['x'] = h5_src[match_path(modelnode.name, 'x')]
                            db['dsets'][dname]['node']['r'] = h5_src[match_path(modelnode.name, 'r')]
                            db['dsets'][dname]['node']['f'] = h5_src[match_path(modelnode.name, 'u[')]
                            db['dsets'][dname]['node']['u'] = h5_src[match_path(modelnode.name, 'u[')]
                            db['dsets'][dname]['node']['tstd'] = h5_src[match_path(modelnode.name, 'u[')]
                            db['dsets'][dname]['node']['cstd'] = h5_src[match_path(modelnode.name, 'u[')]
                            db['dsets'][dname]['node']['tgw8'] = h5_src[match_path(modelnode.name, 'u[')]
                            db['dsets'][dname]['node']['cgw8'] = h5_src[match_path(modelnode.name, 'u[')]
                            db['dsets'][dname]['node']['model'] = modelnode
                            db['dsets'][dname]['node']['data'] = datanode
                            db['dsets'][dname]['node']['avg'] = datanode['avg']
                            db['dsets'][dname]['node']['state'] = statenode
                            db['dsets'][dname]['node']['table'] = statenode.get('tables')
                            db['dsets'][dname]['node']['crono'] = crononode

                            db['dsets'][dname]['tex'] = {
                                'keys': db['tex'],
                                'vals': {
                                    'L': '{}'.format(L),
                                    'J': '{}'.format(J),
                                    'J1': '{}'.format(J[0]),
                                    'J2': '{}'.format(J[1]),
                                    'J3': '{}'.format(J[2]),
                                    'w': '{}'.format(w),
                                    'w1': '{}'.format(w[0]),
                                    'w2': '{}'.format(w[1]),
                                    'w3': '{}'.format(w[2]),
                                    'x': '{:.2f}'.format(x),
                                    'r': '$L$' if r == np.iinfo(np.uint64).max else '{}'.format(r),
                                    'u': '{}'.format(u),
                                    'f': '{:.2f}'.format(f),
                                    'tstd': '{:.2f}'.format(tstd),
                                    'cstd': '{:.2f}'.format(cstd),
                                    'tgw8': '{}'.format(tgw8),
                                    'cgw8': '{}'.format(cgw8),
                                    'num': '{}'.format(num),
                                    'tsim': '{:>.1f}m'.format(tsim),
                                    'bavg': '{}'.format(bavg),
                                    'bmax': '{}'.format(bmax),
                                    'algo': algokey,
                                    'state': statekey,
                                    'crono': cronokey,
                                    'data': datakey,
                                },
                                'eqs': {
                                    'L': '${}{}{}$'.format(db['tex']['L'].strip('$'), '{:}', L),
                                    'J': '$J:{}$'.format(J),
                                    'J1': '$J:{}$'.format(J[0]),
                                    'J2': '$J:{}$'.format(J[1]),
                                    'J3': '$J:{}$'.format(J[2]),
                                    'w': '$w:{}$'.format(w),
                                    'w1': '$w:{}$'.format(w[0]),
                                    'w2': '$w:{}$'.format(w[1]),
                                    'w3': '$w:{}$'.format(w[2]),
                                    'x': '${}{}{:>.4f}$'.format(db['tex']['x'].strip('$'), '{:}', x),
                                    'r': '$r = L$' if r == np.iinfo(np.uint64).max else '$r = {}$'.format(r),
                                    'u': '${}{}{:>.0f}$'.format(db['tex']['u'].strip('$'), '{:}', u),
                                    'f': '${}{}{:>.4f}$'.format(db['tex']['f'].strip('$'), '{:}', f),
                                    'tstd': '${}{}{:>.2f}$'.format(db['tex']['tstd'].strip('$'), '{:}', tstd),
                                    'cstd': '${}{}{:>.2f}$'.format(db['tex']['cstd'].strip('$'), '{:}', cstd),
                                    'tgw8': '{}{}{}'.format(db['tex']['tgw8'].strip('$'), '{:}', tgw8),
                                    'cgw8': '{}{}{}'.format(db['tex']['cgw8'].strip('$'), '{:}', cgw8),
                                },
                            }
                        if debug:
                            print("Finished gather")

    # Sort the keys so that we can iterate through them in order
    db['vals']['L'] = sorted(db['vals']['L'])
    db['vals']['J'] = sorted(db['vals']['J'])
    db['vals']['J1'] = sorted(db['vals']['J1'])
    db['vals']['J2'] = sorted(db['vals']['J2'])
    db['vals']['J3'] = sorted(db['vals']['J3'])
    db['vals']['w'] = sorted(db['vals']['w'])
    db['vals']['w1'] = sorted(db['vals']['w1'])
    db['vals']['w2'] = sorted(db['vals']['w2'])
    db['vals']['w3'] = sorted(db['vals']['w3'])
    db['vals']['x'] = sorted(db['vals']['x'])
    db['vals']['r'] = sorted(db['vals']['r'])
    db['vals']['u'] = sorted(db['vals']['u'])
    db['vals']['f'] = sorted(db['vals']['f'])
    db['vals']['tstd'] = sorted(db['vals']['tstd'])
    db['vals']['cstd'] = sorted(db['vals']['cstd'])
    db['vals']['tgw8'] = sorted(db['vals']['tgw8'])
    db['vals']['cgw8'] = sorted(db['vals']['cgw8'])
    return db


def load_time_database2(h5_src, meta, algo_filter=None, model_filter=None, state_filter=None, debug=False):
    # Gather data

    db = {
        'version': 2,
        'max': {'all': None, 'avg': None, 'mid': None, 'mvg': None},
        'min': {'all': None, 'avg': None, 'mid': None, 'mvg': None},
        'num': {'min': None, 'max': None},
        'keys': {
            'L': [], 'J': [], 'w': [], 'x': [], 'f': [], 'u': [], 'r': [],
            'algo': [],
            'model': [],
            'state': [],
            'crono': [],
            'data': [],
        },
        'tex': {
            'L': '$L$', 'J': '$J$', 'w': '$\omega$', 'x': '$\\xi_J$', 'f': '$f$', 'u': '$d_u$', 'r': '$r$',
            't': '$t$',
            'algo': 'algo',
            'model': 'model',
            'state': 'state',
            'crono': 'crono',
            'data': 'data',
            'num': '$n$',
            'bavg': '$\\bar \chi$',
            'bmax': '$\lceil \chi \\rceil$',
            'time': '$\\bar t_\mathrm{sim}$',
            'tsim': '$\\bar t_\mathrm{sim}$',
        },
        'vals': {
            'L': [], 'J': [], 'w': [], 'x': [], 'f': [], 'u': [], 'r': []
        },
        'dsets': {},
    }
    for Lkey, Lnode in sorted(h5_src.items(), key=natural_items):
        for Jkey, Jnode in sorted(Lnode.items(), key=natural_items):
            for wkey, wnode in sorted(Jnode.items(), key=natural_items):
                for xkey, xnode in sorted(wnode.items(), key=natural_items):
                    for fkey, fnode in sorted(xnode.items(), key=natural_items):
                        for ukey, unode in sorted(fnode.items(), key=natural_items):
                            for rkey, rnode in sorted(unode.items(), key=natural_items):

                                # Skip if not asked for
                                if incl := meta.get('common').get('include'):
                                    key_requested = True
                                    for tag, key in zip(['L', 'J', 'w', 'x', 'f', 'u', 'r'], [Lkey, Jkey, wkey, xkey, fkey, ukey, rkey]):
                                        if tag in incl and not key in incl.get(tag):
                                            key_requested = False
                                            break
                                    if not key_requested:
                                        continue
                                    else:
                                        print(Lkey, Jkey, wkey, xkey, fkey, ukey, rkey)

                                L = rnode.attrs['model_size']
                                J = [rnode.attrs['J1_mean'], rnode.attrs['J2_mean'], rnode.attrs['J3_mean']]
                                w = [rnode.attrs['J1_wdth'], rnode.attrs['J2_wdth'], rnode.attrs['J3_wdth']]
                                r = rnode.attrs['J2_span']
                                x = rnode.attrs['xi_Jcls']
                                f = rnode.attrs['f_mixer']
                                u = rnode.attrs['u_layer']
                                if debug:
                                    print("Adding keys")
                                db['keys']['L'].append(Lkey) if not Lkey in db['keys']['L'] else db['keys']['L']
                                db['keys']['J'].append(Jkey) if not Jkey in db['keys']['J'] else db['keys']['J']
                                db['keys']['w'].append(wkey) if not wkey in db['keys']['w'] else db['keys']['w']
                                db['keys']['x'].append(xkey) if not xkey in db['keys']['x'] else db['keys']['x']
                                db['keys']['f'].append(fkey) if not fkey in db['keys']['f'] else db['keys']['f']
                                db['keys']['u'].append(ukey) if not ukey in db['keys']['u'] else db['keys']['u']
                                # r is special. It can either be r < L or r == L
                                # When r < L we have decided to truncate some interactions, in which case that is interesting.
                                # When r == L we are not interested in varying r. So then merging all r is a good idea
                                # For instance we can set
                                if r == L:
                                    db['keys']['r'].append('r_') if not 'r_' in db['keys']['r'] else db['keys']['r']
                                else:
                                    db['keys']['r'].append(rkey) if not rkey in db['keys']['r'] else db['keys']['r']

                                if debug:
                                    print("Adding vals")
                                db['vals']['L'].append(L) if not L in db['vals']['L'] else db['vals']['L']
                                db['vals']['J'].append(J) if not J in db['vals']['J'] else db['vals']['J']
                                db['vals']['w'].append(w) if not w in db['vals']['w'] else db['vals']['w']
                                db['vals']['x'].append(x) if not x in db['vals']['x'] else db['vals']['x']
                                db['vals']['f'].append(f) if not f in db['vals']['f'] else db['vals']['f']
                                db['vals']['u'].append(u) if not u in db['vals']['u'] else db['vals']['u']
                                if r == np.iinfo(np.uint64).max:
                                    db['vals']['r'].append('$L$') if not r in db['vals']['r'] else db['vals']['r']
                                else:
                                    db['vals']['r'].append(r) if not r in db['vals']['r'] else db['vals']['r']

                                for algokey, algopath, algonode in h5py_group_iterator(node=rnode, keypattern=algo_filter, dep=1):
                                    db['keys']['algo'].append(algokey) if not algokey in db['keys']['algo'] else db['keys']['algo']

                                    for modelkey, modelpath, modelnode in h5py_group_iterator(node=algonode, keypattern=model_filter, dep=1):
                                        db['keys']['model'].append(modelkey) if not modelkey in db['keys']['model'] else db['keys']['model']
                                        for metakey, descr in meta.items():
                                            if metakey == 'include' or metakey == 'common':
                                                continue
                                            if not 'dsetname' in descr:
                                                continue
                                            if debug:
                                                print('Looking for dataset {} in {}'.format(descr['dsetname'],
                                                                                            modelnode.name))
                                            for datakey, datapath, datanode in h5py_node_iterator(node=modelnode,
                                                                                                  keypattern=descr[
                                                                                                      'dsetname'],
                                                                                                  dep=4,
                                                                                                  nodeType=h5py.Group):
                                                print("Loading dset database version 2: {}".format(datapath))

                                                num = np.max(datanode['num'][()])
                                                for dname in [datanode.name]:
                                                    db['dsets'][dname] = {}
                                                    db['dsets'][dname]['version'] = db['version']
                                                    db['dsets'][dname]['keys'] = db['keys']
                                                    db['dsets'][dname]['vals'] = {}
                                                    db['dsets'][dname]['vals']['L'] = L
                                                    db['dsets'][dname]['vals']['J'] = J
                                                    db['dsets'][dname]['vals']['w'] = w
                                                    db['dsets'][dname]['vals']['x'] = x
                                                    db['dsets'][dname]['vals']['f'] = f
                                                    db['dsets'][dname]['vals']['u'] = u
                                                    db['dsets'][dname]['vals']['r'] = L if r == np.iinfo(np.uint64).max else r
                                                    db['dsets'][dname]['vals']['num'] = num

                                                    db['dsets'][dname]['node'] = {}
                                                    db['dsets'][dname]['node']['L'] = Lnode
                                                    db['dsets'][dname]['node']['J'] = Jnode
                                                    db['dsets'][dname]['node']['w'] = wnode
                                                    db['dsets'][dname]['node']['x'] = xnode
                                                    db['dsets'][dname]['node']['f'] = fnode
                                                    db['dsets'][dname]['node']['u'] = unode
                                                    db['dsets'][dname]['node']['model'] = modelnode
                                                    db['dsets'][dname]['node']['data'] = datanode
                                                    db['dsets'][dname]['node']['avg'] = datanode['avg']

                                                    db['dsets'][dname]['tex'] = {
                                                        'keys': db['tex'],
                                                        'vals': {
                                                            'L': '{}'.format(L),
                                                            'J': '{}'.format(J),
                                                            'w': '{}'.format(w),
                                                            'x': '{}'.format(x),
                                                            'f': '{}'.format(f),
                                                            'u': '{}'.format(u),
                                                            'r': '{}'.format(r),
                                                            'num': '{}'.format(num),
                                                            'algo': algokey,
                                                        },
                                                        'eqs': {
                                                            'L': '${}{}{}$'.format(db['tex']['L'].strip('$'), '{:}', L),
                                                            'J': '$J:{}$'.format(J),
                                                            'w': '$w:{}$'.format(w),
                                                            'x': '${}{}{:>.4f}$'.format(db['tex']['x'].strip('$'),
                                                                                        '{:}', x),
                                                            'f': '${}{}{:>.4f}$'.format(db['tex']['f'].strip('$'),
                                                                                        '{:}', f),
                                                            'u': '${}{}{:>.0f}$'.format(db['tex']['u'].strip('$'),
                                                                                        '{:}', u),
                                                            'r': '$r = L$' if r == np.iinfo(
                                                                np.uint64).max else '$r = {}$'.format(r),
                                                        },
                                                    }

                                    for statekey, statepath, statenode in h5py_group_iterator(node=algonode, keypattern=state_filter, dep=1):
                                        db['keys']['state'].append(statekey) if not statekey in db['keys']['state'] else db['keys']['state']

                                        for cronokey, cronopath, crononode in h5py_group_iterator(node=statenode, keypattern='cronos', dep=1):
                                            db['keys']['crono'].append(cronokey) if not cronokey in db['keys']['crono'] else db['keys']['crono']

                                            for metakey, descr in meta.items():
                                                if metakey == 'include' or metakey == 'common':
                                                    continue
                                                objname = None
                                                if 'colname' in descr:  # Consider only tables
                                                    objname = descr['colname']
                                                elif 'dsetname' in descr:  # ... or datasets
                                                    objname = descr['dsetname']
                                                else:
                                                    continue
                                                print("Loading time database version 2: {}".format(cronopath))

                                                for datakey, datapath, datanode in h5py_node_iterator(node=crononode,
                                                                                                      keypattern=descr['groupname'],
                                                                                                      dep=6,
                                                                                                      nodeType=h5py.Group):
                                                    if datanode.name in db['dsets']:
                                                        continue
                                                    if debug:
                                                        print(Lkey, Jkey, wkey, xkey, fkey, ukey, algokey, statekey, cronokey, datakey)
                                                    num = 0
                                                    tsim = 0.0
                                                    bavg_key = 0.0
                                                    if 'colname' in descr:  # We have a table!
                                                        num = np.max(datanode['max']['num'][()])
                                                        tsim = np.max(datanode['avg']['algorithm_time'][()]) / 60
                                                        bavg_key, bavg, bmax_key, bmax = get_bond_info(statenode, datanode)
                                                    else:
                                                        mmntnode = datanode.parent['measurements']
                                                        num = np.shape(datanode['data'])[1]
                                                        tsim = np.max(mmntnode['avg']['algorithm_time'][()]) / 60
                                                        bavg_key, bavg, bmax_key, bmax = get_bond_info(statenode, mmntnode)

                                                    if debug:
                                                        print("Adding node data")
                                                    for dname in [datanode.name]:
                                                        db['dsets'][dname] = {}
                                                        db['dsets'][dname]['version'] = db['version']
                                                        db['dsets'][dname]['keys'] = db['keys']
                                                        db['dsets'][dname]['vals'] = {}
                                                        db['dsets'][dname]['vals']['L'] = L
                                                        db['dsets'][dname]['vals']['J'] = J
                                                        db['dsets'][dname]['vals']['w'] = w
                                                        db['dsets'][dname]['vals']['x'] = x
                                                        db['dsets'][dname]['vals']['f'] = f
                                                        db['dsets'][dname]['vals']['u'] = u
                                                        db['dsets'][dname]['vals']['r'] = L if r == np.iinfo(np.uint64).max else r
                                                        db['dsets'][dname]['vals']['num'] = num
                                                        db['dsets'][dname]['vals']['tsim'] = tsim
                                                        db['dsets'][dname]['vals']['bavg'] = bavg
                                                        db['dsets'][dname]['vals']['bmax'] = bmax

                                                        db['dsets'][dname]['node'] = {}
                                                        db['dsets'][dname]['node']['L'] = Lnode
                                                        db['dsets'][dname]['node']['J'] = Jnode
                                                        db['dsets'][dname]['node']['w'] = wnode
                                                        db['dsets'][dname]['node']['x'] = xnode
                                                        db['dsets'][dname]['node']['f'] = fnode
                                                        db['dsets'][dname]['node']['u'] = unode
                                                        db['dsets'][dname]['node']['state'] = statenode
                                                        db['dsets'][dname]['node']['table'] = statenode.get('tables')
                                                        db['dsets'][dname]['node']['crono'] = crononode
                                                        db['dsets'][dname]['node']['data'] = datanode
                                                        db['dsets'][dname]['node']['avg'] = datanode['avg']

                                                        db['dsets'][dname]['tex'] = {
                                                            'keys': db['tex'],
                                                            'vals': {
                                                                'L': '{}'.format(L),
                                                                'J': '{}'.format(J),
                                                                'w': '{}'.format(w),
                                                                'x': '{}'.format(x),
                                                                'f': '{}'.format(f),
                                                                'u': '{}'.format(u),
                                                                'r': '{}'.format(r),
                                                                'num': '{}'.format(num),
                                                                'tsim': '{:>.1f}m'.format(tsim),
                                                                'bavg': '{}'.format(bavg),
                                                                'bmax': '{}'.format(bmax),
                                                                'algo': algokey,
                                                                'state': statekey,
                                                                'crono': cronokey,
                                                                'data': datakey,
                                                            },
                                                            'eqs': {
                                                                'L': '${}{}{}$'.format(db['tex']['L'].strip('$'), '{:}', L),
                                                                'J': '$J:{}$'.format(J),
                                                                'w': '$w:{}$'.format(w),
                                                                'x': '${}{}{:>.4f}$'.format(db['tex']['x'].strip('$'), '{:}', x),
                                                                'f': '${}{}{:>.4f}$'.format(db['tex']['f'].strip('$'), '{:}', f),
                                                                'u': '${}{}{:>.0f}$'.format(db['tex']['u'].strip('$'), '{:}', u),
                                                                'r': '$r = L$' if r == np.iinfo(np.uint64).max else '$r = {}$'.format(r),
                                                            },
                                                        }
                                                    if debug:
                                                        print("Finished gather")

    # Sort the keys so that we can iterate through them in order
    # sortL = np.argsort(np.asarray(db['vals']['L']))
    # sortb = np.argsort(np.asarray(db['vals']['b']))
    # sortf = np.argsort(np.asarray(db['vals']['f']))
    # sortu = np.argsort(np.asarray(db['vals']['u']))
    #
    # db['vals']['L']   = np.asarray(db['vals']['L'])[sortL]
    # db['keys']['L']   = np.asarray(db['keys']['L'])[sortL]
    # db['vals']['b']  = np.asarray(db['vals']['b'])[sortb]
    # db['keys']['b']  = np.asarray(db['keys']['b'])[sortb]
    # db['vals']['f']  = np.asarray(db['vals']['f'])[sortf]
    # db['keys']['f']  = np.asarray(db['keys']['f'])[sortf]
    # db['vals']['u']  = np.asarray(db['vals']['u'])[sortu]
    # db['keys']['u']  = np.asarray(db['keys']['u'])[sortu]
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
