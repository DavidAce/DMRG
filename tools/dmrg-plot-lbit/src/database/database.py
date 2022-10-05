from src.io.h5ops import *
from src.plotting.filter import *
import itertools


def get_bond_info(statenode, datanode):
    bavg_key = None
    bmax_key = None
    bavg = None
    bmax = None
    try:
        bavg_key = next(k for k in ['bond_mid', 'bond_dimension_midchain'] if k in datanode['avg'].dtype.fields.keys())
        bavg = np.max(datanode['avg'][bavg_key][()])
    except:
        pass
    try:
        bmax_key = next(k for k in ['bond_lim', 'chi_lim'] if k in statenode['tables']['status']['max'].dtype.fields.keys())
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
                                                    print(Lkey, Jkey, wkey, bkey, fkey, ukey, algokey, statekey, pointkey, datakey)

                                                L = rnode.attrs['model_size']
                                                J = [rnode.attrs['J1_mean'], rnode.attrs['J2_mean'], rnode.attrs['J3_mean']]
                                                w = [rnode.attrs['J1_wdth'], rnode.attrs['J2_wdth'], rnode.attrs['J3_wdth']]
                                                b = rnode.attrs['J2_base']
                                                f = rnode.attrs['f_mixer']
                                                u = rnode.attrs['u_layer']
                                                r = rnode.attrs['J2_span']

                                                if debug:
                                                    print("Adding keys")
                                                db['keys']['L'].append(Lkey) if not Lkey in db['keys']['L'] else db['keys']['L']
                                                db['keys']['J'].append(Jkey) if not Jkey in db['keys']['J'] else db['keys']['J']
                                                db['keys']['w'].append(wkey) if not wkey in db['keys']['w'] else db['keys']['w']
                                                db['keys']['b'].append(bkey) if not bkey in db['keys']['b'] else db['keys']['b']
                                                db['keys']['f'].append(fkey) if not fkey in db['keys']['f'] else db['keys']['f']
                                                db['keys']['u'].append(ukey) if not ukey in db['keys']['u'] else db['keys']['u']
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


def load_time_database2(h5_src, meta, algo_filter=None, model_filter=None, state_filter=None, debug=False):
    # Gather data

    db = {
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
            'L': '$L$', 'J': '$J$', 'w': '$w$', 'x': '$\\xi$', 'f': '$f$', 'u': '$u$', 'r': '$r$',
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
                                x = rnode.attrs['J2_xcls']
                                f = rnode.attrs['f_mixer']
                                u = rnode.attrs['u_layer']
                                r = rnode.attrs['J2_span']
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
                                            for datakey, datapath, datanode in h5py_node_iterator(node=modelnode,
                                                                                                  keypattern=descr['dsetname'],
                                                                                                  dep=2,
                                                                                                  nodeType=h5py.Group):
                                                print("Loading dset database version 2: {}[{}]".format(modelkey,
                                                                                                       descr['dsetname']))

                                                num = np.max(datanode['num'][()])
                                                for dname in [datanode.name]:
                                                    db['dsets'][dname] = {}
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
                                                            'x': '${}{}{:>.4f}$'.format(db['tex']['x'].strip('$'), '{:}', x),
                                                            'f': '${}{}{:>.4f}$'.format(db['tex']['f'].strip('$'), '{:}', f),
                                                            'u': '${}{}{:>.0f}$'.format(db['tex']['u'].strip('$'), '{:}', u),
                                                            'r': '$r = L$' if r == np.iinfo(np.uint64).max else '$r = {}$'.format(
                                                                r),
                                                        },
                                                    }

                                    for statekey, statepath, statenode in h5py_group_iterator(node=algonode, keypattern=state_filter, dep=1):
                                        db['keys']['state'].append(statekey) if not statekey in db['keys']['state'] else db['keys']['state']

                                        for cronokey, cronopath, crononode in h5py_group_iterator(node=statenode, keypattern='cronos', dep=1):
                                            db['keys']['crono'].append(cronokey) if not cronokey in db['keys']['crono'] else db['keys']['crono']

                                            for metakey, descr in meta.items():
                                                if metakey == 'include' or metakey == 'common':
                                                    continue
                                                if not 'colname' in descr:  # Consider only tables
                                                    continue

                                                print("Loading time database version 2: {}[{}]".format(descr['groupname'], descr['colname']))

                                                for datakey, datapath, datanode in h5py_node_iterator(node=crononode,
                                                                                                      keypattern=descr['groupname'],
                                                                                                      dep=6,
                                                                                                      nodeType=h5py.Group):
                                                    if datanode.name in db['dsets']:
                                                        continue
                                                    if debug:
                                                        print(Lkey, Jkey, wkey, xkey, fkey, ukey, algokey, statekey, cronokey, datakey)

                                                    num = np.max(datanode['max']['num'][()])
                                                    tsim = np.max(datanode['avg']['algorithm_time'][()]) / 60
                                                    bavg_key, bavg, bmax_key, bmax = get_bond_info(statenode, datanode)

                                                    if debug:
                                                        print("Adding node data")
                                                    for dname in [datanode.name]:
                                                        db['dsets'][dname] = {}
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
                                                        db['dsets'][dname]['node']['table'] = statenode['tables']
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

    # Sort the keys so that we can iterate through them in order
    sortL = np.argsort(np.asarray(db['vals']['L']))
    sortl = np.argsort(np.asarray(db['vals']['lambda']))
    sortd = np.argsort(np.asarray(db['vals']['delta']))

    db['vals']['L'] = np.asarray(db['vals']['L'])[sortL]
    db['keys']['L'] = np.asarray(db['keys']['L'])[sortL]
    db['vals']['lambda'] = np.asarray(db['vals']['lambda'])[sortl]
    db['keys']['lambda'] = np.asarray(db['keys']['lambda'])[sortl]
    db['vals']['delta'] = np.asarray(db['vals']['delta'])[sortd]
    db['keys']['delta'] = np.asarray(db['keys']['delta'])[sortd]

    return db
