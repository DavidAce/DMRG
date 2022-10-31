from src.io.h5ops import *
from src.plotting.filter import *


def load_database(h5_src, meta, debug=False):
    db = {
        'max': {'all': None, 'avg': None, 'mid': None, 'mvg': None},
        'min': {'all': None, 'avg': None, 'mid': None, 'mvg': None},
        'keys': {
            'L': [],
            'l': [],
            'd': [],
            'algo': [],
            'state': [],
            'point': [],
            'data': [],
        },
        'vals': {
            'L': [],
            'l': [],
            'd': [],
            'algo': [], 'state': [], 'point': [], 'data': [],
        },
        'tex': {
            'L': '$L$',
            'l': '$\lambda$',
            'd': '$\Delta$',
            'num': '$n$',
            'time': '$\\bar t$',
            'algo': '$\mathrm{Alg}$',
            'state': '$\psi_n$',
            'point': '$-$',
            'data': '$.$'

        },
        'texdef': {
            'd': '$\Delta = \log \\bar J - \log \\bar h$',
        },
        'dsets': {},
    }

    digit = lambda x: int(''.join(list(filter(str.isdigit, x) or 0)))
    for Lkey, Lnode in sorted(h5_src.items(), key=natural_items):
        for lkey, lnode in sorted(Lnode.items(), key=natural_items):
            for dkey, dnode in sorted(lnode.items(), key=natural_items):
                for algokey, algonode in sorted(dnode.items(), key=natural_items):
                    for statekey, statenode in sorted(algonode.items(), key=natural_items):
                        for pointkey, pointnode in sorted(statenode.items(), key=natural_items):

                            if 'common' in meta and 'include' in meta['common']:
                                if 'algo' in meta['common']['include'] and not algokey in meta['common']['include']['algo']:
                                    continue
                                if 'state' in meta['common']['include'] and not statekey in meta['common']['include']['state']:
                                    continue
                                if 'point' in meta['common']['include'] and not pointkey in meta['common']['include']['point']:
                                    continue
                                if 'L' in meta['common']['include'] and not Lkey in meta['common']['include']['L']:
                                    continue
                                if 'l' in meta['common']['include'] and not lkey in meta['common']['include']['l']:
                                    continue
                                if 'g' in meta['common']['include'] and not lkey in meta['common']['include']['g']:
                                    continue
                                if 'd' in meta['common']['include'] and not dkey in meta['common']['include']['d']:
                                    continue

                            print("Loading databases in {}".format(pointnode.name))
                            for metakey, descr in meta.items():
                                if metakey == 'include' or metakey == 'common':
                                    continue
                                if debug:
                                    print("Starting gather | L [{}] | l [{}] | d [{}] | algo [{}] | state [{}] | point [{}]".format(Lkey, lkey, dkey, algokey,
                                                                                                                                    statekey, pointkey))
                                    print("Found group matching [{}]: [{}]".format(descr['load'],
                                                                                   h5py_node_finder(node=pointnode,
                                                                                                    keypattern=descr['load'],
                                                                                                    dep=3, nodeType=h5py.Group)))
                                for datakey, datapath, datanode in h5py_node_finder(node=pointnode,
                                                                                    keypattern=descr['load'],
                                                                                    dep=3, nodeType=h5py.Group):
                                    if datapath in db['dsets']:
                                        continue
                                    print(" -- Loading node: {}".format(datapath))

                                    tavg = None
                                    if not 'ed' in algokey:
                                        dsetnode = datanode.parent  # Go back one group
                                        tavg = np.max(dsetnode['measurements']['avg']['algorithm_time'][()])

                                    L = dnode.attrs['model_size']
                                    l = dnode.attrs['lambda'] if 'lambda' in dnode.attrs else dnode.attrs['g']
                                    d = dnode.attrs['delta']
                                    n = np.min(datanode['num'][()])
                                    midx = int(L / 2)
                                    mlbl = 'L_{}'.format(midx)

                                    if debug:
                                        print("Adding keys")
                                    db['keys']['L'].append(Lkey) if not Lkey in db['keys']['L'] else db['keys']['L']
                                    db['keys']['l'].append(lkey) if not lkey in db['keys']['l'] else db['keys']['l']
                                    db['keys']['d'].append(dkey) if not dkey in db['keys']['d'] else db['keys']['d']
                                    db['keys']['algo'].append(algokey) if not algokey in db['keys']['algo'] else db['keys']['algo']
                                    db['keys']['state'].append(statekey) if not statekey in db['keys']['state'] else db['keys']['state']
                                    db['keys']['point'].append(pointkey) if not pointkey in db['keys']['point'] else db['keys']['point']
                                    db['keys']['data'].append(datakey) if not datakey in db['keys']['data'] else db['keys']['data']

                                    if debug:
                                        print("Adding vals")
                                    db['vals']['L'].append(L) if not L in db['vals']['L'] else db['vals']['L']
                                    db['vals']['l'].append(l) if not l in db['vals']['l'] else db['vals']['l']
                                    db['vals']['d'].append(d) if not d in db['vals']['d'] else db['vals']['d']

                                    if debug:
                                        print("Adding tex")
                                    if not 'texeq' in db:
                                        db['texeq'] = {}
                                    if not 'tex' in db:
                                        db['texeq'] = {}
                                    db['texeq'][Lkey] = '${}{}{}$'.format(db['tex']['L'].strip('$'), '{:}', L)
                                    db['texeq'][lkey] = '${}{}{:.4f}$'.format(db['tex']['l'].strip('$'), '{:}', l)
                                    db['texeq'][dkey] = '${}{}{:.4f}$'.format(db['tex']['d'].strip('$'), '{:}', d)
                                    db['tex'][Lkey] = "${}$".format(L)
                                    db['tex'][lkey] = "${:.4f}$".format(l)
                                    db['tex'][dkey] = "${:.4f}$".format(d)

                                    if debug:
                                        print("Adding node data")
                                    for dname in [datanode.name]:
                                        db['dsets'][dname] = {}
                                        db['dsets'][dname]['vals'] = {}
                                        db['dsets'][dname]['path'] = {}
                                        db['dsets'][dname]['keys'] = db['keys']
                                        db['dsets'][dname]['num'] = n
                                        db['dsets'][dname]['vals']['L'] = L
                                        db['dsets'][dname]['vals']['l'] = l
                                        db['dsets'][dname]['vals']['d'] = d
                                        db['dsets'][dname]['vals']['midx'] = midx
                                        db['dsets'][dname]['vals']['mlbl'] = mlbl
                                        db['dsets'][dname]['vals']['time'] = tavg / 60

                                        db['dsets'][dname]['path']['L'] = Lnode.name
                                        db['dsets'][dname]['path']['l'] = lnode.name
                                        db['dsets'][dname]['path']['d'] = dnode.name
                                        db['dsets'][dname]['path']['state'] = statenode.name
                                        db['dsets'][dname]['path']['point'] = pointnode.name
                                        db['dsets'][dname]['path']['data'] = datanode.name

                                        if 'tables' in statenode:
                                            db['dsets'][dname]['path']['table'] = statenode['tables'].name
                                            if 'measurements' in statenode['tables']:
                                                db['dsets'][dname]['path']['mmnt'] = statenode['tables']['measurements'].name

                                        db['dsets'][dname]['tex'] = {
                                            'keys': db['tex'],
                                            'vals': {
                                                'L': '{}'.format(L),
                                                'l': '{:.4f}'.format(l),
                                                'd': '{:.4f}'.format(d),
                                                'num': '{}'.format(n),
                                                'algo': 'ED' if 'ed' in algokey else algokey,
                                                'state': '-' if 'ed' in algokey else datakey,
                                                'point': '-' if 'ed' in algokey else ('fes-dec' if pointkey == 'fes' else 'fes-inc'),
                                                'data': '-' if 'ed' in algokey else datakey,
                                                'time': '-' if 'ed' in algokey else '{:>.1f}m'.format(tavg / 60),

                                            },
                                            'eqs': {
                                                'L': '${}{}{}$'.format(db['tex']['L'].strip('$'), '{:}', L),
                                                'l': '${}{}{:.4f}$'.format(db['tex']['l'].strip('$'), '{:}', l),
                                                'd': '${}{}{:.4f}$'.format(db['tex']['d'].strip('$'), '{:}', d),
                                            },
                                        }

                                        if debug:
                                            print("Adding styles")

                                        db['dsets'][datanode.name]['style'] = {}
                                        if "states" in statekey:
                                            print("Found ed: ", statekey)
                                            db['dsets'][dname]['style']['lwidth'] = 3.0
                                            db['dsets'][dname]['style']['lalpha'] = 1.0
                                            db['dsets'][dname]['style']['mstyle'] = None
                                            db['dsets'][dname]['style']['lstyle'] = 'solid'
                                        else:
                                            db['dsets'][dname]['style']['lwidth'] = 1.4
                                            db['dsets'][dname]['style']['lalpha'] = 1.0
                                            db['dsets'][dname]['style']['mstyle'] = '.'
                                            db['dsets'][dname]['style']['lstyle'] = 'dotted'
                                    if debug:
                                        print("Finished gather")

    # Sort the keys so that we can iterate through them in order

    for key in ['L', 'l', 'd', 'b']:
        if key in db['keys'] and key in db['vals'] and len(db['vals'][key]) > 1:
            sort = np.argsort(db['vals'][key])
            db['vals'][key] = list(np.asarray(db['vals'][key])[sort])
            if len(db['keys'][key]) == len(db['vals'][key]):
                db['keys'][key] = list(np.asarray(db['keys'][key])[sort])

    return db
