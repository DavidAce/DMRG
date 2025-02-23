from itertools import product
from pathlib import Path
import h5py
import matplotlib.transforms as transforms
import matplotlib.patheffects as pe
import scipy.stats
import seaborn as sns
from scipy.optimize import curve_fit
from scipy.stats import norm
import logging
from .tools import *

logger = logging.getLogger('plot-divg')

# def floglog(x, a, b,c):
#     return a + b*np.log(np.log(x + c))
def floglog_v2(x, a, b, c):
    with np.errstate(invalid='ignore'):
        return a + b * np.log(np.log(x - c))


def flinear(x, a, b):
    with np.errstate(invalid='ignore'):
        return a + b * x


def fpower(x, a, b):
    with np.errstate(invalid='ignore'):
        return a * x ** b


def plot_divg_v3_fig_sub_line(db, meta, figspec, subspec, linspec, algo_filter=None, state_filter=None,
                              point_filter=None, figs=None, palette_name=None,dbidx=0,dbnum=1):
    path_effects = [pe.SimpleLineShadow(offset=(0.35, -0.35), alpha=0.1), pe.Normal()]
    if 'mplstyle' in meta:
        plt.style.use(meta['mplstyle'])
        if 'slack' in meta['mplstyle']:
            path_effects = [pe.SimpleLineShadow(offset=(0.5, -0.5), alpha=0.3), pe.Normal()]

    legend_col_keys = meta.get('legendcols')
    # if legendcols := meta.get('legendcols'):
    #     for col in legendcols:
    #         speclist = [l.split(':')[0] for l in figspec + subspec + linspec]
    #         if not col in speclist and not col in legend_col_keys:
    #             legend_col_keys.append(col)

    figprod = list(product(*get_vals(db=db, keyfmt=figspec, filter=meta.get('filter'))))  # All combinations of figspecs values
    subprod = list(product(*get_vals(db=db, keyfmt=subspec, filter=meta.get('filter'))))  # All combinations of subspecs values
    linprod = list(product(*get_vals(db=db, keyfmt=linspec, filter=meta.get('filter'))))  # All combinations of linspecs values
    numfigs = len(figprod)
    numsubs = len(subprod)

    numsubs *=dbnum
    subprod *=dbnum
    if figs is None:
        figs = [get_fig_meta(numsubs, meta=meta) for _ in range(numfigs)]

    for fidx, (figvals, f) in enumerate(zip(figprod, figs)):
        logger.debug('- plotting figs: {}'.format(figvals))
        idx_palette = fidx % len(palette_name)
        dbval = None
        luitz = []
        for idx, (subvals, ax) in enumerate(zip(subprod, f['ax'])):
            idx_palette = idx % len(palette_name) if numsubs > 1 else idx_palette # Iterate the palette if there are more than one subplots
            if dbidx != int(idx/dbnum):
                continue
            popt = None
            pcov = None
            dbval = None
            legendrow = None
            logger.debug('-- plotting subs: {}'.format(subvals))
            palette, lstyles = get_colored_lstyles(db, linspec, palette_name, filter=None, idx=idx_palette)
            for linvals, color, lstyle in zip(linprod, palette, lstyles):
                logger.debug('--- plotting lins: {}'.format(linvals))
                datanodes = match_datanodes(db=db, meta=meta, specs=figspec + subspec + linspec,
                                            vals=figvals + subvals + linvals)

                if len(datanodes) != 1:
                    logger.warning(f"match: \n"
                                   f"\tspec:{[figspec + subspec + linspec]}\n"
                                   f"\tvals:{[figvals + subvals + linvals]}")
                    logger.warning(f"found {len(datanodes)} datanodes: {datanodes=}")
                    # continue
                for datanode in datanodes:
                    mmntnode = datanode.parent['measurements']
                    dbval = db['dsets'][datanode.name]
                    ndata = datanode['avg']['num'][()]
                    if np.min(ndata) < 10:
                        logger.warning("ndata < 10: in {}\n\t continuing".format(datanode))
                        continue
                    t = mmntnode['avg']['physical_time'][()].astype(float)
                    y = datanode[meta['dsetname']][()]
                    n = datanode['avg']['num'][()][0]
                    tp = get_timepoints(t, dbval)
                    idx_num = tp.idx_num_saturated
                    idx_ent = tp.idx_ent_saturated
                    idx_sat = idx_num if 'number' in meta['dsetname'] else idx_ent
                    print(f'L={dbval["vals"]["L"]:2} f={dbval["vals"]["f"]:3} n={n} t◀={t[idx_num]:.3e} t●{t[idx_ent]:.3e} sat: [num {len(t[idx_num:])} ent {len(t[idx_ent:])}]')
                    if meta.get('normpage'):
                        y /= page_entropy(dbval['L'])
                    if normalize := meta.get('normalize'):
                        y /= normalize
                    # if t[-1] == t[idx_sat]:
                    #     print('Time window is not long enough: saturation index = {} / {}', idx_sat, len(t))
                    #     continue
                    legendrow = get_legend_row(db=db, datanode=datanode, legend_col_keys=legend_col_keys)

                    if meta.get('plotpinfty'):
                        # ytavg = datanode.parent.parent['pinfty_number_entropies']['data'][()]
                        # Take the shannon-entropy of infinite-time averaged number probabilities
                        if 'pinfty_number_entropies' in  datanode.parent.parent:
                            ytavg = datanode.parent.parent['pinfty_number_entropies']['data'][()]
                        elif 'number_probabilities' in datanode.parent.parent:
                            print('Warning: Could not find pinfty_number entropies. Calculating from scratch!')
                            p = datanode.parent.parent['number_probabilities'][()]
                            ptavg = np.mean(p[:, idx_sat:, :], axis=1)
                            ptavg = np.ma.masked_invalid(np.ma.masked_equal(np.abs(ptavg), 0))
                            with np.errstate(divide='ignore'):
                                ytavg = -np.sum(ptavg * np.log(ptavg), axis=0)
                        else:
                            raise LookupError('Could not find data to generate pinfty entropies')
                    elif meta.get('plotPosX_neel0'):
                        L = dbval['vals']['L']
                        var0 = datanode.parent.parent['nth_particle_position/pos_expvalue_neel0'][int(L // 4) - 1, :,:]  # 1010|1010 index L//4 = 2 is nearest the midchain boundary
                        # var1 = datanode.parent.parent['nth_particle_position/pos_variance_neel1'][int(L // 4) - 0, :,:]  # 1010|1010 index L//4 = 2 is nearest the midchain boundary
                        # y = np.concatenate([var0, var1], axis=1)
                        idx_sat = idx_num # Saturates later than num?
                        ytavg = np.mean(var0[idx_sat:, :], axis=0)
                    elif meta.get('plotPosX_neel1'):
                        L = dbval['vals']['L']
                        # var0 = datanode.parent.parent['nth_particle_position/pos_expvalue_neel0'][int(L // 4) - 1, :,:]  # 1010|1010 index L//4 = 2 is nearest the midchain boundary
                        var1 = datanode.parent.parent['nth_particle_position/pos_expvalue_neel1'][int(L // 4) - 0, :,:]  # 1010|1010 index L//4 = 2 is nearest the midchain boundary
                        # y = np.concatenate([var0, var1], axis=1)
                        idx_sat = idx_num # Saturates later than num?
                        ytavg = np.mean(var1[idx_sat:, :], axis=0)
                    elif meta.get('plotVarX'):
                        L = dbval['vals']['L']
                        idx_sat = idx_ent
                        posnode = datanode.parent.parent['nth_particle_position']
                        N = int(L // 4) # Number of particles
                        r2 = range(N - 1, N + 1)  # Two central particles
                        y0 = np.atleast_2d(np.mean(posnode['pos_variance_neel0'][r2, idx_sat:, :],axis=1))  # 1010|1010 index L//4 = 2 is nearest the midchain boundary
                        y1 = np.atleast_2d(np.mean(posnode['pos_variance_neel1'][r2, idx_sat:, :],axis=1))  # 1010|1010 index L//4 = 2 is nearest the midchain boundary
                        ytavg = np.ravel(np.concatenate([y0,y1], axis=1)) # Do not average over particles
                        print(f'shape ytavg: {np.shape(ytavg)}')

                        # L = dbval['vals']['L']
                        # var0 = datanode.parent.parent['nth_particle_position/pos_variance_neel0'][int(L // 4) - 1, :,:]  # 1010|1010 index L//4 = 2 is nearest the midchain boundary
                        # var1 = datanode.parent.parent['nth_particle_position/pos_variance_neel1'][int(L // 4) - 0, :,:]  # 1010|1010 index L//4 = 2 is nearest the midchain boundary
                        # var_davg_evn = datanode.parent.parent['pos_variance_davg_evn'][int(L // 4), :]   # 1010|1010 index L//4 = 2 is nearest the midchain boundary
                        # var_davg_odd = datanode.parent.parent['pos_variance_davg_odd'][int(L // 4)-1, :] # 0101|0101 index L//4 -1 = is nearest the midchain boundary
                        # var_dste_evn = np.std(var_evn, axis=2)/np.sqrt(np.shape(var_evn)[2])
                        # var_dste_odd = np.std(var_odd, axis=2)/np.sqrt(np.shape(var_odd)[2])
                        # y = np.concatenate([var0, var1], axis=1)
                        # idx_sat = idx_num # Saturates later than num?
                        # ytavg = np.mean(var0[idx_sat:, :], axis=0)
                    # elif meta.get('plotVarX_neel1'):
                    #     L = dbval['vals']['L']
                    #     # var0 = datanode.parent.parent['nth_particle_position/pos_variance_neel0'][int(L // 4) - 1, :,:]  # 1010|1010 index L//4 = 2 is nearest the midchain boundary
                    #     var1 = datanode.parent.parent['nth_particle_position/pos_variance_neel1'][int(L // 4) - 0, :,:]  # 1010|1010 index L//4 = 2 is nearest the midchain boundary
                    #     # var_davg_evn = datanode.parent.parent['pos_variance_davg_evn'][int(L // 4), :]   # 1010|1010 index L//4 = 2 is nearest the midchain boundary
                    #     # var_davg_odd = datanode.parent.parent['pos_variance_davg_odd'][int(L // 4)-1, :] # 0101|0101 index L//4 -1 = is nearest the midchain boundary
                    #     # var_dste_evn = np.std(var_evn, axis=2)/np.sqrt(np.shape(var_evn)[2])
                    #     # var_dste_odd = np.std(var_odd, axis=2)/np.sqrt(np.shape(var_odd)[2])
                    #     # y = np.concatenate([var0, var1], axis=1)
                    #     idx_sat = idx_num # Saturates later than num?
                    #     ytavg = np.mean(var1[idx_sat:, :], axis=0)
                    else:
                        # Calculate the infinite time average (1/T) integral_0^T y(t) dt in the saturated interval
                        ytavg = np.mean(y[idx_sat:, :], axis=0)
                    hist, edges = np.histogram(ytavg, bins=meta['bins'],range=meta['range'], density=True)
                    bincentres = [(edges[j] + edges[j + 1]) / 2. for j in range(len(edges) - 1)]
                    # line, = ax.plot(bincentres, hist,label=None,
                    #                 color=color, path_effects=path_effects)
                    line, = ax.step(x=bincentres, y=hist, where='mid', label=None,
                                    # linewidth=0.85,
                                    color=color, path_effects=path_effects)
                    # ax.set_yscale('log')
                    # mu, std = norm.fit(ytavg)
                    # res = scipy.stats.normaltest(ytavg)
                    # print(res)
                    # ax.plot(bincentres, norm.pdf(bincentres, mu,std),label=None,
                    #
                    #                 linestyle='--',
                    #                 color=color, path_effects=path_effects)
                    if meta.get('markL/2'):
                        L = dbval['vals']['L']
                        ax.axvline(x=L // 2 - 0.5, color='black', alpha=0.85)
                        trans = transforms.blended_transform_factory(ax.transData, ax.transAxes)
                        ax.text(L // 2 - 0.5, 0.90, 'halfchain', fontsize='small', color='black', alpha=0.75, ha='right',
                                va='center', rotation='vertical', transform=trans)
                    icol = 0
                    for (col, key) in zip(legendrow, legend_col_keys):
                        key, fmt = key.split(':') if ':' in key else [key, '']
                        f['legends'][idx][icol]['handle'].append(line)
                        f['legends'][idx][icol]['title'] = db['tex'][key]
                        f['legends'][idx][icol]['label'].append(col)
                        icol += 1
                    # f['legends'][idx][icol]['handle'].append(line)
                    # f['legends'][idx][icol]['title'] = '$w_i$'
                    # f['legends'][idx][icol]['label'].append('off' if f['figcount'] == 0 else 'on')
            if ( meta.get('plotluitz')
                    and f['figcount'] == 0
                    and 'number' in meta['dsetname']
                    and legendrow is not None
                    and not meta.get('plotPosX_neel0')
                    and not meta.get('plotPosX_neel1')
                    and not meta.get('plotVarX')):
                # Plot Luitz's data
                # Let's only do this for the largest system size (16)
                # and only after the last plot, so that this legend entry is last
                lenval = 14 # Picks out the largest system size in Luitz's data
                lenkey = f"L{lenval}"
                with h5py.File('external/raw_EE_NE_CE_distributions_random_XXX_chain.h5', 'r') as h5ext:
                    graypalette = sns.color_palette('Greys', n_colors=5)[1:-1] # Skip edge colors (too bright and too dark)
                    for W, color in zip([3.0, 6.0, 10.0], graypalette):
                    # for W,color in zip([10.0, 6.0, 3.0], graypalette):
                        lwpath = f'{lenkey}/W{W}'
                        if h5ext.get(lwpath) is None:
                            continue
                        if (idx,lwpath) in luitz: # Already plotted on this ax
                            continue
                        luitz.append((idx, lwpath))
                        hist = h5ext[lwpath]['hist[NE1][100]'][()]
                        edges = h5ext[lwpath]['binedges[NE1][100]'][()]
                        bincentres = [(edges[j] + edges[j + 1]) / 2. for j in range(len(edges) - 1)]
                        # line_ext, = ax.plot(bincentres,hist, label=None, color=color, alpha=1.0,
                        #                     path_effects=path_effects,
                                            # zorder=0 )
                        line_ext, = ax.step(x=bincentres, y=hist, where='mid', label=None, color=color, alpha=1.0,
                                            path_effects=path_effects,
                                            linewidth=0.65,
                                            zorder=0)
                        if idx == 0: # Only on the leftmost subplot
                            # key, fmt = key.split(':') if ':' in key else [key, '']
                            f['legends2'][idx][0]['handle'].append(line_ext)
                            f['legends2'][idx][0]['label'] .append(f'{lenval}')
                            if f['legends2'][idx][0]['title'] is None:
                                # f['legends'][idx][0]['header'] = 'PRB\n102.100202'
                                # f['legends'][idx][0]['header'] = 'Luitz et al.'
                                f['legends2'][idx][0]['title'] = '$L$'

                            f['legends2'][idx][1]['handle'].append(line_ext)
                            f['legends2'][idx][1]['label'].append(f'{W}')
                            if f['legends2'][idx][1]['title'] is None:
                                f['legends2'][idx][1]['title'] = '$W$'


                            # f['legends'][idx][0]['label'].append('\makebox[3ex][l]{' + f'$L$=${lenval}$,$W$=${W}$' + '}')
                            for icol, (col, key) in enumerate(zip(legendrow[2:], legend_col_keys[2:])):
                                key, fmt = key.split(':') if ':' in key else [key, '']
                                f['legends2'][idx][icol]['handle'].append(line_ext)
                                f['legends2'][idx][icol]['title'] = db['tex'][key]
                                f['legends2'][idx][icol]['label'].append('')
            if 'entanglement' in meta['dsetname'] and legendrow is not None:
                # Plot Luitz's data
                # Let's only do this for the largest system size (16)
                # and only after the last plot, so that this legend entry is last
                lenval = 14 # Picks out the largest system size in Luitz's data
                lenkey = f"L{lenval}"
                with h5py.File('external/raw_EE_NE_CE_distributions_random_XXX_chain.h5', 'r') as h5ext:
                    graypalette = sns.color_palette('Greys', n_colors=5)[1:-1] # Skip edge colors (too bright and too dark)
                    for W, color in zip([3.0, 6.0, 10.0], graypalette):
                    # for W,color in zip([10.0, 6.0, 3.0], graypalette):
                        lwpath = f'{lenkey}/W{W}'
                        if h5ext.get(lwpath) is None:
                            continue
                        if (idx,lwpath) in luitz: # Already plotted on this ax
                            continue
                        luitz.append((idx, lwpath))
                        hist = h5ext[lwpath]['hist[EE1][100]'][()]
                        edges = h5ext[lwpath]['binedges[EE1][100]'][()]
                        bincentres = [(edges[j] + edges[j + 1]) / 2. for j in range(len(edges) - 1)]
                        # line_ext, = ax.plot(bincentres,hist, label=None, color=color, alpha=1.0,
                        #                     path_effects=path_effects,
                                            # zorder=0 )
                        line_ext, = ax.step(x=bincentres, y=hist, where='mid', label=None, color=color, alpha=1.0,
                                            path_effects=path_effects,
                                            linewidth=0.85,
                                            zorder=0)
                        if idx == 0: # Only on the leftmost subplot
                            # key, fmt = key.split(':') if ':' in key else [key, '']
                            f['legends2'][idx][0]['handle'].append(line_ext)
                            f['legends2'][idx][0]['label'] .append(f'{lenval}')
                            if f['legends2'][idx][0]['title'] is None:
                                # f['legends'][idx][0]['header'] = 'PRB\n102.100202'
                                # f['legends'][idx][0]['header'] = 'Luitz et al.'
                                f['legends2'][idx][0]['title'] = '$L$'

                            f['legends2'][idx][1]['handle'].append(line_ext)
                            f['legends2'][idx][1]['label'].append(f'{W}')
                            if f['legends2'][idx][1]['title'] is None:
                                f['legends2'][idx][1]['title'] = '$W$'


                            # f['legends'][idx][0]['label'].append('\makebox[3ex][l]{' + f'$L$=${lenval}$,$W$=${W}$' + '}')
                            for icol, (col, key) in enumerate(zip(legendrow[2:], legend_col_keys[2:])):
                                key, fmt = key.split(':') if ':' in key else [key, '']
                                f['legends2'][idx][icol]['handle'].append(line_ext)
                                f['legends2'][idx][icol]['title'] = db['tex'][key]
                                f['legends2'][idx][icol]['label'].append('')

            if meta.get('marklog2'):
                ax.axvline(x=np.log(2), color='black', alpha=0.85)
                trans = transforms.blended_transform_factory(ax.transData, ax.transAxes)
                ax.text(np.log(2), 0.90, '$\ln 2$', fontsize='small', color='black', alpha=0.75, ha='right',
                        va='center', rotation='vertical', transform=trans)
            if meta.get('marklog3'):
                ax.axvline(x=np.log(3), color='black', alpha=0.85)
                trans = transforms.blended_transform_factory(ax.transData, ax.transAxes)
                ax.text(np.log(3), 0.90, '$\ln 3$', fontsize='small', color='black', alpha=0.75, ha='right',
                        va='center', rotation='vertical', transform=trans)

            # ax.axvline(x=np.log(3), color='darkseagreen')
            # ax.text(np.log(3), 0.8, '$\ln 3$', fontsize='small', color='darkseagreen', ha='right', va='center', rotation='vertical',
            #         transform=trans)


            if not idx in f['axes_used']:
                f['axes_used'].append(idx)

            if axtitle := get_default(meta, 'axtitle'):
                if dbidx == 0:
                    if isinstance(axtitle, bool):
                        axtitle = get_title(dbval, subspec)
                    ax.set_title(axtitle,horizontalalignment='center', x=0.5,fontstretch="ultra-condensed")

        if figspec_title := get_figspec_title(meta, dbval, figspec):
            f['fig'].suptitle(figspec_title)

        # prettify_plot4(fmeta=f, lgnd_meta=axes_legends)
        suffix = ''
        suffix = suffix + '_normpage' if 'normpage' in meta and meta['normpage'] else suffix
        if filename := meta.get('filename'):
            f['filename'] = f"{meta['plotdir']}/{filename}"
        else:
            f['filename'] = "{}/{}_fig({})_sub({}){}".format(meta['plotdir'], meta['plotprefix'],
                                                       get_specvals(db, figspec, figvals),
                                                       get_specvals(db, subspec), suffix)
        f['figcount'] += 1
    return figs


def plot_divg_v2_fig3_sub3_line1(db, meta, figspec, subspec, linspec, algo_filter=None, state_filter=None,
                                 point_filter=None, figs=None, palette_name=None):
    if len(figspec) + len(subspec) + len(linspec) != 7:
        raise AssertionError("Must add to 7 elems: \n figspec {}\n subspec {}\n linespec {}")
    if 'mplstyle' in meta:
        plt.style.use(meta['mplstyle'])
    if 'mplstyle' in meta and 'slack' in meta['mplstyle']:
        # palette_name = "Spectral"
        if not palette_name:
            palette_name = "colorblind"
        # path_effects = [pe.SimpleLineShadow(offset=(0.5, -0.5), alpha=0.3), pe.Normal()]
        path_effects = None
    else:
        if not palette_name:
            palette_name = "colorblind"
        # path_effects = None
        path_effects = [pe.SimpleLineShadow(offset=(0.5, -0.5), alpha=0.3), pe.Normal()]

    prb_style = 'prb' in meta['mplstyle'] if 'mplstyle' in meta else False

    # legend_col_keys = list(itertools.chain(l1, [col for col in meta['legendcols'] if 'legendcols' in meta]))
    legend_col_keys = []
    if legendcols := meta.get('legendcols'):
        for col in legendcols:
            if not col in [l.split(':')[0] for l in subspec]:
                legend_col_keys.append(col)

    figprod = list(product(*get_keys(db, figspec)))  # All combinations of figspecs (names of parameters that iterate figures)
    subprod = list(product(*get_keys(db, subspec)))  # All combinations of subspecs (names of parameters that iterate subplots)
    linprod = list(product(*get_keys(db, linspec)))  # All combinations of linspecs (names of parameters that iterate lines)
    dirprod = list(product(db['keys']['algo'], db['keys']['state'], db['keys']['crono']))
    numfigs = len(figprod)
    numsubs = len(subprod)
    if figs is None:
        figs = [get_fig_meta(numsubs, meta=meta) for _ in range(numfigs)]

    for figkeys, f in zip(figprod, figs):
        print('- plotting figkeys: {}'.format(figkeys))
        dbval = None
        for idx, (subkeys, ax) in enumerate(zip(subprod, f['ax'])):
            popt = None
            pcov = None
            dbval = None
            print('-- plotting subkeys: {}'.format(subkeys))
            for dirkeys in dirprod:
                palette, lstyles = get_colored_lstyles(db, linspec, palette_name)
                for linkeys, color, lstyle in zip(linprod, palette, lstyles):
                    findlist = list(figkeys) + list(subkeys) + list(dirkeys) + list(linkeys) + [meta['groupname']]
                    datanode = [value['node']['data'] for key, value in db['dsets'].items() if
                                all(k in key for k in findlist)]
                    if len(datanode) != 1:
                        print("ERROR: found", len(datanode), "datanodes: ", datanode, " | findlist: ", findlist)
                        continue
                        # raise LookupError("Found incorrect number of datanodes")
                    print('-- plotting linkeys: {}'.format(linkeys))
                    datanode = datanode[0]
                    mmntnode = datanode.parent['measurements']
                    dbval = db['dsets'][datanode.name]
                    ndata = datanode['avg']['num'][()]
                    if np.min(ndata) < 10:
                        print("ndata < 10: in {}\n\t continuing".format(datanode))
                        continue
                    t = mmntnode['avg']['physical_time'][()]
                    y = datanode[meta['dsetname']][()]
                    n = datanode['avg']['num'][()][0]
                    idx_sat = find_saturation_idx3(t, dbval)
                    if meta.get('normpage'):
                        y /= page_entropy(dbval['L'])
                    if normalize := meta.get('normalize'):
                        y /= normalize
                    # if t[-1] == t[idx_sat]:
                    #     print('Time window is not long enough: saturation index = {} / {}', idx_sat, len(t))
                    #     continue
                    legendrow = get_legend_row(db=db, datanode=datanode, legend_col_keys=legend_col_keys)
                    print('linkeys:', linkeys)
                    if 'number' in meta['dsetname'] and linkeys == linprod[-1] and plot_divg_fig_sub_line.prb is None:
                        # Plot Luitz data
                        with h5py.File('external/raw_EE_NE_CE_distributions_random_XXX_chain.h5', 'r') as h5ext:
                            lenval = "{}".format(get_vals(db=dbval, keyfmt='L'))
                            lenkey = f"L{lenval}"
                            hist = h5ext[f'{lenkey}/W4.0']['hist[NE1][100]'][()]
                            edges = h5ext[f'{lenkey}/W4.0']['binedges[NE1][100]'][()]
                            bincentres = [(edges[j] + edges[j + 1]) / 2. for j in range(len(edges) - 1)]
                            line_ext, = ax.step(x=bincentres, y=hist, where='mid', label=None, color='gray', alpha=1.0,
                                                zorder=0)
                            for icol, (col, key) in enumerate(zip(legendrow, legend_col_keys)):
                                key, fmt = key.split(':') if ':' in key else [key, '']
                                f['legends'][idx][icol]['handle'].append(line_ext)
                                f['legends'][idx][icol]['title'] = db['tex'][key]
                                f['legends'][idx][icol]['label'].append(
                                    '\makebox[3ex][l]{PRB:102.100202 ' + 'L={},W=4.0'.format(lenval) + '}')
                                # if key in linkeys[0]:
                                #     f['legends'][idx][icol]['label'].append('\makebox[3ex][l]{PRB:102.100202 ' + 'L={},W=4.0'.format(lenval) + '}')
                                # else:
                                #     f['legends'][idx][icol]['label'].append('')

                            hist = h5ext[f'{lenkey}/W6.0']['hist[NE1][100]'][()]
                            edges = h5ext[f'{lenkey}/W6.0']['binedges[NE1][100]'][()]
                            bincentres = [(edges[j] + edges[j + 1]) / 2. for j in range(len(edges) - 1)]
                            line_ext, = ax.step(x=bincentres, y=hist, where='mid', label=None, color='gray', alpha=0.70,
                                                zorder=0)
                            for icol, (col, key) in enumerate(zip(legendrow, legend_col_keys)):
                                key, fmt = key.split(':') if ':' in key else [key, '']
                                f['legends'][idx][icol]['handle'].append(line_ext)
                                f['legends'][idx][icol]['title'] = db['tex'][key]
                                f['legends'][idx][icol]['label'].append(
                                    '\makebox[3ex][l]{PRB:102.100202 ' + 'L={},W=6.0'.format(lenval) + '}')
                                # if key in linkeys[0]:
                                #     f['legends'][idx][icol]['label'].append('\makebox[3ex][l]{PRB:102.100202 ' + 'L={},W=6.0'.format(lenval) + '}')
                                # else:
                                #     f['legends'][idx][icol]['label'].append('')

                            hist = h5ext[f'{lenkey}/W8.0']['hist[NE1][100]'][()]
                            edges = h5ext[f'{lenkey}/W8.0']['binedges[NE1][100]'][()]
                            bincentres = [(edges[j] + edges[j + 1]) / 2. for j in range(len(edges) - 1)]
                            line_ext, = ax.step(x=bincentres, y=hist, where='mid', label=None, color='gray', alpha=0.40,
                                                zorder=0)
                            for icol, (col, key) in enumerate(zip(legendrow, legend_col_keys)):
                                key, fmt = key.split(':') if ':' in key else [key, '']
                                f['legends'][idx][icol]['handle'].append(line_ext)
                                f['legends'][idx][icol]['title'] = db['tex'][key]
                                f['legends'][idx][icol]['label'].append(
                                    '\makebox[3ex][l]{PRB:102.100202 ' + 'L={},W=8.0'.format(lenval) + '}')

                                # if key in linkeys[0]:
                                # f['legends'][idx][icol]['label'].append('\makebox[3ex][l]{PRB:102.100202 L={},W=8.0}'.format(lenval))
                                # else:
                                #     f['legends'][idx][icol]['label'].append('')

                            plot_divg_fig_sub_line.prb = True

                    # Calculate the infinite time average (1/T) integral_0^T y(t) dt in the saturated interval
                    ytavg = np.mean(y[idx_sat:, :], axis=0)
                    hist, edges = np.histogram(ytavg, bins=meta['bins'], density=True)
                    bincentres = [(edges[j] + edges[j + 1]) / 2. for j in range(len(edges) - 1)]
                    line, = ax.step(x=bincentres, y=hist, where='mid', label=None, linewidth=1.25,
                                    color=color, path_effects=path_effects)
                    # line = ax.scatter(x=bincentres, y=hist, label=None,
                    #                    color=color, path_effects=path_effects)
                    ax.axvline(x=np.log(2), color='black', alpha=0.75)
                    trans = transforms.blended_transform_factory(ax.transData, ax.transAxes)
                    ax.text(np.log(2), 0.35, '$\ln 2$', fontsize='small', color='black', alpha=0.75, ha='right',
                            va='center', rotation='vertical',
                            transform=trans)
                    # ax.axvline(x=np.log(3), color='darkseagreen')
                    # ax.text(np.log(3), 0.8, '$\ln 3$', fontsize='small', color='darkseagreen', ha='right', va='center', rotation='vertical',
                    #         transform=trans)
                    for icol, (col, key) in enumerate(zip(legendrow, legend_col_keys)):
                        key, fmt = key.split(':') if ':' in key else [key, '']
                        f['legends'][idx][icol]['handle'].append(line)
                        f['legends'][idx][icol]['title'] = db['tex'][key]
                        f['legends'][idx][icol]['label'].append(col)

                    if not idx in f['axes_used']:
                        f['axes_used'].append(idx)

            if dbval:
                ax.set_title(get_title(dbval, subspec),
                             horizontalalignment='left', x=0.05,
                             fontstretch="ultra-condensed",
                             # bbox=dict(boxstyle='square,pad=0.15', facecolor='white', alpha=0.6)
                             )

        if not prb_style and dbval:
            f['fig'].suptitle('{} distribution\n{}'.format(meta['titlename'], get_title(dbval, figspec)))

        # prettify_plot4(fmeta=f, lgnd_meta=axes_legends)
        suffix = ''
        suffix = suffix + '_normpage' if 'normpage' in meta and meta['normpage'] else suffix
        suffix = suffix + '_loglog' if meta.get('timeselection') == 'lnlnt' else suffix
        f['filename'] = "{}/{}_divg_fig({})_sub({}){}".format(meta['plotdir'], meta['plotprefix'],
                                                              '-'.join(map(str, figkeys)),
                                                              '-'.join(map(str, get_keys(db, subspec))),
                                                              suffix)

    return figs


def plot_divg_fig_sub_line(db, meta, figspec, subspec, linspec, algo_filter=None, state_filter=None, point_filter=None,
                           figs=None, palette_name=None,dbidx=0,dbnum=1):
    if db['version'] == 2:
        specs = figspec + subspec + linspec
        nonv3spec = lambda x: x not in ['cstd', 'tstd', 'tgw8', 'cgw8']
        specs = list(filter(nonv3spec, specs))
        fig3 = specs[:3]
        sub3 = specs[3:6]
        lin1 = specs[6:7]
        return plot_divg_v2_fig3_sub3_line1(db=db, meta=meta, figspec=fig3, subspec=sub3, linspec=lin1,
                                            algo_filter=algo_filter, state_filter=state_filter,
                                            point_filter=point_filter, figs=figs, palette_name=palette_name)
    elif db['version'] == 3:
        return plot_divg_v3_fig_sub_line(db=db, meta=meta, figspec=figspec, subspec=subspec, linspec=linspec,
                                         algo_filter=algo_filter, state_filter=state_filter, point_filter=point_filter,
                                         figs=figs, palette_name=palette_name,dbidx=dbidx,dbnum=dbnum)
    else:
        raise NotImplementedError('database version not implemented:' + db['version'])
