from itertools import product
from pathlib import Path
from matplotlib.ticker import FixedLocator,MaxNLocator
import matplotlib.patheffects as pe
import seaborn as sns
from scipy.optimize import curve_fit
from scipy.stats.mstats import gmean
import h5py
from .tools import *


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


def plot_v3_tavg_fig_sub_line(db, meta, figspec, subspec, linspec, xaxspec, algo_filter=None, state_filter=None, point_filter=None, figs=None, palette_name=None,dbidx=0,dbnum=1):
    path_effects = [pe.SimpleLineShadow(offset=(0.35, -0.35), alpha=0.4), pe.Normal()]
    if 'mplstyle' in meta:
        plt.style.use(meta['mplstyle'])
        if 'slack' in meta['mplstyle']:
            path_effects = [pe.SimpleLineShadow(offset=(0.5, -0.5), alpha=0.3), pe.Normal()]

    # legend_col_keys = list(itertools.chain(l1, [col for col in meta['legendcols'] if 'legendcols' in meta]))
    legend_col_keys = linspec.copy()
    if legendcols := meta['legendcols']:
        for col in legendcols:
            c = col.split(':')[0] # Put a ? in the fmt string to skip
            if not c in [l.split(':')[0] for l in figspec + subspec + linspec + xaxspec]:
                legend_col_keys.append(col)


    figprod = list(product(*get_vals(db=db, keyfmt=figspec, filter=meta.get('filter'))))  # All combinations of figspecs values
    subprod = list(product(*get_vals(db=db, keyfmt=subspec, filter=meta.get('filter'))))  # All combinations of subspecs values
    linprod = list(product(*get_vals(db=db, keyfmt=linspec, filter=meta.get('filter'))))  # All combinations of linspecs values
    xaxprod = list(product(*get_vals(db=db, keyfmt=xaxspec, filter=meta.get('filter'))))  # All combinations of linspecs values
    numfigs = len(figprod)
    numsubs = len(subprod)

    numsubs *= dbnum
    subprod *= dbnum

    if figs is None:
        figs = [get_fig_meta(numsubs, meta=meta) for _ in range(numfigs)]

    for figvals, f in zip(figprod, figs):
        logger.debug('- plotting figs: {}'.format(figvals))
        dbval = None
        luitz = []
        for idx, (subvals, ax) in enumerate(zip(subprod, f['ax'])):
            idx_palette = idx % len(palette_name)
            if dbidx != int(idx / dbnum):
                continue
            popt = None
            pcov = None
            dbval = None
            logger.debug('-- plotting subs: {}'.format(subvals))
            for lidx, linvals in enumerate(linprod):
                palette, lstyles = get_colored_lstyles(db, linspec, palette_name, filter=None, idx=idx_palette)
                color = palette[-1]
                lstyle = lstyles[-1]
                logger.debug('--- plotting lins: {}'.format(linvals))
                xvals, yvals, pvals,mvals,gvals, evals,qvals, nvals = [], [],[], [], [], [], [], []
                legendrow = None
                for xaxvals in xaxprod:
                    logger.debug('--- plotting xaxs: {}'.format(xaxvals))
                    datanodes = match_datanodes(db=db, meta=meta, specs=figspec + subspec + linspec + xaxspec,
                                                vals=figvals + subvals + linvals + xaxvals)
                    if len(datanodes) != 1:
                        logger.warning(f"match: \n"
                                       f"\tspec:{[figspec + subspec + linspec]}\n"
                                       f"\tvals:{[figvals + subvals + linvals]}")
                        logger.warning(f"found {len(datanodes)} datanodes: {datanodes=}")
                    for datanode in datanodes:
                        # print(f'found datanode: {datanode}')
                        dbval = db['dsets'][datanode.name]
                        # ydata = datanode['avg'][meta['colname']][()]
                        # edata = datanode['ste'][meta['colname']][()]  # Use standard deviation for curve_fit
                        # tdata = datanode['avg']['physical_time'][()].astype(float)
                        # ndata = datanode['avg']['num'][()]
                        # raise LookupError("Found incorrect number of datanodes")
                        mmntnode = datanode.parent['measurements']
                        # dbval = db['dsets'][datanode.name]
                        # Get this datapoint
                        t = mmntnode['avg']['physical_time'][()].astype(float)
                        y = datanode[meta['dsetname']][()]
                        n = datanode['avg']['num'][()][0]
                        t = get_timepoints(t,dbval)
                        idx_num, idx_ent = t.idx_num_saturated, t.idx_ent_saturated
                        idx_sat = idx_num if 'number' in meta['dsetname'] else idx_ent
                        if meta.get('normpage'):
                            y /= midchain_page_entropy(dbval['vals']['L'])
                        elif meta.get('normsize'):
                            y /= dbval['vals']['L']
                        elif normalize := meta.get('normalize'):
                            y /= normalize
                        if meta.get('plotpinfty'):
                            if 'pinfty_number_entropies' in datanode.parent.parent:
                                p = datanode.parent.parent['pinfty_number_entropies']['data'][()]
                            elif 'number_probabilities' in datanode.parent.parent:
                                print('Warning: Could not find pinfty_number entropies. Calculating from scratch!')

                                # Take the shannon-entropy of infinite-time averaged number probabilities
                                probs = datanode.parent.parent['number_probabilities'][()]
                                ptavg = np.mean(probs[:, idx_sat:, :], axis=1)
                                ptavg = np.ma.masked_invalid(np.ma.masked_equal(np.abs(ptavg), 0))
                                with np.errstate(divide='ignore'):
                                    p = -np.sum(ptavg * np.log(ptavg), axis=0)

                            else:
                                raise LookupError('Could not find data to generate pinfty entropies')

                        # Calculate the infinite time average (1/T) integral_0^T y(t) dt
                        ytavg = np.mean(y[idx_sat:, :], axis=0)
                        if 'mean' in meta['ystats']:
                            yvals.append(np.mean(ytavg))
                            if meta.get('plotpinfty'):
                                pvals.append(np.mean(p))
                        if 'median' in meta['ystats']:
                            mvals.append(np.median(ytavg))
                            if meta.get('plotpinfty'):
                                pvals.append(np.median(p))
                        if 'gmean' in meta['ystats']:
                            gvals.append(gmean(np.where(ytavg < 1e-10, 1e-10, ytavg)))
                            if meta.get('plotpinfty'):
                                pvals.append(gmean(np.where(p < 1e-10, 1e-10, p)))
                        if meta.get('plotpinfty'):
                            qvals.append(np.std(p) / np.sqrt(len(p)))
                        xvals.append(get_vals(dbval, xaxspec))
                        evals.append(np.std(ytavg) / np.sqrt(len(ytavg)))
                        nvals.append(n)
                yvals = np.ravel(yvals)
                pvals = np.ravel(pvals)
                qvals = np.ravel(qvals)
                mvals = np.ravel(mvals)
                gvals = np.ravel(gvals)
                xvals = np.ravel(xvals)
                evals = np.ravel(evals)
                nvals = np.ravel(nvals)
                axn = ax
                if meta.get('ytwinx'):
                    if lidx == 0:
                        axn.annotate('', xy=(0.35, 0.85), xytext=(0.50, 0.85), xycoords='axes fraction',
                                     arrowprops=dict(color=color, arrowstyle='->'))
                    else :
                        axn = ax.twinx()
                        axn.set_box_aspect(1)
                        # axn.set_visible(True)
                        # axn.axis('off')
                        axn.set_frame_on(False)
                        axn.yaxis.set_visible(True)
                        axn.annotate('', xy=(0.65, 0.85), xytext=(0.50,0.85), xycoords='axes fraction',
                                arrowprops=dict(color=color, arrowstyle='->'))

                        # l = ax.get_ylim()
                        # ln = yvals.min(), yvals.max()
                        # fun = lambda x: ln[0] + (x - l[0]) / (l[1] - l[0]) * (ln[1] - ln[0])
                        # ticks = fun(ax.get_yticks())
                        axn.yaxis.set_major_locator(MaxNLocator(nbins=3))

                line2 = None
                if 'median' in meta['ystats']:
                    line, = axn.plot(xvals, mvals, marker=None,linestyle=':', color=color, path_effects=path_effects)
                if 'gmean' in meta['ystats']:
                    # line, = ax.plot(xvals, gvals, marker=None, linestyle='--', color=color, path_effects=path_effects)
                    line = axn.errorbar(x=xvals, y=gvals, yerr=evals, linestyle='--', color=color, path_effects=path_effects)
                if 'mean' in meta['ystats']:
                    line = axn.errorbar(x=xvals, y=yvals, yerr=evals, color=color, path_effects=path_effects)
                if meta.get('plotpinfty'):
                    line2 = axn.errorbar(x=xvals, y=pvals, yerr=qvals, linestyle='--', color=color, path_effects=path_effects)

                    # ax.fill_between(x=xvals, y1=yvals - evals, y2=yvals + evals, alpha=0.10, color=color)
                    # line, = ax.plot(xvals, yvals, marker=None, color=color, path_effects=path_effects)
                if not datanode:
                    continue
                legendrow = get_legend_row(db=db, datanode=datanode, legend_col_keys=legend_col_keys)
                icol = 0
                for icol, (col, key) in enumerate(zip(legendrow, legend_col_keys)):
                    key, fmt = key.split(':') if ':' in key else [key, '']
                    if '?' in fmt:
                        continue
                    f['legends'][idx][icol]['handle'].append(line)
                    f['legends'][idx][icol]['title'] = db['tex'][key]
                    f['legends'][idx][icol]['label'].append(col)
                f['legends'][idx][icol]['handle'].append(line)
                f['legends'][idx][icol]['title'] = 'RPS' if dbidx==0 else 'Néel'
                f['legends'][idx][icol]['label'].append('$\overline S_\mathrm{N}^\infty$')

                f['legends'][idx][icol]['handle'].append(line2)
                f['legends'][idx][icol]['title'] = 'RPS' if dbidx==0 else 'Néel'
                f['legends'][idx][icol]['label'].append('$\overline S(p_\mathrm{N}^\infty)$')


                if 'number' in meta['dsetname'] and meta.get('plotluitz') and legendrow is not None:
                    with h5py.File('external/raw_EE_NE_CE_distributions_random_XXX_chain.h5', 'r') as h5ext:
                        # Plot Luitz's data
                        # Let's only do this for the largest system size (16)
                        # and only after the last plot, so that this legend entry is last
                        lenvals = [4, 6, 8, 10, 12, 14, 16]  # Picks out some system sizes
                        lenkeys = [f"L{lenval}" for lenval in lenvals]
                        Wlist = []
                        Wlist = [8.0] if idx == 0 else Wlist
                        Wlist = [6.0] if idx == 1 else Wlist
                        Wlist = [] if idx == 2 else Wlist
                        Wlist = [] if idx == 3 else Wlist
                        graypalette = sns.color_palette('Greys', n_colors=len(Wlist)+2)[1:-1]  # Skip edge colors (too bright and too dark)
                        for W, color in zip(Wlist, graypalette):
                            size_NE,mean_NE, errr_NE, mean_NL,errr_NL = [],[],[],[],[]
                            for lenval, lenkey in zip(lenvals, lenkeys):
                                lwpath = f'{lenkey}/W{W}'
                                if h5ext.get(lwpath) is None:
                                    continue
                                if (idx, lwpath) in luitz:  # Already plotted on this ax
                                    continue
                                luitz.append((idx, lwpath))
                                mean_NE.append(h5ext[lwpath]['mean[NE1]'][()])
                                errr_NE.append(h5ext[lwpath]['error[NE1]'][()])
                                mean_NL.append(h5ext[lwpath]['mean[NE1_limit]'][()])
                                errr_NL.append(h5ext[lwpath]['error[NE1_limit]'][()])
                                size_NE.append(lenval)
                            line_ext = ax.errorbar(x=size_NE, y=mean_NE, yerr=errr_NE, label=None, color=color, alpha=1.0,
                                                path_effects=path_effects,
                                                linewidth=0.85,
                                                zorder=0)
                            ax.errorbar(x=size_NE, y=mean_NL, yerr=errr_NL, label=None, color=color, alpha=1.0,
                                                path_effects=path_effects,
                                                linewidth=0.85,
                                                linestyle='--',
                                                zorder=0)
                            if idx == 0 or idx == 1:  # Only on the leftmost subplot
                                # key, fmt = key.split(':') if ':' in key else [key, '']
                                f['legends2'][idx][0]['handle'].append(line_ext)
                                f['legends2'][idx][0]['label'].append(f'{W}')
                                if f['legends2'][idx][0]['title'] is None:
                                    # f['legends'][idx][0]['header'] = 'PRB\n102.100202'
                                    # f['legends'][idx][0]['header'] = 'Luitz et al.'
                                    f['legends2'][idx][0]['title'] = '$W$'

                                # f['legends2'][idx][1]['handle'].append(line_ext)
                                # f['legends2'][idx][1]['label'].append(f'{W}')
                                # if f['legends2'][idx][1]['title'] is None:
                                #     f['legends2'][idx][1]['title'] = '$W$'

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
                    lenval = 14  # Picks out the largest system size in Luitz's data
                    lenkey = f"L{lenval}"
                    with h5py.File('external/raw_EE_NE_CE_distributions_random_XXX_chain.h5', 'r') as h5ext:
                        graypalette = sns.color_palette('Greys', n_colors=5)[
                                      1:-1]  # Skip edge colors (too bright and too dark)
                        for W, color in zip([3.0, 6.0, 10.0], graypalette):
                            # for W,color in zip([10.0, 6.0, 3.0], graypalette):
                            lwpath = f'{lenkey}/W{W}'
                            if h5ext.get(lwpath) is None:
                                continue
                            if (idx, lwpath) in luitz:  # Already plotted on this ax
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
                            if idx == 0:  # Only on the leftmost subplot
                                # key, fmt = key.split(':') if ':' in key else [key, '']
                                f['legends2'][idx][0]['handle'].append(line_ext)
                                f['legends2'][idx][0]['label'].append(f'{lenval}')
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

                if not idx in f['axes_used']:
                    f['axes_used'].append(idx)

            if axtitle := get_default(meta, 'axtitle'):
                if dbidx == 0:
                    if isinstance(axtitle, bool):
                        axtitle = get_title(dbval, subspec)
                    ax.set_title(axtitle,horizontalalignment='center', x=0.5,fontstretch="ultra-condensed")
            # if dbval:
            #     ax.set_title(get_title(dbval, subspec),
            #                  horizontalalignment='left', x=0.05,
            #                  fontstretch="ultra-condensed",
            #                  # bbox=dict(boxstyle='square,pad=0.15', facecolor='white', alpha=0.6)
            #                  )


        if figspec_title := get_figspec_title(meta, dbval, figspec):
            f['fig'].suptitle(figspec_title)
        # prettify_plot4(fmeta=f, lgnd_meta=axes_legends)
        if not f['filename']:
            suffix = ''
            suffix = suffix + '_normpage' if 'normpage' in meta and meta['normpage'] else suffix
            f['filename'] = "{}/{}_tavg_fig({})_sub({}){}".format(meta['plotdir'], meta['plotprefix'],
                                                                  '-'.join(map(str, figvals)),
                                                                  '-'.join(map(str, get_keys(db, subspec))),
                                                                  suffix)
        suffix = ''
        suffix = suffix + '_normpage' if 'normpage' in meta and meta['normpage'] else suffix
        f['filename'] = "{}/{}_fig({})_sub({}){}".format(meta['plotdir'], meta['plotprefix'],
                                                       get_specvals(db, figspec, figvals),
                                                       get_specvals(db, subspec), suffix)
    return figs



def plot_tavg_fig_sub_line(db, meta, figspec, subspec, linspec, xaxspec, algo_filter=None, state_filter=None,
                           point_filter=None, figs=None, palette_name=None,dbidx=0,dbnum=1):
    if db['version'] == 2:
        raise ValueError("version 2 not supported")
    elif db['version'] == 3:
        return plot_v3_tavg_fig_sub_line(db=db, meta=meta, figspec=figspec, subspec=subspec, linspec=linspec,
                                         xaxspec=xaxspec,
                                         algo_filter=algo_filter, state_filter=state_filter, point_filter=point_filter,
                                         figs=figs, palette_name=palette_name,dbidx=dbidx,dbnum=dbnum)
    else:
        raise NotImplementedError('database version not implemented:' + db['version'])
