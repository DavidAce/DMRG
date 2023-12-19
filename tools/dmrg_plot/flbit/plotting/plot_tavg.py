from itertools import product
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator,MaxNLocator
import matplotlib.patheffects as pe
import seaborn as sns
from scipy.optimize import curve_fit
from scipy.stats.mstats import gmean
import h5py
from .tools import *
from itertools import islice
from scipy.signal import savgol_filter,find_peaks, find_peaks_cwt
from scipy import interpolate
from dataclasses import asdict

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

def batched(iterable, n):
    # batched('ABCDEFG', 3) --> ABC DEF G
    if n < 1:
        raise ValueError('n must be at least one')
    it = iter(iterable)
    while batch := list(islice(it, n)):
        yield batch

@dataclass
class Plots:
    x : np.ndarray       = None
    y : np.ndarray       = None
    e : np.ndarray       = None
    y_tavgs = []
    e_tavgs = []
    x_peaks : np.ndarray  = None
    y_peaks : np.ndarray  = None
    y_pinfty : np.ndarray  = None
    e_pinfty : np.ndarray  = None
    y_varx : np.ndarray    = None
    e_varx : np.ndarray    = None
    y_hartley : np.ndarray = None
    e_hartley : np.ndarray = None


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
                logger.debug('--- plotting lins: {}'.format(linvals))
                plots = Plots()
                xvals, yvals, pvals, mvals, gvals, evals, qvals, nvals, vvals, vevals, speak, ybavg, ebavg, ybmed = [], [], [], [], [], [], [], [], [], [], [], [], [], []
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
                        L = dbval['vals']['L']
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
                        tdata = t
                        t = get_timepoints(t,dbval)
                        idx_num, idx_ent = t.idx_num_saturated, t.idx_ent_saturated
                        idx_sat = idx_num if 'number' in meta['dsetname'] else idx_ent
                        # idx_sat = idx_ent
                        idx_len = len(t.time)
                        print(f'Saturation for {L=}: at {idx_num=} {idx_ent=}  ({idx_len}) total): time averaging {idx_len-idx_sat} points')
                        if meta.get('normpage'):
                            y /= midchain_page_entropy(dbval['vals']['L'])
                        elif meta.get('normsize'):
                            y /= dbval['vals']['L']
                        elif normalize := meta.get('normalize'):
                            y /= normalize

                        if meta.get('plotpeakln2'):
                            ytavg = np.nanmean(y[idx_sat:, :], axis=0)  # y dims: time and disorder realizations
                            hist, edges = np.histogram(ytavg, bins=meta['plotpeakln2_bins'], range=meta['plotpeakln2_range'], density=True)
                            bincentres = np.array([(edges[j] + edges[j + 1]) / 2. for j in range(len(edges) - 1)])
                            hist_smooth = savgol_filter(hist, 20, 4)  # window size 51, polynomial order 3
                            peak_idx,_ = find_peaks(hist_smooth, width=(10,None), prominence=(1.0,None))
                            print(peak_idx)
                            print(speak[-1])

                            if len(peak_idx)!=1:
                                raise AssertionError(f"Could not find exactly 1 peak:\n {peak_idx}")

                            fig, axs = plt.subplots()
                            axs.step(x=bincentres, y=hist, where='mid', label=None,color='black', path_effects=path_effects)
                            axs.step(x=bincentres, y=hist_smooth, where='mid', label=None,color='red', path_effects=path_effects)
                            axs.scatter(bincentres[peak_idx],hist_smooth[peak_idx], s=40,color='green', marker='v',zorder=50)
                            speak.append(bincentres[peak_idx])
                            # plt.show()

                        if meta.get('plotpinfty'):
                            if 'pinfty_number_entropies' in datanode.parent.parent:
                                # p = datanode.parent.parent['pinfty_number_entropies']['avg'][()]
                                # rhalf = int(n/2) if L==16 else int(n)
                                p = np.nanmean(datanode.parent.parent['pinfty_number_entropies']['data'][idx_sat:,:], axis=1)
                                print(f'shape p: {np.shape(p)}')
                                print(f'{p=}')
                            elif 'number_probabilities' in datanode.parent.parent:
                                print('Warning: Could not find pinfty_number entropies. Calculating from scratch!')
                                Lhalf = int(L/2)
                                rhalf = int(n/2) if L==16 else int(n)
                                realz = range(0, rhalf)
                                # realz = range(rhalf, n)
                                print('Reading shape')
                                d0,d1,d2,d3 = np.shape(datanode.parent.parent['number_probabilities'])
                                print('Allocating probs')
                                probs = np.empty(shape=(d0, d2, rhalf))
                                for r in batched(realz, 10000):
                                    print(f'reading number_probabilities {r[0]}/{rhalf}')
                                    probs[:, :, r] = datanode.parent.parent['number_probabilities'][:, Lhalf, :, r]

                                # Take the shannon-entropy of infinite-time averaged number probabilities
                                # probs = datanode.parent.parent['number_probabilities'][:, Lhalf, idx_sat:,:]
                                ptavg = np.nanmean(probs[:, idx_sat:, :], axis=1)
                                ptavg = np.ma.masked_invalid(np.ma.masked_equal(np.abs(ptavg), 0))
                                with np.errstate(divide='ignore'):
                                    p = -np.sum(ptavg * np.log(ptavg), axis=0)
                            elif 'number_probabilities2' in datanode.parent.parent:
                                print('Warning: Could not find pinfty_number entropies. Calculating from scratch!')
                                Lhalf = int(L/2)
                                rhalf = int(n/2) if L==16 else int(n)
                                realz = range(0, rhalf)
                                # realz = range(rhalf, n)
                                print('Reading shape')
                                d0,d1,d2,d3 = np.shape(datanode.parent.parent['number_probabilities'])
                                print('Allocating probs')
                                probs = np.empty(shape=(d0, d2, rhalf))
                                for r in batched(realz, 10000):
                                    print(f'reading number_probabilities {r[0]}/{rhalf}')
                                    probs[:, :, r] = datanode.parent.parent['number_probabilities'][:, Lhalf, :, r]

                                # Take the shannon-entropy of infinite-time averaged number probabilities
                                # probs = datanode.parent.parent['number_probabilities'][:, Lhalf, idx_sat:,:]
                                ptavg = np.nanmean(probs[:, idx_sat:, :], axis=(1,2))
                                ptavg = np.ma.masked_invalid(np.ma.masked_equal(np.abs(ptavg), 0))
                                with np.errstate(divide='ignore'):
                                    p = -np.sum(ptavg * np.log(ptavg), axis=0,keepdims=True)
                            else:
                                raise LookupError('Could not find data to generate pinfty entropies')
                        if meta.get('plotVarX'):
                            L = dbval['vals']['L']
                            posnode = datanode.parent.parent['nth_particle_position']
                            r0 = range(int(L // 4) - 1, int(L // 4) + 0)  # One particle to the left o the half-chain
                            r1 = range(int(L // 4) - 0, int(L // 4) + 1)  # One particle o the right of the half-chain
                            r2 = range(int(L // 4) - 1, int(L // 4) + 1)  # Two central particles on either side of the half-chain
                            y0 = np.atleast_3d(posnode['pos_variance_neel0'][r2, :, :])   # 1010|1010 index L//4 = 2 is nearest the midchain boundary
                            y1 = np.atleast_3d(posnode['pos_variance_neel1'][r2, :, :])   # 1010|1010 index L//4 = 2 is nearest the midchain boundary
                            # y = np.atleast_2d(np.concatenate([y0, y1], axis=2))
                            # Now interpret each particle is just more realizations at the midchain
                            y = np.concatenate([y0[0, :, :], y0[1, :, :], y1[0, :, :], y1[1, :, :]], axis=1)
                            # y = np.concatenate([y0[0, :, :], y1[1, :, :]], axis=1)
                            # y = np.concatenate([y0[0, :, :], y1[0, :, :]], axis=1)

                        if meta.get('plothartley'):
                            y = np.atleast_2d(datanode.parent.parent['hartley_number_entropies/data'][()]).T

                        # Calculate the infinite time average (1/T) integral_0^T y(t) dt
                        # reals = np.shape(y)[1]
                        # off = int(reals / 2)
                        # ext = int(reals / 5)
                        # esb = find_entropy_inftime_saturation_value_from_bootstrap(sdata=y, tdata=tdata, nbs=100, dsetname=meta['dsetname'])
                        # esb = find_entropy_inftime_saturation_value_from_bootstrap(sdata=y, tdata=tdata, nbs=meta.get('num-bootstraps', 100), dsetname=meta['dsetname'])
                        dsetname = 'varx' if meta.get('plotVarX') else meta['dsetname']
                        filenamejson = "{}/tsat_{}_L[{}]_x[{}]_w[{}]_f[{}].json".format(dbval['vals']['cachedir'],
                                                                                  dsetname,
                                                                                  dbval['vals']['L'],
                                                                                  dbval['vals']['x'],
                                                                                  dbval['vals']['w'],
                                                                                  dbval['vals']['f']
                                                                                  )
                        if meta.get('loadjson') and os.path.isfile(filenamejson):
                            with open(filenamejson, 'r') as fp:
                                tsb_json = json.load(fp)
                                tsb = entropy_saturation_bootstrap(**tsb_json)
                        else:
                            # tboot_idx_avg, tboot_idx_err, sboot_avg, sboot_err = get_entropy_saturation_from_bootstrap(ydata=s, nbs=100)
                            tsb = find_entropy_inftime_saturation_value_from_bootstrap(sdata=y, tdata=tdata, nbs=meta.get('num-bootstraps', 100), dsetname=dsetname)

                            if meta.get('savejson'):
                                with open(filenamejson, 'w') as fp:
                                    json.dump(asdict(tsb), fp, indent=4)



                        if 'mean' in meta['ystats']:
                            # if meta.get('plotVarX'):
                            #     yvals.append(np.mean(y))
                            # else:
                            yvals.append(tsb.sinf_full_avg)
                            if meta.get('plotpinfty'):
                                pvals.append(np.mean(p))
                        if 'median' in meta['ystats']:
                            mvals.append(tsb.sinf_full_med)
                            if meta.get('plotpinfty'):
                                pvals.append(np.median(p))
                        if 'gmean' in meta['ystats']:
                            gvals.append(gmean(np.where(ytavg < 1e-10, 1e-10, ytavg)))
                            if meta.get('plotpinfty'):
                                pvals.append(gmean(np.where(p < 1e-10, 1e-10, p)))
                        if meta.get('plotpinfty'):
                            qvals.append(np.std(p) / np.sqrt(len(p)))
                        # if meta.get('plotVarX'):
                        #     vevals.append(np.std(v) / np.sqrt(len(v)))
                        xvals.append(get_vals(dbval, xaxspec))
                        evals.append(tsb.sinf_full_err)
                        nvals.append(n)
                        ybavg.append(tsb.sinf_boot_avg)
                        ybmed.append(tsb.sinf_boot_med)
                        ebavg.append(tsb.sinf_boot_err)

                yvals = np.ravel(yvals)
                pvals = np.ravel(pvals)
                qvals = np.ravel(qvals)
                mvals = np.ravel(mvals)
                gvals = np.ravel(gvals)
                xvals = np.ravel(xvals)
                evals = np.ravel(evals)
                nvals = np.ravel(nvals)
                speak = np.ravel(speak)
                ybavg = np.ravel(ybavg)
                ybmed = np.ravel(ybmed)
                ebavg = np.ravel(ebavg)
                if meta.get('ytwinx'):
                    if lidx == 0:
                        ax.annotate('', xy=(0.35, 0.85), xytext=(0.50, 0.85), xycoords='axes fraction',
                                     arrowprops=dict(color=color, arrowstyle='->'))
                    else :
                        ax = ax.twinx()
                        ax.set_box_aspect(1)
                        # ax.set_visible(True)
                        # ax.axis('off')
                        ax.set_frame_on(False)
                        ax.yaxis.set_visible(True)
                        ax.annotate('', xy=(0.65, 0.85), xytext=(0.50,0.85), xycoords='axes fraction',
                                arrowprops=dict(color=color, arrowstyle='->'))

                        # l = ax.get_ylim()
                        # ln = yvals.min(), yvals.max()
                        # fun = lambda x: ln[0] + (x - l[0]) / (l[1] - l[0]) * (ln[1] - ln[0])
                        # ticks = fun(ax.get_yticks())
                        ax.yaxis.set_major_locator(MaxNLocator(nbins=3))

                line2 = None
                palette, lstyles = get_colored_lstyles(db, xaxspec, palette_name, filter=None, idx=idx_palette)
                # color = palette[-1]
                # lstyle = lstyles[-1]
                if 'median' in meta['ystats']:
                    line, = ax.plot(xvals, mvals, marker=None,linestyle=':', color=color, path_effects=path_effects)
                if 'gmean' in meta['ystats']:
                    # line, = ax.plot(xvals, gvals, marker=None, linestyle='--', color=color, path_effects=path_effects)
                    line = ax.errorbar(x=xvals, y=gvals, yerr=evals, linestyle='--', color=color, path_effects=path_effects)
                if 'mean' in meta['ystats']:
                    line, = ax.plot(xvals, yvals, marker=None, color='black', linewidth=1.0, path_effects=None)
                    # ax.plot(xvals, ybmed, color='red', linestyle=':', linewidth=1.0, path_effects=path_effects)
                    # ax.errorbar(x=xvals, y=ybavg, yerr=ebavg, color='blue', linestyle='-', capsize=1.0, path_effects=path_effects)
                    for xval,yval,eval,color in zip(xvals, yvals, evals,palette):
                        ax.errorbar(x=xval, y=yval, yerr=eval, color=color,linestyle='none',capsize=2.0, path_effects=path_effects,zorder=10)

                if meta.get('plotpeakln2'):
                    ax.plot(xvals, speak, marker=None, color='gray', path_effects=None)
                    for xval, yval, color in zip(xvals, speak, palette):
                        ax.scatter(x=xval, y=yval, color=color, marker='^', path_effects=path_effects, zorder=10)

                if meta.get('plotpinfty'):
                    line2, = ax.plot(xvals, pvals, marker=None, color='black', linestyle=':', path_effects=None)
                    for xval,pval,qval,color in zip(xvals, pvals, qvals,palette):
                        ax.errorbar(x=xval, y=pval, yerr=qval, color=color,linestyle='none',capsize=2.0, path_effects=path_effects,zorder=10)
                    # line2 = ax.errorbar(x=xvals, y=pvals, yerr=qvals, linestyle='--', color=color, path_effects=path_effects)

                # if meta.get('plotVarX'):
                #     line2 = ax.errorbar(x=xvals, y=vvals, yerr=vevals, linestyle='--', color=color,
                #                          path_effects=path_effects)

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
                # f['legends'][idx][icol]['handle'].append(line)
                # f['legends'][idx][icol]['title'] = 'RPS' if dbidx==0 else 'Néel'
                # f['legends'][idx][icol]['label'].append('$\overline S_\mathrm{N}^\infty$')
                #
                # f['legends'][idx][icol]['handle'].append(line2)
                # f['legends'][idx][icol]['title'] = 'RPS' if dbidx==0 else 'Néel'
                # f['legends'][idx][icol]['label'].append('$\overline S(p_\mathrm{N}^\infty)$')


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
                if 'entanglement' in meta['dsetname'] and meta.get('plotluitz') and legendrow is not None:
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
