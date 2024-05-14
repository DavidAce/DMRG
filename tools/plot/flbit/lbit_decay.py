import seaborn as sns
from matplotlib import transforms

from common.io.h5ops import *
from plotting.tools import *
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import matplotlib.patheffects as pe
from pathlib import Path
import logging
from numba import njit
from matplotlib.ticker import MaxNLocator
from matplotlib.ticker import LogLocator, \
    LogFormatter, LogFormatterExponent, LogFormatterSciNotation, \
    LogFormatterMathtext, NullFormatter, MultipleLocator, MaxNLocator, \
    FixedLocator, LinearLocator, ScalarFormatter

@njit(parallel=True, cache=True)
def fit_exp_stretched(x, C, xi, beta):
    return C * np.exp(-(x / xi) ** beta)


@njit(parallel=True, cache=True)
def fit_log_stretched(x, C, xi, beta):
    return np.log(C) - (x / xi) ** beta


@njit(parallel=True, cache=True)
def fit_exp(x, C, xi):
    return C * np.exp(-(x / xi))


@njit(parallel=True, cache=True)
def fit_log(x, C, xi):
    return np.log(C) - (x / xi)


datafiles = [
    # 'data/mbl_10003-uniform500.h5',
    # 'data/mbl_10003-normal500.h5',
    # 'data/mbl_10003-choked10.h5',
    # 'data/mbl_10003-choked100.h5',
    # 'data/mbl_10003-squared500.h5',
    # 'data/mbl_10003-constricted50.h5',
    # 'data/mbl_10003-unconstricted50.h5',
    # 'data/mbl_10003-constricted50-L20-u5-f0.45.h5',
    # 'data/mbl_10003-constricted50-L24-u5-f0.45.h5',
    # 'data/mbl_10003-constricted50-random-L12-u[2-5]-f0.45.h5',
    # 'data/mbl_10003-constricted50-random-L16-u[2-5]-f0.45.h5',
    # 'data/mbl_10003-constricted50-random-L20-u[2-5]-f0.45.h5',
    # 'data/mbl_10003-constricted50-random-L24-u[2-5]-f0.45.h5',
    # 'data/mbl_10003-constricted50-random-L24-u[2-5]-f[0.25-0.45].h5',
    # 'data/mbl_10003-constricted50-random-L32-u[2-8]-f[0.25-0.45].h5',
    # 'data/mbl_10003-constricted10-random-L32-u[8]-f[0.2-0.5].h5',
    # 'data/mbl_10003-constricted50-L24-u8-f[0.2-0.5].h5',
    # 'data/mbl_10003-constricted50-L24-u8-f0.5.h5', # <----- This one is in the draft
    # 'data/mbl_10004-constricted1000-L24-u8-f0.5.h5',
    'data/mbl_10004-constricted1000-random-L24-u8-f0.5.h5',
    # 'data/mbl_10004-uniform50-L24-u8-f0.5.h5',
    # 'data/mbl_10004-uniform50-L24-u4-f0.25.h5',
    # 'data/mbl_10004-uniform50-L24-u6-f0.3.h5', # <----- This one is in the draft
    'data/mbl_10004-uniform50-L24-u8-f0.5.h5',
]

lbitpath = 'fLBIT/model/lbits'
mplstyle = Path('common/stylesheets/prl.mplstyle', )
# mplstyle = '../common/stylesheets/slack.mplstyle',
plotdir = Path("plots/{}".format(mplstyle.stem))
plotdir.mkdir(parents=True, exist_ok=True)


def get_lbit_decay(x, y, fit, p0, bounds):
    ylog = np.log(np.abs(y))
    popt, pcov = curve_fit(fit, xdata=x, ydata=ylog, p0=p0, bounds=bounds)
    pstd = np.sqrt(np.diag(pcov))
    return popt, pstd


def plot_decay_fit(ax, xdata, ydata, p0, bounds, color, path_effects):
    popt, pstd = get_lbit_decay(x=xdata, y=ydata, fit=fit_log, p0=p0, bounds=bounds)
    ax.plot(xdata, fit_exp(xdata, *popt), label=None, color=color, path_effects=path_effects)
    print(popt, pstd)
    return np.append(popt, [np.nan]), np.append(pstd, [np.nan])


def plot_decay_fit_stretched(ax, xdata, ydata, p0, bounds, color, path_effects):
    popt, pstd = get_lbit_decay(x=xdata, y=ydata, fit=fit_log_stretched, p0=p0, bounds=bounds)
    ax.plot(xdata, fit_exp_stretched(xdata, *popt), label=None, color=color, path_effects=path_effects)
    print(popt, pstd)
    return popt, pstd


def plot_decay(meta, figs=None, palette=None):
    logging.basicConfig(level=logging.WARN)
    logger_lbit = logging.getLogger('lbit').setLevel(level=logging.WARN)
    logger_tools = logging.getLogger('tools').setLevel(level=logging.INFO)

    if palette is None:
        palette = itertools.cycle(sns.color_palette())
    #
    if 'slack' in meta['filename']:
        palette_name = "Spectral_r"
        pe_expfit = [pe.SimpleLineShadow(offset=(-0.3, -0.5), alpha=0.4), pe.Normal()]
    else:
        palette_name = "colorblind"
        pe_expfit = None

    h5data = h5open(meta['datafile'], 'r')
    frange = h5data["{}".format(lbitpath)].attrs["u_fmix"][()]
    urange = h5data["{}".format(lbitpath)].attrs["u_depth"][()]  # formerly u_depth
    data = h5data["{}/pata".format(lbitpath)]

    fmix_range = range(0, len(frange), 2)
    ucol_range = range(0, len(urange), 2)
    numplots = len(ucol_range)
    if figs is None:
        figs = [get_fig_meta(numplots, meta=meta)]
    f = figs[-1]
    figin = None
    fits = {'u': [], 'color': [], 'f': [], 'y0': [], 'xi': [], 'beta': []}
    for idx, (col, ax) in enumerate(zip(ucol_range, f['ax'])):
        current_palette = itertools.cycle(sns.color_palette())
        fit = {'u': urange[col], 'color': [], 'f': [], 'y0': [], 'xi': [], 'beta': []}
        for i, row in enumerate(fmix_range):
            color = next(current_palette)
            # decay  = np.mean(data[row,col,:, :, :], axis=(0,1))
            pupps = np.abs(data[row, col, :, :, :])  # Permuted supports
            decay = np.mean(pupps, axis=(0, 1))
            error = np.std(pupps, axis=(0, 1))
            e = error[np.abs(decay) > 1e-16]
            y = decay[np.abs(decay) > 1e-16]
            x = np.array(range(len(y)))
            print(e)

            # popt, pstd = plot_decay_fit_stretched(ax=ax, xdata=x, ydata=y,
            #                              p0=(1.0, 1.0, 1.0), bounds=(None, None),
            #                              color='red', path_effects=pe_expfit)
            # popt, pstd = plot_decay_fit(ax=ax, xdata=x, ydata=y,
            #                              p0=(1.0, 1.0), bounds=(None, None),
            #                              color='blue', path_effects=pe_expfit)
            popt, pstd = plot_decay_fit(ax=ax, xdata=x[2:], ydata=y[2:],
                                        p0=(1.0, 1.0), bounds=(None, None),
                                        color='red', path_effects=pe_expfit)
            # popt, pstd = plot_decay_fit(ax=ax, xdata=x[4:], ydata=y[4:],
            #                              p0=(1.0, 1.0), bounds=(None, None),
            #                              color='green', path_effects=pe_expfit)
            # popt, pstd = print_decay_fit(ax, fit=fit_exp_stretched, xdata=x[4:], ydata=y[4:], p0=(0.1, 0.5, 1.0),
            #                              bounds=([0.0, 0.0, 0.5], [1.0, np.inf, np.inf]), color=color, path_effects=pe_expfit)
            line, = ax.plot(x, y, marker='o', linestyle='None', color=color, path_effects=pe_expfit)
            # line = ax.errorbar(x=x, y=y,yerr=e, marker='o', linestyle='None', color=color, path_effects=pe_expfit)

            # popt,pstd = get_lbit_decay_stretched(y)
            fit['color'].append(color)
            fit['f'].append(frange[row])
            fit['y0'].append(popt[0])
            fit['xi'].append(popt[1])
            fit['beta'].append(popt[2])
            # cls_fit[row, col] = popt[1]
            # cls_std[row, col] = pstd[1]
            # cls_bet[row, col] = popt[2]
            # cls_1pc[row, col] = popt[1] * (np.log(popt[0] / (1.0 / 100.0)) ** (1.0 / popt[2]))
            # ax.set_title("$d_u={}$".format(urange[col]), bbox={'facecolor': 'white', 'pad': 2})
            # if insetpos := meta.get('inset').get('pos'):
            # if not axin:
            #     axin = ax.inset_axes([0.175,0.10,0.35,0.35])
            # axin.plot(popt[1], popt[2], marker='o',color=color)
            # axin.tick_params(labelsize=6, length=2, pad=1.0)
            # axin.set_yticks([1.0, 1.4])
            # axin.set_xticks([0.2, 1.0])
            # bbox = {'edgecolor': 'white', 'facecolor': 'white', 'pad': 2, 'alpha':0.0}
            # axin.set_ylabel('$\\beta$', labelpad=-10, fontsize=7, bbox=bbox)
            # axin.set_xlabel('$\\xi$', labelpad=-5, fontsize=7, bbox=bbox)

            # axin.set_xticklabels([],fontsize=5)
            # axin.set_yticklabels([],fontsize=5)

            icol = 0
            if 'f' in meta['legendcols']:
                f['legends'][idx][icol]['handle'].append(line)
                f['legends'][idx][icol]['label'].append('{:.2f}'.format(frange[row]))
                f['legends'][idx][icol]['title'] = '$f$'
                icol += 1
            if 'xi' in meta['legendcols']:
                f['legends'][idx][icol]['handle'].append(line)
                f['legends'][idx][icol]['label'].append('{:.2f}'.format(popt[1]))
                f['legends'][idx][icol]['title'] = '$\\xi$'
                icol += 1
            if 'beta' in meta['legendcols']:
                f['legends'][idx][icol]['handle'].append(line)
                f['legends'][idx][icol]['label'].append('{:.2f}'.format(popt[2]))
                f['legends'][idx][icol]['title'] = '$\\beta$'
                icol += 1
        fits['color'].append(fit['color'])
        fits['f'].append(fit['f'])
        fits['y0'].append(fit['y0'])
        fits['xi'].append(fit['xi'])
        fits['beta'].append(fit['beta'])

        f['axes_used'].append(idx)
    if meta['insetfit']:
        if not figin:
            figin = f['fig'].add_axes([0.50, 0.20, 0.15, 0.15])
        for uidx, (xi, beta, color) in enumerate(zip(fits['xi'], fits['beta'], fits['color'])):
            if uidx + 1 == len(fits['xi']):
                figin.scatter(xi, beta, color=color, marker='o')
        figin.tick_params(labelsize=7, length=2, pad=2.0)
        figin.set_yticks([1.0, 1.4])
        figin.set_xticks([0.1, 1.0])
        bbox = {'edgecolor': 'white', 'facecolor': 'white', 'pad': 2, 'alpha': 0.0}
        figin.set_ylabel('$\\beta$', labelpad=-10, fontsize=7, bbox=bbox)
        figin.set_xlabel('$\\xi$', labelpad=-5, fontsize=7, bbox=bbox)

    return figs


def plot_lbits(meta, figs=None, color=None):
    h5data = h5open(meta['datafile'], 'r')
    if color is None:
        color = next(itertools.cycle(sns.color_palette()))

    frange = h5data["{}".format(lbitpath)].attrs["f_mixer"][()]
    urange = h5data["{}".format(lbitpath)].attrs["u_layer"][()]
    lbits = h5data["{}".format(lbitpath)]["data"][()]
    decay = h5data["{}/decay".format(lbitpath)]  # formerly "curves"
    fields = h5data["fLBIT/model/hamiltonian"]['J1_rand']
    print('frange: {}'.format(frange))
    print('urange: {}'.format(urange))

    # lbits = lbits[-7::7, -1::1, :, :, :]
    # frange = frange[-7::7]
    # urange = urange[-1::1]
    # lbits = lbits[0::4, 0::1, :, :, :]
    # frange = frange[0::4]
    # urange = urange[::1]

    print('frange: {}'.format(frange))
    print('urange: {}'.format(urange))

    flen, ulen, reps, rows, cols = np.shape(lbits)
    numsubs = flen * ulen
    if figs is None:
        figs = [get_fig_meta(numsubs, meta=meta)]
    f = figs[-1]

    if 'slack' in meta['filename']:
        palette_name = "Spectral_r"
        pe_lbits = [pe.SimpleLineShadow(offset=(-0.3, -0.5), alpha=0.4), pe.Normal()]
    else:
        palette_name = "colorblind"
        pe_lbits = None

    pe_expfit = [pe.Stroke(linewidth=2, foreground='white'), pe.Normal()]
    pe_fields = [pe.Stroke(linewidth=2, foreground='white'), pe.Normal()]

    for fidx, fval in enumerate(frange):
        for uidx, uval in enumerate(urange):
            idx = np.ravel_multi_index([[fidx], [uidx]], dims=(flen, ulen), )[0]
            print('idx {} | fidx {} | uidx {}'.format(idx, fidx, uidx))
            ax = f['ax'][idx]
            x = np.array(range(cols))
            xmid = int(rows / 2)
            maxreps = np.min([reps, 50])
            # C, xi, beta, yfit, _, _ = get_lbit_fit(x, decay[fidx, uidx], beta=True, ymin=1e-10)
            # fit = get_lbit_fit(x, decay[fidx, uidx], beta=True, ymin=1e-10)
            fit = get_lbit_fit_data(x=x, y=np.atleast_2d(decay[fidx, uidx]).T, e=None,
                              beta=meta.get('fit-beta', True),
                              ymin=meta.get('fit-ymin', 1e-16),
                              )
            # print(np.shape(decay))
            # print(popt, pstd)
            #
            # ax.plot(x + xmid, fit_exp_stretched(x, popt[0], popt[1], popt[2]), label=None, color='green',
            #         path_effects=pe_expfit, zorder=60, )
            # ax.plot(x + xmid, exp_fit(x, popt[0], popt[1]), label=None, color=color, path_effects=pe_expfit)

            for site in [xmid]:
                l = lbits[fidx, uidx, range(maxreps), [site], :].T
                y = np.mean(l, axis=1)
                x = range(len(y))
                e = np.std(l, axis=1)

                line, = ax.plot(x, y, linewidth=1.0, color=color, alpha=0.0)
                # ax.xaxis.set_major_locator(MaxNLocator(integer=True))
                # f['legends'][idx][0]['handle'].append(line)
                # f['legends'][idx][0]['label'].append('{}'.format(cols))
                # f['legends'][idx][0]['title'] = '$L$'
                if 'f' in meta['legendcols']:
                    f['legends'][idx][1]['handle'].append(line)
                    f['legends'][idx][1]['label'].append('${:.2f}$'.format(fval))
                    f['legends'][idx][1]['title'] = '$f$'
                if 'u' in meta['legendcols']:
                    f['legends'][idx][2]['handle'].append(line)
                    f['legends'][idx][2]['label'].append('${}$'.format(uval))
                    f['legends'][idx][2]['title'] = '$d_u$'
                if 'xi' in meta['legendcols']:
                    f['legends'][idx][3]['handle'].append(line)
                    f['legends'][idx][3]['label'].append('$\quad{:.2f}$'.format(fit.xi))
                    f['legends'][idx][3]['title'] = '$\quad\\xi_\\tau$'
                if 'xi' in meta['legendcols']:
                    f['legends'][idx][4]['handle'].append(line)
                    f['legends'][idx][4]['label'].append('${:.2f}$'.format(fit.beta))
                    f['legends'][idx][4]['title'] = '$\\beta$'

                for ridx in range(maxreps):
                    y = lbits[fidx, uidx, ridx, [site], :].T
                    y[y < 1e-15] = np.nan
                    ax.plot(x, y, alpha=0.30, color=color)

                    if not idx in f['axes_used']:
                        f['axes_used'].append(idx)

            if meta['fieldinset']:
                trans = transforms.blended_transform_factory(ax.transData, ax.transAxes)
                fieldcolor = 'steelblue'
                for site in range(rows):
                    ax.arrow(site, 0.12, 0.0, fields[site] / 15.0, width=0.025, length_includes_head=False,
                             head_width=0.00, head_length=0.00,
                             path_effects=pe_fields, zorder=50,
                             transform=trans,
                             color=fieldcolor,
                             # color='salmon'
                             )
                ax.arrow(0, 0.12, rows, 0, width=0.001, length_includes_head=False,
                         head_width=0.00, head_length=0.00,
                         zorder=60, transform=trans,
                         color=fieldcolor
                         # color='salmon'
                         )

    return figs
def plot_lbits_new(meta, figs=None, color=None):
    h5data = h5open(meta['datafile'], 'r')
    if color is None:
        color = next(itertools.cycle(sns.color_palette()))
        print(color)
    u_fmix  = h5data[lbitpath].attrs["u_fmix"][0]
    u_depth = h5data[lbitpath].attrs["u_depth"][0]
    corrmat   = np.abs(h5data[lbitpath]["corrmat"][()])
    decay   = np.abs(h5data[f"{lbitpath}/corravg"][()])  # formerly "curves"
    fields  = h5data["fLBIT/model/hamiltonian"]['J1_rand'][()].astype(float)
    print(f'{u_fmix=}')
    print(f'{u_depth=}')

    reps, rows, cols = np.shape(corrmat)
    if figs is None:
        figs = [get_fig_meta(1, meta=meta)]
    f = figs[-1]

    pe_expfit = [pe.Stroke(linewidth=2, foreground='white'), pe.Normal()]
    pe_fields = [pe.Stroke(linewidth=2, foreground='white'), pe.Normal()]
    ax = f['ax'][0]
    # ax.yaxis.set_label_coords(-0.17, 0.4)

    x = np.arange(cols)
    maxreps = np.min([reps, 25])
    lbavg = get_lbit_avg(corrmat=corrmat, site=meta.get('lbit-site'),
                         mean=meta.get('lbit-mean'))
    fit = get_lbit_fit_data(x=x, y=np.atleast_2d(lbavg.full).T, e=np.atleast_2d(lbavg.stdv).T,
                      beta=meta.get('fit-beta'),
                      # beta=meta.get('fit-beta', True),
                      ymin=meta.get('fit-ymin', 1e-20),
                      )
    # print(np.shape(decay))
    # print(popt, pstd)
    #
    # ax.plot(x + xmid, fit_exp_stretched(x, popt[0], popt[1], popt[2]), label=None, color='green',
    #         path_effects=pe_expfit, zorder=60, )
    # ax.plot(x + xmid, exp_fit(x, popt[0], popt[1]), label=None, color=color, path_effects=pe_expfit)

    for site in [meta.get('lbit-site')]:
        l = corrmat[range(maxreps), [site], :].T
        y = np.mean(l, axis=1)
        x = np.arange(len(y))
        e = np.std(l, axis=1)
        print(fit.yfit)
        line, = ax.plot([], [], linewidth=1.5, color=color, alpha=1.0)
        # line, = ax.plot(x, fit.yfit, linewidth=1.5, color=color, alpha=1.0, path_effects=pe_expfit, zorder=4)
        # line, = ax.plot(x, y, linewidth=1.0, color=color, alpha=1.0)

        # ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        # f['legends'][idx][0]['handle'].append(line)
        # f['legends'][idx][0]['label'].append('{}'.format(cols))
        # f['legends'][idx][0]['title'] = '$L$'
        idx = 0
        if 'wi' in meta['legendcols']:
            f['legends'][idx][1]['handle'].append(line)
            f['legends'][idx][1]['label'].append(meta['wilabel'])
            f['legends'][idx][1]['title'] = '$w_i$'
        if 'f' in meta['legendcols']:
            f['legends'][idx][2]['handle'].append(line)
            f['legends'][idx][2]['label'].append('${:.1f}$'.format(u_fmix))
            f['legends'][idx][2]['title'] = '$f$'
        if 'u' in meta['legendcols']:
            f['legends'][idx][3]['handle'].append(line)
            f['legends'][idx][3]['label'].append('${}$'.format(u_depth))
            f['legends'][idx][3]['title'] = '$d_u$'
        if 'xi' in meta['legendcols']:
            f['legends'][idx][4]['handle'].append(line)
            f['legends'][idx][4]['label'].append('${:.1f}$'.format(fit.xi))
            f['legends'][idx][4]['title'] = '$\\xi_\\tau$'
        if 'beta' in meta['legendcols']:
            f['legends'][idx][5]['handle'].append(line)
            f['legends'][idx][5]['label'].append('${:.1f}$'.format(fit.beta))
            f['legends'][idx][5]['title'] = '$\\beta$'

        for ridx in range(maxreps):
            y = corrmat[ridx, [site], :].T
            ax.plot(x, y, alpha=0.5, linewidth=0.4, color=color)

            if not idx in f['axes_used']:
                f['axes_used'].append(idx)

    if meta['fieldinset']:
        trans = transforms.blended_transform_factory(ax.transData, ax.transAxes)
        fieldcolor = 'steelblue'
        for site in range(rows):
            ax.arrow(site, 0.12, 0.0, fields[site] / 10.0, width=0.025, length_includes_head=False,
                     head_width=0.00, head_length=0.00,
                     path_effects=pe_fields, zorder=50,
                     transform=trans,
                     color=fieldcolor,
                     # color='salmon'
                     )
        ax.arrow(0, 0.12, rows, 0, width=0.0005, length_includes_head=False,
                 head_width=0.00, head_length=0.00,
                 zorder=60, transform=trans,
                 color=fieldcolor
                 # color='salmon'
                 )

    return figs


if __name__ == '__main__':
    subplots1x1 = {
        'top': 1.0,
        'bottom': 0.28,
        'left': 0.20,
        'right': 1.0,
        'wspace': 0,
        'hspace': 0,
    }
    default = {
        'box_aspect': 0.46,
        'figsize': (3.404 * 0.5, 3.404*0.25),
        'subspec_title': False,
        'figspec_title': False,
        'legendoutside' : False,
        'legendcollect' : False,
        'subplots': subplots1x1,
    }
    subplots_halfheight = {
        'top': 1.0,
        'bottom': 0.15,
        'left': 0.17,
        'right': 1.0,
        'wspace': 0,
        'hspace': 0,
    }
    default_halfheight = {
        'box_aspect': 0.46,
        'figsize': (3.404 * 5.0/9, 3.404*0.25),
        'subspec_title': False,
        'figspec_title': False,
        'legendoutside' : False,
        'legendcollect' : False,
        'subplots': subplots_halfheight,
    }



    meta_lbit_log_fit_identity_bigc = {
        'default': default,
        'figsize': (1.702, 1.702),
        'filename': "{}/lbit-log-fit".format(plotdir),
        'datafile': 'data/mbl_10004-constricted50-random-bigC-L16-u8-f0.5.h5',
        'suptitle': None,  # 'l-bit fit $y=y_0 \exp[-(x/\\xi)^\\beta]$',
        'ylabel': '$\\bar O(i,j)$',
        'xlabel': '$|i-j|$',
        'yticks': [1e-0, 1e-4, 1e-8, 1e-12],
        'xticks': [0, 2, 4, 6, 8],
        'sharex': 'all',
        'sharey': 'all',
        'yscale': 'log',
        'xscale': 'linear',
        'insetfit': True,
        'legendcols': ['f'],
        'legendoutside': True,
        'legendcollect': True,
        'legendlocation': 'center',
        'owr_pad': 1.5,  # Make the left-most subplot wider by this factor to make space for ylabel and ticklabels
        # 'legendcols': ['f', 'xi', 'beta'],
        # 'legend_in_subfig': False
    }

    meta_lbit_log_fit_expdecay = {
        'default': default,
        'figsize': (1.702, 1.702),
        'filename': "{}/lbit-log-fit".format(plotdir),
        'datafile': 'data/mbl_10004-constricted50-random-L16-u8-f0.5.h5',
        'suptitle': None,  # 'l-bit fit $y=y_0 \exp[-(x/\\xi)^\\beta]$',
        'ylabel': '$\\bar O(i,j)$',
        'xlabel': '$|i-j|$',
        'yticks': [1e-0, 1e-4, 1e-8, 1e-12],
        'xticks': [0, 2, 4, 6, 8],
        'sharex': 'all',
        'sharey': 'all',
        'yscale': 'log',
        'xscale': 'linear',
        'insetfit': True,
        'legendcols': ['f'],
        'legendoutside': True,
        'legendcollect': True,
        'legendlocation': 'center',
        'owr_pad': 1.5,  # Make the left-most subplot wider by this factor to make space for ylabel and ticklabels
        # 'legendcols': ['f', 'xi', 'beta'],
        # 'legend_in_subfig': False
    }

    meta_lbit_exp_fit_blocked = {
        'default': default,
        'figsize': (1.702, 1.702),
        'filename': "{}/lbit-exp-fit".format(plotdir),
        'datafile': 'data/mbl_10004-constricted500-random-L32-u8-f0.5.h5',
        'suptitle': None,  # 'l-bit fit $y=y_0 \exp[-(x/\\xi)^\\beta]$',
        'ylabel': '$\\bar O(i,j)$',
        'xlabel': '$|i-j|$',
        'yticks': [1e-0, 1e-4, 1e-8, 1e-12],
        'xticks': [0, 2, 4, 6, 8],
        'sharex': 'all',
        'sharey': 'all',
        'yscale': 'log',
        'xscale': 'linear',
        'insetfit': True,
        'legendcols': ['f'],
        'legendoutside': True,
        'legendcollect': True,
        'legendlocation': 'center',
        'owr_pad': 1.5,  # Make the left-most subplot wider by this factor to make space for ylabel and ticklabels
        # 'legendcols': ['f', 'xi', 'beta'],
        # 'legend_in_subfig': False
    }
    meta_lbit_exp_fit_uniform = {
        'default': default,
        'figsize': (1.702, 1.702),
        'filename': "{}/lbit-exp-fit".format(plotdir),
        'datafile': 'data/mbl_10004-uniform50-L24-u8-f0.5.h5',
        'suptitle': None,  # 'l-bit fit $y=y_0 \exp[-(x/\\xi)^\\beta]$',
        'ylabel': '$\\bar O(i,j)$',
        'xlabel': '$|i-j|$',
        'yticks': [1e-0, 1e-4, 1e-8, 1e-12],
        'xticks': [0, 2, 4, 6, 8, 10, 12, 14, 16],
        'sharex': 'all',
        'sharey': 'all',
        'yscale': 'log',
        'xscale': 'linear',
        'insetfit': True,
        'legendcols': ['f'],
        'legendoutside': True,
        'legendcollect': True,
        'legendlocation': 'center',
        'owr_pad': 1.5,  # Make the left-most subplot wider by this factor to make space for ylabel and ticklabels
        # 'legendcols': ['f', 'xi', 'beta'],
        # 'legend_in_subfig': False
    }

    meta_lbits_expdecay = {
        'default': default,
        'figsize': (1.702, 1.702),
        'filename': "{}/lbits-expedecay".format(plotdir),
        # 'datafile': 'data/mbl_10004-constricted50-L24-u8-f0.5.h5',
        # 'datafile': 'data/mbl-expdecay50-L20-u16-f0_51.h5',# Pretty good
        # 'datafile': 'data/mbl-expdecay50-L20-u16-f0_53.h5',# Even better
        'datafile': 'data/mbl-expdecay20-L20-u16-f0.2.h5',
        'suptitle': None,  # 'l-bit fit $y=y_0 \exp[-(x/\\xi)^\\beta]$',
        'ylabel': '$|O(L/2,j)|$',
        'xlabel': '$j$',
        'yticklabels': [1e-0, 1e-4, 1e-8, 1e-12],
        'yticks': [1e-0, 1e-4, 1e-8, 1e-12],
        'xticks': [0, 5, 10, 15 ,20],
        # 'xticks': [0, 8, 16, 24],
        # 'xticks': [0, 6, 12, 18, 24],
        'ymaloc': LogLocator(base=10, numticks=4, numdecs=4),
        # 'xmiloc': LogLocator(base=10, numticks=10, subs=(.1, .2, .3, .4, .5, .6, .7, .8, .9)),
        # 'ymax': 2.38, # To make room for an axes label with a box
        'ymafmt': LogFormatterMathtext(labelOnlyBase=False),
        'box_aspect': 1,
        'sharex': 'all',
        'sharey': 'all',
        'yscale': 'log',
        'ynopos': 'mask',
        # 'box_aspect': 1,
        # 'xscale': 'linear',
        # 'ymin': 1e-20,
        'ymax': 2e0,
        'ymin': 1e-20,
        'xmin': -1,
        'xmax': 21,
        'fieldinset': True,
        'plotfit': True,
        'lbit-site': 10,
        'lbit-mean': 'arithmetic',
        'legendoutside': False,
        'legendcollect': False,
        'legendlocation': (0.25, 0.22),
        'owr_pad': 1.2,  # Make the left-most subplot wider by this factor to make space for ylabel and ticklabels
        'ohr_pad': 1.0,  # Make the left-most subplot wider by this factor to make space for ylabel and ticklabels
        # 'legendcols': ['f', 'u', 'xi', 'beta'],
        # 'legendcols': ['f', 'u'],
        'legendcols': ['wi', 'xi', 'beta'],
        'wilabel': 'on',
        # 'legend_in_subfig': False

        #     auto fields = std::vector<double>{
        #         -0.19585106363434682, // -0.19585106363434682, // -0.19585106363434682,
        #         0.9613590045929011,   // 0.9613590045929011,   // 0.9613590045929011,
        #         -0.521851938675176,   // -0.521851938675176,   // -0.521851938675176
        #         0.8162082991856783,   // 0.2162082991856783,   // 0.8162082991856783
        #         -0.01399945249037011, // -0.11399945249037011, // -0.01399945249037011
        #         1.8519104076563720,   // 1.8519104076563720,    // 1.8519104076563720
        #         1.8753286843228538,   // 1.8753286843228538,   // 1.8753286843228538
        #         2.1889052226219166,   // 1.1889052226219166,   // 2.1889052226219166
        #         -0.5914534639884793,  // 0.5914534639884793,   // -0.5914534639884793
        #         0.4715913075170785,   // 0.4715913075170785,   // 0.4715913075170785
        #         -0.889838217949121,   // -0.889838217949121,   // -0.889838217949121
        #         -1.3603722511496172,  // -1.3603722511496172,  // -1.3603722511496172
        #         -0.6781471711782956,  // -0.9781471711782956,  // -0.6781471711782956
        #         -0.7558612412400791,  // -0.7558612412400791,  // -0.7558612412400791
        #         0.9543219068987924,   // 0.9543219068987924,   // 0.9543219068987924
        #         0.8885578827760061,   // 0.5885578827760061,   // 0.3885578827760061
        #         1.1576139995533495,   // 1.9576139995533495,   // 1.9576139995533495
        #         0.7620903099433074,   // 0.7620903099433074,   // 0.7620903099433074
        #         -1.4662507329159788,  // -1.4662507329159788,  // -1.4662507329159788
        #         -0.26967264195989105, // -1.26967264195989105, // -0.26967264195989105
        #     };


    }
    meta_lbits_identity = {
        'default': default,
        'figsize': (1.702, 1.702),
        'filename': "{}/lbits-expedecay".format(plotdir),
        'datafile': 'data/mbl-identity20-L20-u16-f0.2.h5',
        'suptitle': None,  # 'l-bit fit $y=y_0 \exp[-(x/\\xi)^\\beta]$',
        'ylabel': '$|O(L/2,j)|$',
        'xlabel': '$j$',
        'yticklabels': [1e-0, 1e-4, 1e-8, 1e-12],
        'yticks': [1e-0, 1e-4, 1e-8, 1e-12],
        'xticks': [0, 5, 10, 15, 20],
        # 'xticks': [0, 8, 16, 24],
        # 'xticks': [0, 6, 12, 18, 24],
        'ymaloc': LogLocator(base=10, numticks=4, numdecs=4),
        # 'xmiloc': LogLocator(base=10, numticks=10, subs=(.1, .2, .3, .4, .5, .6, .7, .8, .9)),
        # 'ymax': 2.38, # To make room for an axes label with a box
        'ymafmt': LogFormatterMathtext(labelOnlyBase=False),
        'box_aspect': 1,
        'sharex': 'all',
        'sharey': 'all',
        'yscale': 'log',
        'ynopos': 'mask',
        # 'box_aspect': 1,
        # 'xscale': 'linear',
        # 'ymin': 1e-20,
        'ymax': 2e0,
        'ymin': 1e-20,
        'xmin': -1,
        'xmax': 21,
        'lbit-site': 10,
        'lbit-mean': 'arithmetic',
        'fieldinset': False,
        'plotfit': True,
        'legendoutside': False,
        'legendcollect': False,
        'legendlocation': (0.25, 0.22),
        'owr_pad': 1.2,  # Make the left-most subplot wider by this factor to make space for ylabel and ticklabels
        'ohr_pad': 1.0,  # Make the left-most subplot wider by this factor to make space for ylabel and ticklabels
        # 'legendcols': ['f', 'u', 'xi', 'beta'],
        'legendcols': ['wi'],
        'wilabel': 'off',
        # 'legend_in_subfig': False
    }


    meta_lbits_identity_old = {
        'default': default,
        'figsize': (1.702, 1.702),
        'filename': "{}/lbits-identity".format(plotdir),
        'datafile': 'data/mbl_10004-uniform50-L24-u6-f0.3.h5',
        'suptitle': None,  # 'l-bit fit $y=y_0 \exp[-(x/\\xi)^\\beta]$',
        'ylabel': '$\log_{10} O(\\frac{L}{2},j)$',
        'xlabel': '$j$',
        'yticks': [1e-0, 1e-4, 1e-8, 1e-12],
        'xticks': [0, 8, 16, 24],
        # 'xticks': [0, 6, 12, 18, 24],
        'box_aspect': 1,
        'sharex': 'all',
        'sharey': 'all',
        'yscale': 'log',
        'ynopos': 'mask',
        # 'box_aspect': 1,
        # 'xscale': 'linear',
        'ymin': 1e-12,
        'ymax': 2e0,
        # 'ymin': 1e-16,
        'xmin': -1,
        'xmax': 25,
        'fieldinset': False,
        'legendoutside': False,
        'legendcollect': False,
        'legendlocation': (0.54, 0.75),
        'owr_pad': 1.0,  # Make the left-most subplot wider by this factor to make space for ylabel and ticklabels
        'ohr_pad': 1.0,  # Make the left-most subplot wider by this factor to make space for ylabel and ticklabels
        # 'legendcols': ['f', 'u', 'xi', 'beta'],
        'legendcols': ['f', 'u'],
        # 'legend_in_subfig': False
    }

    meta_lbits_expdecay_halfheight = {
        'default': default_halfheight,
        # 'figsize': (1.702, 1.702*0.5),
        'filename': "{}/lbits-expedecay".format(plotdir),
        # 'datafile': 'data/mbl_10004-constricted50-L24-u8-f0.5.h5',
        # 'datafile': 'data/mbl-expdecay50-L20-u16-f0_51.h5',# Pretty good
        # 'datafile': 'data/mbl-expdecay50-L20-u16-f0_53.h5',# Even better
        'datafile': 'data/mbl-expdecay20-L20-u16-f0.2.h5',
        'suptitle': None,  # 'l-bit fit $y=y_0 \exp[-(x/\\xi)^\\beta]$',
        'xticks': [0, 20],
        'xlabel': '$j$',
        'xlabelpad': -8,
        'ylabel': '$O(\\frac{L}{2},j)$',
        'yticks': [1e-2, 1e-20],
        # 'yticklabels': ['$-2$', '$-8$', '$-14$', '$-20$'],
        'ylabelpad': -16,
        # 'xticks': [0, 8, 16, 24],
        # 'xticks': [0, 6, 12, 18, 24],
        # 'ymaloc': LogLocator(base=10, numticks=4, numdecs=4),
        # 'xmiloc': LogLocator(base=10, numticks=10, subs=(.1, .2, .3, .4, .5, .6, .7, .8, .9)),
        # 'ymax': 2.38, # To make room for an axes label with a box
        # 'ymafmt': LogFormatterMathtext(labelOnlyBase=False),
        'sharex': 'all',
        'sharey': 'all',
        'yscale': 'log',
        'ynopos': 'mask',
        # 'xscale': 'linear',
        # 'ymin': 1e-20,
        'ymax': 2e0,
        'ymin': 1e-20,
        'xmin': -1,
        'xmax': 21,
        'fieldinset': True,
        'plotfit': True,
        'lbit-site': 10,
        'lbit-mean': 'arithmetic',
        'fit-beta': True,
        'legendoutside': False,
        'legendcollect': False,
        'legendlocation': (0.44, 0.20),
        'owr_pad': 1.2,  # Make the left-most subplot wider by this factor to make space for ylabel and ticklabels
        'ohr_pad': 1.0,  # Make the left-most subplot wider by this factor to make space for ylabel and ticklabels
        # 'legendcols': ['f', 'u', 'xi', 'beta'],
        # 'legendcols': ['f', 'u'],
        'legendcols': ['wi'],
        'wilabel': 'on',
        # 'legend_in_subfig': False

        #     auto fields = std::vector<double>{
        #         -0.19585106363434682, // -0.19585106363434682, // -0.19585106363434682,
        #         0.9613590045929011,   // 0.9613590045929011,   // 0.9613590045929011,
        #         -0.521851938675176,   // -0.521851938675176,   // -0.521851938675176
        #         0.8162082991856783,   // 0.2162082991856783,   // 0.8162082991856783
        #         -0.01399945249037011, // -0.11399945249037011, // -0.01399945249037011
        #         1.8519104076563720,   // 1.8519104076563720,    // 1.8519104076563720
        #         1.8753286843228538,   // 1.8753286843228538,   // 1.8753286843228538
        #         2.1889052226219166,   // 1.1889052226219166,   // 2.1889052226219166
        #         -0.5914534639884793,  // 0.5914534639884793,   // -0.5914534639884793
        #         0.4715913075170785,   // 0.4715913075170785,   // 0.4715913075170785
        #         -0.889838217949121,   // -0.889838217949121,   // -0.889838217949121
        #         -1.3603722511496172,  // -1.3603722511496172,  // -1.3603722511496172
        #         -0.6781471711782956,  // -0.9781471711782956,  // -0.6781471711782956
        #         -0.7558612412400791,  // -0.7558612412400791,  // -0.7558612412400791
        #         0.9543219068987924,   // 0.9543219068987924,   // 0.9543219068987924
        #         0.8885578827760061,   // 0.5885578827760061,   // 0.3885578827760061
        #         1.1576139995533495,   // 1.9576139995533495,   // 1.9576139995533495
        #         0.7620903099433074,   // 0.7620903099433074,   // 0.7620903099433074
        #         -1.4662507329159788,  // -1.4662507329159788,  // -1.4662507329159788
        #         -0.26967264195989105, // -1.26967264195989105, // -0.26967264195989105
        #     };


    }
    meta_lbits_identity_halfheight = {
        'default': default_halfheight,
        'filename': "{}/lbits-expedecay".format(plotdir),
        'datafile': 'data/mbl-identity20-L20-u16-f0.2.h5',
        'suptitle': None,  # 'l-bit fit $y=y_0 \exp[-(x/\\xi)^\\beta]$',
        'xticks': [0, 19],
        'xlabel': '$j$',
        'xlabelpad': -8,
        'ylabel': '$O(\\frac{L}{2},j)$',
        'yticks': [1e-2, 1e-20],
        # 'yticklabels': ['$-2$', '$-8$', '$-14$', '$-20$'],
        'ylabelpad': -16,
        # 'xticks': [0, 8, 16, 24],
        # 'xticks': [0, 6, 12, 18, 24],
        # 'ymaloc': LogLocator(base=10, numticks=4, numdecs=4),
        # 'xmiloc': LogLocator(base=10, numticks=10, subs=(.1, .2, .3, .4, .5, .6, .7, .8, .9)),
        # 'ymax': 2.38, # To make room for an axes label with a box
        # 'ymafmt': LogFormatterMathtext(labelOnlyBase=False),
        'sharex': 'all',
        'sharey': 'all',
        'yscale': 'log',
        'ynopos': 'mask',
        # 'xscale': 'linear',
        # 'ymin': 1e-20,
        'ymax': 2e0,
        'ymin': 1e-20,
        'xmin': 0,
        'xmax': 20,
        'lbit-site': 10,
        'lbit-mean': 'arithmetic',
        'fit-beta':True,
        'fieldinset': False,
        'plotfit': True,
        'legendoutside': False,
        'legendcollect': False,
        'legendlocation': (0.44, 0.20),
        'owr_pad': 1.2,  # Make the left-most subplot wider by this factor to make space for ylabel and ticklabels
        'ohr_pad': 1.0,  # Make the left-most subplot wider by this factor to make space for ylabel and ticklabels
        # 'legendcols': ['f', 'u', 'xi', 'beta'],
        'legendcols': ['wi'],
        'wilabel': '$w_i=1$',
        # 'legend_in_subfig': False
    }

    with plt.style.context(mplstyle):
        # fig_exp_fit = plot_decay(meta_lbit_exp_fit_blocked)
        # fig_exp_fit = plot_decay(meta_lbit_exp_fit_uniform, figs=fig_exp_fit)
        # save_figure(fig_exp_fit)

        # fig_identity = plot_lbits_new(meta_lbits_identity, color='coral')
        # fig_expdecay = plot_lbits_new(meta=meta_lbits_expdecay, figs=fig_identity, color='forestgreen')

        fig_identity = plot_lbits_new(meta_lbits_identity_halfheight, color='gray')
        fig_expdecay = plot_lbits_new(meta=meta_lbits_expdecay_halfheight, figs=fig_identity, color='black')
        save_figure(fig_identity)
        save_figure(fig_expdecay)

        plt.show()

# def plot_lbits(figs=None):
#     meta = {
#         'filename': "{}/lbits".format(plotdir),
#         'suptitle': None,  # 'l-bit fit $y=y_0 \exp[-(x/\\xi)^\\beta]$',
#         'ylabel': '$O(L/2,j)$',
#         'xlabel': '$j$',
#         # 'yticks': [1e-0, 1e-4, 1e-8, 1e-12],
#         # 'xticks': [0, 8, 16, 24],
#         'xticks': [0, 6, 12, 18, 24],
#         'sharex': 'all',
#         'sharey': 'all',
#         'yscale': 'log',
#         'ynopos': 'mask',
#         # 'box_aspect': 1,
#         # 'xscale': 'linear',
#         'ymin': 1e-10,
#         'ymax': 2e0,
#         # 'ymin': 1e-16,
#         'xmin': -1,
#         'xmax': 25,
#         'legendcols': ['f'],
#         'legendoutside': False,
#         'legendcollect': False,
#         'legendlocation': 'upper right',
#         'owr_pad': 1.0,  # Make the left-most subplot wider by this factor to make space for ylabel and ticklabels
#         'ohr_pad': 1.0,  # Make the left-most subplot wider by this factor to make space for ylabel and ticklabels
#         # 'legendcols': ['f', 'xi', 'beta'],
#         # 'legend_in_subfig': False
#     }
#     if 'slack' in meta['filename']:
#         palette_name = "Spectral_r"
#         path_effects = [pe.SimpleLineShadow(offset=(-0.3, -0.5), alpha=0.4), pe.Normal()]
#     else:
#         palette_name = "colorblind"
#         path_effects = None
#
#     f = None
#     # sns.set_palette(sns.color_palette(palette_name, n_colors=reps))
#     print('Starting plots colors')
#     print('Defining colors')
#     current_palette = itertools.cycle(sns.color_palette())
#     print(current_palette)
#     for datafile in datafiles:
#         h5data = h5open(datafile, 'r')
#         color = next(current_palette)
#         frange = h5data["{}".format(lbitpath)].attrs["u_fmix"][()]
#         urange = h5data["{}".format(lbitpath)].attrs["u_depth"][()]
#         lbits = h5data["{}".format(lbitpath)]["data"][()]
#         print('frange: {}'.format(frange))
#         print('urange: {}'.format(urange))
#
#         # lbits = lbits[-7::7, -1::1, :, :, :]
#         # frange = frange[-7::7]
#         # urange = urange[-1::1]
#         # lbits = lbits[0::4, 0::1, :, :, :]
#         # frange = frange[0::4]
#         # urange = urange[::1]
#
#         print('frange: {}'.format(frange))
#         print('urange: {}'.format(urange))
#
#         flen, ulen, reps, rows, cols = np.shape(lbits)
#         numsubs = flen * ulen
#         if figs is None:
#             figs = [get_fig_meta(numsubs, meta=meta)]
#         f = figs[-1]
#         fields = h5data["fLBIT/model/hamiltonian"]['J1_rand']
#         print(fields)
#         # path_effects = [pe.SimpleLineShadow(offset=(0.5, -0.5), alpha=0.3), pe.Normal()]
#         path_effects = [pe.Stroke(linewidth=2, foreground='white'), pe.Normal()]
#
#         for fidx, fval in enumerate(frange):
#             for uidx, uval in enumerate(urange):
#                 idx = np.ravel_multi_index([[fidx], [uidx]], dims=(flen, ulen), )[0]
#                 print('idx {} | fidx {} | uidx {}'.format(idx, fidx, uidx))
#                 ax = f['ax'][idx]
#                 maxreps = np.min([reps, 50])
#
#                 for site in [int(rows / 2)]:
#                     x = range(cols)
#                     l = lbits[fidx, uidx, range(maxreps), [site], :].T
#                     y = np.mean(l, axis=1)
#                     e = np.std(l, axis=1)
#                     line, = ax.plot(x, y, linewidth=1.0, color=color, alpha=0.5)
#                     # ax.xaxis.set_major_locator(MaxNLocator(integer=True))
#                     # f['legends'][idx][0]['handle'].append(line)
#                     # f['legends'][idx][0]['label'].append('{}'.format(cols))
#                     # f['legends'][idx][0]['title'] = '$L$'
#                     f['legends'][idx][1]['handle'].append(line)
#                     f['legends'][idx][1]['label'].append('{:.2f}'.format(fval))
#                     f['legends'][idx][1]['title'] = '$f$'
#                     f['legends'][idx][2]['handle'].append(line)
#                     f['legends'][idx][2]['label'].append('{}'.format(uval))
#                     f['legends'][idx][2]['title'] = '$d_u$'
#                     for ridx in range(maxreps):
#                         y = lbits[fidx, uidx, ridx, [site], :].T
#                         y[y < 1e-15] = np.nan
#                         ax.plot(x, y, alpha=0.30, color=color)
#
#                         if not idx in f['axes_used']:
#                             f['axes_used'].append(idx)
#
#                 # trans = transforms.blended_transform_factory(ax.transData, ax.transAxes)
#                 # fieldcolor = 'steelblue'
#                 # for site in range(rows):
#                 #     ax.arrow(site,0.12,0.0,fields[site]/15.0, width=0.025, length_includes_head=False,
#                 #              head_width=0.00, head_length=0.00,
#                 #              path_effects=path_effects, zorder=50,
#                 #              transform=trans,
#                 #              color=fieldcolor,
#                 #              #color='salmon'
#                 #              )
#                 # ax.arrow(0, 0.12, rows, 0, width=0.001, length_includes_head=False,
#                 #          head_width=0.00, head_length=0.00,
#                 #          zorder=60,transform=trans,
#                 #          color=fieldcolor
#                 #          #color='salmon'
#                 #          )
#
#
#     return figs
