import seaborn as sns
from src.plotting.tools import *
import matplotlib.pyplot as plt
from src.io.h5ops import *
from scipy.optimize import curve_fit
import matplotlib.patheffects as pe
from pathlib import Path
import logging
from numba import njit
from matplotlib.ticker import MaxNLocator


@njit(parallel=True, cache=True)
def stretched_exp(x, C, xi, beta):
    return C * np.exp(-(x / xi) ** beta)


@njit(parallel=True, cache=True)
def stretched_log(x, C, xi, beta):
    return np.log(C) - (x / xi) ** beta


plotdir = 'plots/lbit-decay'
datafiles = [
    '../../output/mbl_10003-uniform500.h5',
    '../../output/mbl_10003-normal500.h5',
    '../../output/mbl_10003-squared500.h5',
]

lbitpath = 'fLBIT/model/lbits'


def plot_decay():
    logging.basicConfig(level=logging.WARN)
    logger_lbit = logging.getLogger('lbit').setLevel(level=logging.WARN)
    logger_tools = logging.getLogger('tools').setLevel(level=logging.INFO)


    meta = {
        'suptitle': None,  # 'l-bit fit $y=y_0 \exp[-(x/\\xi)^\\beta]$',
        'ylabel': '$\\bar O(i,j)$',
        'xlabel': '$|i-j|$',
        'yticks': [1e-0, 1e-4, 1e-8, 1e-12],
        'xticks': [0, 2, 4, 6, 8],
        'sharex': 'all',
        'sharey': 'all',
        'yscale': 'log',
        'xscale': 'linear',
        'plotprefix': 'lbit-decay',
        'plotdir': plotdir,
        'mplstyle': './src/plotting/stylesheets/prb.mplstyle',
        # 'mplstyle': './src/plotting/stylesheets/slack.mplstyle',
        'legendcols': ['f'],
        'legendoutside': True,
        'legendcollect': True,
        'legendlocation': 'center',
        'owr_pad': 1.5,  # Make the left-most subplot wider by this factor to make space for ylabel and ticklabels
        # 'legendcols': ['f', 'xi', 'beta'],
        # 'legend_in_subfig': False
    }


    if 'mplstyle' in meta:
        plt.style.use(meta['mplstyle'])
    if 'plotdir' in meta and 'mplstyle' in meta:
        if Path(meta['plotdir']).stem != Path(meta['mplstyle']).stem:
            meta['plotdir'] = Path(meta['plotdir'], Path(meta['mplstyle']).stem)
            Path(meta['plotdir']).mkdir(parents=True, exist_ok=True)
            print("Setting plotdir: ", meta['plotdir'])
    if 'slack' in meta.get('mplstyle'):
        palette_name = "Spectral_r"
        path_effects = [pe.SimpleLineShadow(offset=(-0.3, -0.5), alpha=0.4), pe.Normal()]
    else:
        palette_name = "colorblind"
        path_effects = None
    meta['filename'] = "{}/{}_liom-decay".format(meta['plotdir'], meta['plotprefix'])
    f = None
    for datafile in datafiles:
        h5data = h5open(datafile, 'r')

        frange = h5data["{}".format(lbitpath)].attrs["f_mixer"][()]
        urange = h5data["{}".format(lbitpath)].attrs["u_layer"][()]  # formerly u_depth
        cls_avg = h5data["{}/cls_avg".format(lbitpath)][()]
        sse_avg = h5data["{}/sse_avg".format(lbitpath)][()]
        decay = h5data["{}/decay".format(lbitpath)]  # formerly "curves"

        print("decay shape", np.shape(decay))
        rows, cols = np.shape(cls_avg)
        cls_fit = np.zeros((rows, cols))
        cls_std = np.zeros((rows, cols))
        cls_1pc = np.zeros((rows, cols))
        cls_bet = np.zeros((rows, cols))
        print(rows, cols)
        fmix_range = range(0, len(frange), 2)
        ucol_range = range(0, len(urange), 2)
        numplots = len(ucol_range)
        print(numplots)
        if f is None:
            f = get_fig_meta(numplots, meta=meta)
        sns.set_palette(sns.color_palette(palette_name, n_colors=len(fmix_range)))
        figin = None
        fits = {'u': [], 'color': [], 'f': [], 'y0': [], 'xi': [], 'beta': []}
        for idx, (col, ax) in enumerate(zip(ucol_range, f['ax'])):
            current_palette = itertools.cycle(sns.color_palette())
            axin = None
            fit = {'u': urange[col], 'color': [], 'f': [], 'y0': [], 'xi': [], 'beta': []}
            for i, row in enumerate(fmix_range):
                color = next(current_palette)
                y = decay[row, col, :]
                y = y[y > 1e-15]
                x = np.array(range(len(y)))
                p0 = 0.9, 0.4, 1.0
                try:
                    ylog = np.log(y)
                    popt, pcov = curve_fit(stretched_log, xdata=x, ydata=ylog, p0=p0)
                    pstd = np.sqrt(np.diag(pcov))
                    print(popt)

                    # bounds = ([-1, -1, -1], [1, np.inf, np.inf])
                    # popt,pcov = curve_fit(stretched_exp, xdata=x, ydata=y, p0=p0)
                    # pstd = np.sqrt(np.diag(pcov))
                    # print(popt)
                except:
                    popt = [1, 0, 1]
                    pstd = [0, 0, 0]
                fit['color'].append(color)
                fit['f'].append(frange[row])
                fit['y0'].append(popt[0])
                fit['xi'].append(popt[1])
                fit['beta'].append(popt[2])
                cls_fit[row, col] = popt[1]
                cls_std[row, col] = pstd[1]
                cls_bet[row, col] = popt[2]
                cls_1pc[row, col] = popt[1] * (np.log(popt[0] / (1.0 / 100.0)) ** (1.0 / popt[2]))
                line, = ax.plot(x, y, marker='o', linestyle='None', color=color, path_effects=path_effects)
                ax.plot(x, stretched_exp(x, popt[0], popt[1], popt[2]), label=None, color=color, path_effects=path_effects)
                ax.set_title("$d_u={}$".format(urange[col]), bbox={'facecolor': 'white', 'pad': 2})
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
    return f


def plot_lbits():
    meta = {
        'suptitle': None,  # 'l-bit fit $y=y_0 \exp[-(x/\\xi)^\\beta]$',
        'ylabel': '$O(L/2,j)$',
        'xlabel': '$j$',
        # 'yticks': [1e-0, 1e-4, 1e-8, 1e-12],
        # 'xticks': [0, 2, 4, 6, 8],
        'sharex': 'all',
        'sharey': 'all',
        'yscale': 'log',
        'ynopos': 'mask',
        # 'box_aspect': 1,
        # 'xscale': 'linear',
        'ymin': 1e-12,
        'plotprefix': 'lbit-decay',
        'plotdir': plotdir,
        # 'mplstyle': './src/plotting/stylesheets/prb.mplstyle',
        'mplstyle': './src/plotting/stylesheets/slack.mplstyle',
        'legendcols': ['f'],
        'legendoutside': False,
        'legendcollect': False,
        'legendlocation': 'upper right',
        'owr_pad': 1.0,  # Make the left-most subplot wider by this factor to make space for ylabel and ticklabels
        'ohr_pad': 1.0,  # Make the left-most subplot wider by this factor to make space for ylabel and ticklabels
        # 'legendcols': ['f', 'xi', 'beta'],
        # 'legend_in_subfig': False
    }

    if 'mplstyle' in meta:
        plt.style.use(meta['mplstyle'])
    if 'plotdir' in meta and 'mplstyle' in meta:
        if Path(meta['plotdir']).stem != Path(meta['mplstyle']).stem:
            meta['plotdir'] = Path(meta['plotdir'], Path(meta['mplstyle']).stem)
            Path(meta['plotdir']).mkdir(parents=True, exist_ok=True)
            print("Setting plotdir: ", meta['plotdir'])
    if 'slack' in meta.get('mplstyle'):
        palette_name = "Spectral_r"
        path_effects = [pe.SimpleLineShadow(offset=(-0.3, -0.5), alpha=0.4), pe.Normal()]
    else:
        palette_name = "colorblind"
        path_effects = None
    meta['filename'] = "{}/{}_lioms".format(meta['plotdir'], meta['plotprefix'])
    f = None
    current_palette = itertools.cycle(sns.color_palette())
    for datafile, gcolor in zip(datafiles, ['lightsteelblue', 'limegreen', 'gray']):
        h5data = h5open(datafile, 'r')
        color = next(current_palette)
        frange = h5data["{}".format(lbitpath)].attrs["f_mixer"][()]
        urange = h5data["{}".format(lbitpath)].attrs["u_layer"][()]  # formerly u_depth
        lbits = h5data["{}".format(lbitpath)]["data"][()]

        lbits = lbits[1::4, 1::1, :, :, :]
        frange = frange[1::4]
        urange = urange[1::1]
        print('frange: {}'.format(frange))
        print('urange: {}'.format(urange))

        flen, ulen, reps, rows, cols = np.shape(lbits)
        numplots = flen * ulen

        if f is None:
            print('Defining f')
            f = get_fig_meta(numplots, meta=meta)
        print('Defining colors')
        sns.set_palette(sns.color_palette(palette_name, n_colors=reps))
        print('Starting plots colors')
        for fidx, fval in enumerate(frange):
            for uidx, uval in enumerate(urange):
                idx = ((fidx + 1) * (uidx + 1)) - 1
                ax = f['ax'][idx]
                mid = int(rows / 2)
                x = range(cols)
                l = lbits[fidx, uidx, range(reps), [mid], :].T
                y = np.mean(l, axis=1)
                e = np.std(l, axis=1)
                line, = ax.plot(x, y, linewidth=2.0, color=color)
                ax.xaxis.set_major_locator(MaxNLocator(integer=True))
                f['legends'][idx][0]['handle'].append(line)
                f['legends'][idx][0]['label'].append('{:.2f}'.format(fval))
                f['legends'][idx][0]['title'] = '$f$'

                for ridx in range(reps):
                    y = lbits[fidx, uidx, ridx, [mid], :].T
                    ax.plot(x, y, alpha=0.05, color=color)

                    if not idx in f['axes_used']:
                        f['axes_used'].append(idx)

    return f


if __name__ == '__main__':
    # f = plot_decay()
    # save_figure(f)
    f = plot_lbits()
    save_figure(f)
    plt.show()
