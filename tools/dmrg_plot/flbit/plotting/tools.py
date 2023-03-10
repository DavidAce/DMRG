import os

import matplotlib.pyplot
import numpy as np
from numba import njit
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.ticker import FormatStrFormatter
from matplotlib import rcParams, legend
from copy import deepcopy
import itertools
import pkg_resources
from packaging import version
import logging
import matplotlib.gridspec as gs
from git import Repo
import seaborn as sns
from itertools import product
from scipy.optimize import curve_fit
from scipy.stats import linregress
from collections.abc import Iterable

logger = logging.getLogger('tools')
import tikzplotlib
from pathlib import Path

matplotlib_version = pkg_resources.get_distribution("matplotlib").version
if version.parse(matplotlib_version) == version.parse("3.5"):
    logger.warning("Matplotlib version {} may not work correctly with gridspec".format(matplotlib_version))


def get_markerlist():
    return itertools.cycle(('^', 'v', '<', '>', 'o'))


def get_linestyles():
    return itertools.cycle(("-", "--", "-.", ":"))


def get_uniform_palette_names(num):
    if num > 6:
        raise ValueError("num == {} must be smaller than 7".format(num))
    palettes = [
        'Blues',  # sequential from matplotlib
        'Greens',  # sequential from matplotlib
        'Oranges',  # sequential from matplotlib
        'Purples',  # sequential from matplotlib
        'Reds',  # sequential from matplotlib
        'crest',  # perceptually uniform from seaborn
        'flare',  # perceptually uniform from seaborn
        'magma_r',  # perceptually uniform from matplotlib (reversed)
        'viridis_r',  # perceptually uniform from matplotlib (reversed)
        'dark:salmon_r',
    ]

    return palettes[0:num]


def get_colored_lstyles(db, linspec, default_palette):
    linprod = list(
        product(*get_vals(db, linspec)))  # All combinations of linspecs (names of parameters that iterate lines)
    palette = sns.color_palette(palette=default_palette, n_colors=len(linprod))
    lstyles = [None] * len(linprod)
    if len(linspec) == 2:
        linkey0 = get_keys(db, linspec[0])  # Sets number of colors
        linkey1 = get_keys(db, linspec[1])  # Sets number of linestyles
        if len(linkey0) == 2:
            palette = reversed(sns.color_palette(palette='tab20', n_colors=len(linprod)))
        if len(linkey0) == 3:
            palette = reversed(sns.color_palette(palette='tab20c', n_colors=len(linprod)))
            # del palette[4 - 1::4]
        if len(linkey0) == 4:
            palette = reversed(sns.color_palette(palette='tab20c', n_colors=len(linprod)))

    return palette, lstyles


def match_datanodes(db, meta, specs, vals):
    nodes = set()
    gname = meta.get('groupname')
    dname = meta.get('dsetname')
    for dsetpath, dset in db['dsets'].items():
        # print('-- checking', dsetpath)
        if gname is not None:
            if isinstance(gname, Iterable) and not any([g in dsetpath for g in gname]):
                continue
            elif not gname in dsetpath:
                continue
        if dname is not None:
            if isinstance(dname, Iterable) and not any([d in dsetpath for d in dname]):
                continue
            elif not dname in dsetpath:
                continue
        equal = True
        for s, v in zip(specs, vals):
            s_noformat = s.split(':')[0]
            if isinstance(s, Iterable) and isinstance(v, Iterable):
                equal = all([x == y for x, y in zip(dset['vals'][s_noformat], v)])
            else:
                equal = dset['vals'][s_noformat] == v
            if not equal:
                break
        if not equal:
            continue
        nodes.add(dset['node']['data'])
    return nodes


def write_attributes(*args, **kwargs):
    for a in args:
        print(a)
    for k, v in kwargs.items():
        print("%s = %s" % (k, v))


def get_optimal_subplot_num(numplots):
    r = np.sqrt(numplots)
    cols = int(np.ceil(r))
    rows = int(np.floor(r))
    while cols * rows < numplots:
        if (cols <= rows):
            cols = cols + 1
        else:
            rows = rows + 1
    return rows, cols


def get_empty_figs_w_subplots(num_figs, num_subfigs, figsize=3.5):
    figures = []
    axes = []
    rows, cols = get_optimal_subplot_num(num_subfigs)
    for i in range(num_figs):
        fig, ax = plt.subplots(nrows=rows, ncols=cols, figsize=(figsize * cols, figsize * rows), num=i)
        figures.append(fig)
        axes.append(ax)
    return figures, axes


def remove_empty_subplots(fig, axes, axes_used):
    for idx, ax in enumerate(np.ravel(axes)):  # Last one is for the legend
        if not idx in axes_used:
            fig.delaxes(ax)


@njit(parallel=True, cache=True)
def find_saturation_idx(ydata, std_threshold):
    # Consider Y vs X: a noisy signal decaying in the shape of a hockey-club, say.
    # We want to identify the point at which the signal stabilizes. We use the fact that the
    # standard deviation is high if it includes parts of the non-stable signal, and low if
    # it includes only the stable part.
    # Here we monitor the standard deviation of the signal between [start_point, end_point],
    # and move "start_point" towards the end. If the standard deviation goes below a certain
    # threshold, i.e. threshold < max_std, then we have found the stabilization point.
    sdata = np.empty_like(ydata)
    for i in range(len(ydata)):
        sdata[i] = np.std(ydata[i:])
    sdiff = -np.log(np.abs(np.diff(sdata)))
    return np.argmax(sdiff)


#
# def find_saturation_idx(ydata, std_threshold):
#     sdata = []
#     for i in range(len(ydata)):
#         sdata.append(np.std(ydata[i:]))
#     sdiff = -np.log(np.abs(np.diff(sdata)))
#     return np.argmax(sdiff)


@njit(parallel=True, cache=True)
def stretched_exp(x, C, xi, beta):
    return C * np.exp(-(x / xi) ** beta)


@njit(parallel=True, cache=True)
def stretched_log(x, C, xi, beta):
    return np.log(C) - (x / xi) ** beta


@njit(parallel=True, cache=True)
def regular_exp(x, C, xi):
    return C * np.exp(-(x / xi))


@njit(parallel=True, cache=True)
def linear_log(x, C, xi):
    return np.log(C) - (x / xi)


@njit(parallel=True, cache=True)
def linear_fit(x, m, k):
    return m + k * x


# def floglog(x, a, b,c):
#     return a + b*np.log(np.log(x + c))
@njit(parallel=True, cache=True)
def floglog_v2(x, a, b, c):
    with np.errstate(invalid='ignore'):
        return a + b * np.log(np.log(x - c))


@njit(parallel=True, cache=True)
def fpower(x, a, b):
    with np.errstate(invalid='ignore'):
        return a * x ** b


# @njit(parallel=True, cache=True)
def flinear(x, a, b):
    with np.errstate(invalid='ignore'):
        return a + b * x


# Returns C, xi, beta, yfit, LinregressResult instance or pstd
def get_lbit_cls(x, y, stretched=False, ymin=1e-5):
    if np.size(y) <= 1:
        print('get_lbit_cls: y is too short:', y)
        return None, None, None, None, None, None
    ydata = np.ndarray.flatten(y)
    xdata = np.ndarray.flatten(x)
    idx0 = 1  # Skip first entry, which usually doesn't obey exp decay
    idxN = 0
    for idx, val in enumerate(ydata):
        if abs(val) < np.max(np.abs(ydata[idx:])):
            break
        if idx >= 2 and abs(val) > 0.5 * (abs(ydata[idx - 1]) + abs(ydata[idx - 2])):
            break
        if abs(val) < ymin:
            break
        if val == 0:
            break

        idxN = idx

    ydata = ydata[idx0:idxN]
    xdata = xdata[idx0:idxN]
    try:
        if idxN >= 0 and idxN <= idx0:
            raise IndexError("Invalid index order: idx0 {} | idxN {}".format(idx0, idxN))
        if idxN >= 0 and idxN < idx0 + 2:
            raise IndexError("Too few datapoints for a fit: idx0 {} | idxN {}".format(idx0, idxN))

        with np.errstate(invalid='ignore'):
            ylogs = np.log(np.abs(ydata))
            if stretched:
                p0 = 0.5, 1.0, 1.0
                popt, pcov = curve_fit(stretched_log, xdata=xdata, ydata=ylogs, p0=p0)
                pstd = np.sqrt(np.diag(pcov))
                return popt[0], popt[1], popt[2], stretched_exp(xdata, *popt), pstd, idxN
            else:
                result = linregress(x=xdata, y=ylogs, alternative='less')
                C = np.exp(result.intercept)
                xi = 1.0 / abs(result.slope)
                yfit = C * np.exp(-xdata / xi)
                return C, xi, None, yfit, result, idxN
    except IndexError as e:
        print("Index error:", e)
        pass
    except ValueError as e:
        print("Fit failed:", e)

    return None, None, None, None, None, None

def find_saturation_idx2(ydata, threshold=1e-2):
    if len(ydata) <= 2:
        return
    ylog = -np.log10(ydata)
    ylog = ylog / ylog[-1]
    sdata = []
    w = 2
    for i, yl in enumerate(ylog):
        min_idx = np.min([len(ylog) - w, i])
        min_idx = np.max([min_idx, 0])
        s = np.std(ylog[min_idx:])
        sdata.append(s)
    idx = np.argwhere(np.asarray(sdata) < threshold)[0, 0]
    return idx


def find_loglog_window(tdata, ydata, J2_width, J2_ctof, threshold1=0.6, threshold2=1e-2):
    if len(ydata) <= 2:
        return
    ylog = -np.log10(ydata)
    ylog = ylog / ylog[-1]
    sdata = []
    w = 2
    for i, yl in enumerate(ylog):
        min_idx = np.min([len(ylog) - w, i])
        min_idx = np.max([min_idx, 0])
        s = np.std(ylog[min_idx:])
        sdata.append(s)
    tdx = np.argwhere(np.asarray(tdata) >= 1.0)
    idx1 = np.argwhere(np.asarray(sdata) < threshold1)
    idx2 = np.argwhere(np.asarray(sdata) < threshold2)

    tdx = tdx[0, 0] if len(tdx) > 0 else 0
    idx1 = idx1[0, 0] if len(idx1) > 0 else 0
    idx2 = idx2[0, 0] if len(idx2) > 0 else 0

    idx1 = np.max([tdx, idx1])
    return idx1, idx2


def find_loglog_window2(tdata, ydata, db, threshold2=1e-2):
    if len(tdata) == 1:
        return 0, 0
    L = db['vals']['L']
    r = db['vals']['r']
    x = db['vals']['x']
    wn = db['vals']['w']

    w1 = wn[0]  # The width of distribution for on-site field.
    w2 = wn[
        1]  # The width of distribution for pairwise interactions. The distribution is either U(J2_mean-w,J2_mean+w) or N(J2_mean,w)
    w3 = wn[2]  # The width of distribution for three-body interactions.

    if r == np.iinfo(np.uint64).max:
        r = L

    # We multiply by np.sqrt(2 / np.pi) because we are interested in <|N(0,w)|>
    w1 *= np.sqrt(1 - 2 / np.pi)
    w2 *= np.sqrt(1 - 2 / np.pi)
    w3 *= np.sqrt(1 - 2 / np.pi)

    J1min = w1 * np.exp(0 / x)
    J1max = w1 * np.exp(0 / x)
    tmin1 = 1.0 / J1max
    tmax1 = 1.0 / J1min

    r2max = np.min([r, int((L - 1))])  # Maximum interaction range, max(|i-j|)
    Jmin2 = w2 * np.exp(
        - r2max / (2 * x))  # Order of magnitude of the smallest 2-body terms (furthest neighbor, up to L/2)
    Jmax2 = w2 * np.exp(- 1 / x)  # Order of magnitude of the largest 2-body terms (nearest neighbor)
    tmax2 = 1.0 / Jmin2  # (0.5 to improve fits) Time that it takes for the most remote site to interact with the middle
    tmin2 = 1.0 / Jmax2  # Time that it takes for neighbors to interact

    J3min = w3 * np.exp(-2 / x)
    J3max = w3 * np.exp(-2 / x)
    tmin3 = 1.0 / J3max
    tmax3 = 1.0 / J3min

    tmin = np.max([tmin1, tmin2, tmin3])
    tmax = np.max([tmax1, tmax2, tmax3])
    # Make sure the time points are ordered
    tmin = np.min([tmin, tmax])
    tmax = np.max([tmin, tmax])

    idx1 = np.argmax(tdata >= tmin) - 1
    idx2 = np.argmax(tdata >= tmax) - 1

    # idx1 = np.where(tdata >= tmin)[0][0]
    # idx2 = np.where(tdata <= tmax)[0][-1]
    # print(tmin1, tmin2, tmin3)
    # print(tmax1, tmax2, tmax3)
    # print(tdata[idx1], tdata[idx2])
    if idx2 == len(tdata) - 1:
        idx2 = np.max([idx1, idx2 - 1])
    return idx1, idx2


def find_saturation_idx3(tdata, ydata, db, threshold2=1e-2):
    if len(tdata) == 1:
        return 0, 0
    L = db['vals']['L']
    r = db['vals']['r']
    x = db['vals']['x']
    wn = db['vals']['w']

    w1 = wn[0]  # The width of distribution for on-site field.
    w2 = wn[1]  # The width of distribution for pairwise interactions. The distribution is either U(J2_mean-w,J2_mean+w) or N(J2_mean,w)
    w3 = wn[2]  # The width of distribution for three-body interactions.

    if r == np.iinfo(np.uint64).max:
        r = L

    tmax1 = 1.0 / w1

    r2max = np.min([r, L])  # Number of sites from the center site to the edge site, max(|i-j|)/2
    Jmin2 = np.exp(-(r2max - 1) / x) * w2 * 2 * np.sqrt(2 / np.pi)  # Order of magnitude of the smallest 2-body terms (furthest neighbor, up to L/2)
    Jmax2 = np.exp(-(1 - 1) / x) * w2  # Order of magnitude of the largest 2-body terms (nearest neighbor)
    tmax2 = 1.0 / Jmin2  # (0.5 to improve fits) Time that it takes for the most remote site to interact with the middle

    tmax3 = 1.0 / w3

    tmax = np.max([tmax1, tmax2, tmax3])
    idx = np.where(tdata <= tmax)[0][-1]
    idx = np.min([idx, len(tdata) - 10])  # Retain at least 10 points
    idx = np.max([idx, 0])  # Make sure its non-negative
    return idx


@njit(parallel=True, cache=True)
def page_entropy(L):
    n = int(2 ** (L / 2))
    S = - float(n - 1) / (2 * n)
    for k in range(1 + n, 1 + n ** 2):  # Include last
        S = S + 1.0 / k
    return S


def get_prop(db, keyfmt, prop):
    if isinstance(keyfmt, list) or isinstance(keyfmt, tuple):
        keys = []
        for kf in keyfmt:
            k = kf.split(':')[0]
            if prop in db:
                if k in db[prop]:
                    keys.append(db[prop][k])
                elif k in db[prop]['keys']:
                    keys.append(db[prop]['keys'][k])
        return keys
    else:
        if prop in db:
            k = keyfmt.split(':')[0]
            if k in db[prop]:
                return db[prop][k]
            elif k in db[prop]['keys']:
                return db[prop]['keys'][k]
        return []


def get_vals(db, keyfmt):
    return get_prop(db, keyfmt, 'vals')


def get_keys(db, keyfmt):
    return get_prop(db, keyfmt, 'keys')


def get_tex(db, keyfmt):
    keys = get_prop(db, keyfmt, 'tex')
    if isinstance(keys, list):
        keys = '$' + ', '.join([k.strip('$') for k in keys]) + '$'
    return keys


def get_title(db, keys, width=20):
    fmtvals = []
    newline = 0
    for s in keys:
        key, fmt = s.split(':') if ':' in s else [s, '']
        if key in db['vals'] and key in db['tex']['keys']:
            if isinstance(db['vals'][key], list):
                fmtvals.append('${}=[{}]$'.format(
                    db['tex']['keys'][key].strip('$'),
                    ','.join('{0:{1}}'.format(val, fmt) for val in db['vals'][key])))
            else:
                if key == 'r' and db['vals'][key] == np.iinfo(np.uint64).max:
                    fmtvals.append(db['tex']['eqs'][key])
                else:
                    fmtvals.append(
                        '${}={}$'.format(db['tex']['keys'][key].strip('$'), '{0:{1}}'.format(db['vals'][key], fmt)))
        elif key in db['tex']['eqs']:
            fmtvals.append(db['tex']['eqs'][key])
        if len(fmtvals) - newline >= width:
            newline = len(fmtvals)
            fmtvals.append('\n')
    return ', '.join(fmtvals)


def get_safe_table_field_name(node, keys):
    # Returns all keys that exist in node
    match = []
    if isinstance(keys, list):
        match = [field for field in keys if field in node.dtype.fields.keys()]
    else:
        if keys in node.dtype.fields.keys():
            match.append(keys)
    if not match:
        raise LookupError('Field not found in table:\n'
                          'Table   : {}\n'
                          'Expected: {}\n'
                          'Existing: {}'.format(node.name, keys, node.dtype.fields.keys()))
    return match


def get_table_data(node, keys=None, dtype='f8', transpose=None):
    if keys == None:
        return node[()].ravel().view((dtype, (1,)))[()], np.array([None], ndmin=2)

    elif node.dtype.fields == None:
        # This is a scalar put inside an ndarray, so that we have
        # a consistent return type
        return np.array(node[()]), np.array([None], ndmin=2)
    else:
        # Find the column names that match keys
        cols = get_safe_table_field_name(node, keys)
        return node.fields(cols)[()].view(dtype=(dtype, (len(cols),)))[()], np.array(cols)


def get_legend_row(db, datanode, legend_col_keys):
    legendrow = [None] * len(legend_col_keys)
    dbval = db['dsets'][datanode.name]
    for idx, key in enumerate(legend_col_keys):
        key, fmt = key.split(':') if ':' in key else [key, '']
        sfx = 'm' if 'time' in key or 'tsim' in key else ''
        if key in dbval['vals'] and fmt:
            legendrow[idx] = '{0:{1}}{2}'.format(dbval['vals'][key], fmt, sfx)
        elif key in dbval['tex']['vals']:
            legendrow[idx] = dbval['tex']['vals'][key]

        if legendrow[idx] is None:
            print(legend_col_keys)
            print(legendrow)
            print(dbval['tex'])
            raise
    return legendrow


def get_fig_meta(numplots: int, meta: dict):
    f = {
        'fig': None,
        'filename': meta.get('filename'),
        'lstyles': get_linestyles(),
        'constrained_layout': meta.get('constrained_layout'),
        'numplots': numplots,
        'nrows': None,
        'ncols': None,
        'figsize': meta.get('figsize'),
        'box_aspect': meta.get('box_aspect'),
        'crows': 1,  # Common row put at the bottom of all subplots
        'ccols': 1,  # Common col put at the right of all subplots
        'irows': 2,  # We make a 2x2 grid where 0,0 is the plot, and 0,1 and 1,0 are legends
        'icols': 2,  # We make a 2x2 grid where 0,0 is the plot, and 0,1 and 1,0 are legends

        'owr': None,
        'ohr': None,
        'iwr': [10000, 1],  # Width ratio between plot and right legend
        'ihr': [10000, 1],  # Height ratio between plot and bottom legend
        'owr_pad': meta.get('owr_pad') if 'owr_pad' in meta else 1.00,
        'ohr_pad': meta.get('ohr_pad') if 'ohr_pad' in meta else 1.00,

        'go': None,  # Outer gridspec
        'gi': [],  # Inner gridspec
        'rc': [],  # List of coordinates "(row,col)" for each ax
        'ax': [],  # List of subplots with plots (i.e. [0,0] in gsi)
        'ix': [],  # List of insets, one for each subplot axis
        'lr': [],  # List of subplots with legend right (i.e. [0,1] in gsi)
        'lb': [],  # List of subplots with legend below (i.e. [1,0] in gsi)
        'lc': [],  # List of subplots with legend common to all subplots
        'suptitle': meta.get('suptitle'),
        'ymax': meta.get('ymax'),
        'ymin': meta.get('ymin'),
        'xmax': meta.get('xmax'),
        'xmin': meta.get('xmin'),
        'sharex': meta.get('sharex') if 'sharex' in meta else 'all',
        'sharey': meta.get('sharey') if 'sharey' in meta else 'all',
        'xscale': meta.get('xscale') if 'xscale' in meta else 'linear',
        'yscale': meta.get('yscale') if 'yscale' in meta else 'linear',
        'xnopos': meta.get('xnopos'),
        'ynopos': meta.get('ynopos'),
        'xmaloc': meta.get('xmaloc'),
        'ymaloc': meta.get('ymaloc'),
        'xmiloc': meta.get('xmiloc'),
        'ymiloc': meta.get('ymiloc'),
        'xmafmt': meta.get('xmafmt'),
        'ymafmt': meta.get('ymafmt'),
        'xmifmt': meta.get('xmifmt'),
        'ymifmt': meta.get('ymifmt'),
        'xformat': meta.get('xformat'),
        'yformat': meta.get('yformat'),
        'xticks': meta.get('xticks'),
        'yticks': meta.get('yticks'),
        'xlabel': meta.get('xlabel'),
        'ylabel': meta.get('ylabel'),
        'axes_used': [],
        'legends': [],
        'legendoutside': meta.get('legendoutside'),
        'legendcollect': meta.get('legendcollect'),
        'legendlocation': meta.get('legendlocation'),
        'legendtitle': meta.get('legendtitle'),
    }
    logger.info('Generated f dict: \n {}'.format(f))
    # Initialize a figure. The size is taken from the stylesheet
    # constrained_layout will make sure to add objects to fill the area
    f['fig'] = plt.figure(constrained_layout=f.get('constrained_layout'), figsize=f.get('figsize'))
    # Set padding between subplots
    # w_pad: Width padding in inches.
    #        This is the pad around Axes and is meant to make sure there is enough room for fonts to look good.
    #        Defaults to 3 pts = 0.04167 inches
    # h_pad: Height padding in inches. Defaults to 3 pts.
    # wspace: Width padding between subplots, expressed as a fraction of the subplot width.
    #         The total padding ends up being w_pad + wspace.
    # hspace: Height padding between subplots, expressed as a fraction of the subplot width.
    #         The total padding ends up being h_pad + hspace.
    # , sharex = True, sharey = True
    # f['fig'].set_constrained_layout_pads(w_pad=1 / 72, h_pad=1 / 72, hspace=0, wspace=0)
    f['fig'].set_constrained_layout_pads(w_pad=0, h_pad=0, hspace=0, wspace=0)

    # Get the outer grid rows and columns
    f['nrows'], f['ncols'] = get_optimal_subplot_num(numplots=numplots)

    # Set width ratios for outer and inner grid rows
    f['owr'] = [1.0] * (f['ncols'] + f['ccols'])
    f['ohr'] = [1.0] * (f['nrows'] + f['crows'])
    f['owr'][0] *= f['owr_pad']  # The left-most subplot needs more space for the ylabel (NEEDS FINE TUNING)
    f['ohr'][-2] *= f['ohr_pad']  # The bottom subplot needs more space for the xlabel (NEEDS FINE TUNING)
    f['owr'][-1] = 0.0001  # The right common area can be small to begin with
    f['ohr'][-1] = 0.0001  # The bottom common area can be small to begin with

    # Create the outer subplot grid  ([g]rid[o]uter)
    f['go'] = f['fig'].add_gridspec(nrows=f['nrows'] + f['crows'], ncols=f['ncols'] + f['ccols'], width_ratios=f['owr'],
                                    height_ratios=f['ohr'], wspace=0.00, hspace=0.00)

    # Generate the outer grid common areas
    f['lc'].append(f['fig'].add_subplot(f['go'][:, -1]))  # The right common area
    f['lc'].append(f['fig'].add_subplot(f['go'][-1, :]))  # The bottom common area
    # Turn off their axes
    for lc in f['lc']:
        lc.axis('off')

    # Keep track of the axes rows/cols
    ax_matrix = np.empty(shape=(f['nrows'], f['ncols']), dtype=plt.Axes)

    # Generate the inner grids
    for ir, ic in np.ndindex(f['nrows'], f['ncols']):
        go = f['go'][ir, ic]
        logger.info('Setting up outer gridspec: row {} col {}'.format(ir, ic))
        f['gi'].append(go.subgridspec(nrows=f['irows'], ncols=f['icols'], width_ratios=f['iwr'], height_ratios=f['ihr'],
                                      wspace=0.02, hspace=0.02))

        gi = f['gi'][-1]

        # Setup axis sharing
        sharex = None
        sharey = None
        if shx := f.get('sharex'):
            if shx == 'col':
                sharex = ax_matrix[0, ic]
            if shx == 'all' or shx is True:
                sharex = ax_matrix[0, 0]
        if shy := f.get('sharey'):
            if shy == 'row':
                sharey = ax_matrix[ir, 0]
            if shy == 'all' or shy is True:
                sharey = ax_matrix[0, 0]

        f['ax'].append(f['fig'].add_subplot(gi[0, 0], sharex=sharex, sharey=sharey))
        f['ix'].append(None)
        f['lr'].append(f['fig'].add_subplot(gi[0, 1]))
        f['lb'].append(f['fig'].add_subplot(gi[1, 0]))
        ax_matrix[ir, ic] = f['ax'][-1]

        # Turn off axis elements on the legend box
        f['lr'][-1].axis('off')
        f['lb'][-1].axis('off')

        if box_aspect := f['box_aspect']:
            f['ax'][-1].set_box_aspect(box_aspect)

        if xscale := f.get('xscale'):
            if xnopos := f.get('xnopos'):
                f['ax'][-1].set_xscale(xscale, nonpositive=xnopos)
            else:
                f['ax'][-1].set_xscale(xscale)
        if yscale := f.get('yscale'):
            if ynopos := f.get('ynopos'):
                f['ax'][-1].set_yscale(yscale, nonpositive=ynopos)
            else:
                f['ax'][-1].set_yscale(yscale)

        if f['ymax'] is not None:
            f['ax'][-1].set_ylim(ymax=f['ymax'])
        if f['ymin'] is not None:
            f['ax'][-1].set_ylim(ymin=f['ymin'])
        if f['xmax'] is not None:
            f['ax'][-1].set_xlim(xmax=f['xmax'])
        if f['xmin'] is not None:
            f['ax'][-1].set_xlim(xmin=f['xmin'])

        if xticks := f.get('xticks'):
            f['ax'][-1].set_xticks(xticks)
        if yticks := f.get('yticks'):
            f['ax'][-1].set_yticks(yticks)
        if xlabel := f.get('xlabel'):
            f['ax'][-1].set_xlabel(xlabel)
        if ylabel := f.get('ylabel'):
            f['ax'][-1].set_ylabel(ylabel)

        # if yformat := f.get('yformat'):
        #     f['ax'][-1].yaxis.set_major_formatter(FormatStrFormatter(yformat))
        # if xformat := f.get('xformat'):
        #     f['ax'][-1].yaxis.set_major_formatter(FormatStrFormatter(xformat))
        # if ymafmt := f.get('ymafmt'):
        #     f['ax'][-1].yaxis.set_major_formatter(ymafmt)
        # if xmafmt := f.get('xmafmt'):
        #     f['ax'][-1].xaxis.set_major_formatter(xmafmt)
        # if ymifmt := f.get('ymifmt'):
        #     f['ax'][-1].yaxis.set_minor_formatter(ymifmt)
        # if xmifmt := f.get('xmifmt'):
        #     f['ax'][-1].xaxis.set_minor_formatter(xmifmt)
        #
        # if ymaloc := f.get('ymaloc'):
        #     f['ax'][-1].yaxis.set_major_locator(ymaloc)
        # if xmaloc := f.get('xmaloc'):
        #     f['ax'][-1].xaxis.set_major_locator(xmaloc)
        # if ymiloc := f.get('ymiloc'):
        #     f['ax'][-1].yaxis.set_minor_locator(ymiloc)
        # if xmiloc := f.get('xmiloc'):
        #     f['ax'][-1].xaxis.set_minor_locator(xmiloc)

        is_last_row = ir + 1 == f['nrows']
        is_first_col = ic == 0
        if not is_last_row and f.get('sharex') in ['col', 'all', True]:
            f['ax'][-1].tick_params(axis='x', which='both', labelbottom=False, zorder=0)
            f['ax'][-1].xaxis.label.set_visible(False)

        if not is_first_col and f.get('sharey') in ['row', 'all', True]:
            f['ax'][-1].tick_params(axis='y', which='both', labelleft=False, zorder=0)
            f['ax'][-1].yaxis.label.set_visible(False)


    # Set title
    if title := f.get('suptitle'):
        f['fig'].suptitle(title)


    legend = {'handle': [], 'label': [], 'title': None, 'header': None}  # One such per legend column

    f['legends'] = {}
    for i, ax in enumerate(f['ax']):
        # We should probably never need more than 10 legend columns per ax
        f['legends'][i] = {}
        for j in range(10):
            f['legends'][i][j] = deepcopy(legend)  # i  for ax, j for legend column

    # m1_legend: {'handle': [], 'label': [], 'ncol': 1, 'unique': True, 'loc': 'upper left', 'insubfig': False,
    #           'title': ["$t_{}$".format('{\ln\ln}')]},

    return f

def get_formatted_columns(columns):
    # columns is a list of legend columns like

    #     [[L, 8, 12, 16 ...],
    #      [n, 2400, 2400, 2400 ...]]
    #
    # The first entry in each legend column is its title
    # This function then returns a list of formatted rows like this
    #
    #       [['L  n   '],
    #        ['8  2400'],
    #        ['12 2400'],
    #        ['16 2400'],
    #
    # Where each column has been formatted to uniform width.

    # First, we check that each list has the same length
    num_cols = len(columns)
    if num_cols == 0:
        return []
    num_rows = len(columns[0])  # Number of labels in each column
    if not all(len(x) == num_rows for x in columns):
        print("columns shape: ", np.shape(columns))
        for i, col in enumerate(columns):
            print("-- col {}: shape {}: {}".format(i, np.shape(col), col))
        raise ValueError("All columns must be of the same length:\n" + str(columns))

    # Now we get the maximum string length of each column
    column_longest = [0] * num_cols
    column_phantom = [''] * num_cols  # The longest string in each column, not including the title
    for icol, col in enumerate(columns):  # Iterate columns
        maxlenstr = max([str(c) for c in col[1:]], key=len)  # Get the longest string in the column
        column_longest[icol] = len(maxlenstr)
        column_phantom[icol] = maxlenstr
    # Now we know how wide each column needs to be, So we can now format each row
    fmt_rows = []
    for irow in range(num_rows):  # Iterate the rows in each column
        fmt_label = ''
        for icol in range(num_cols):  # Iterate columns
            c = columns[icol][irow]  # An entry at irow,icol
            # Let n == len(str(c)) and m == phantom_length[i] be the number of characters in the longest column word
            # Then we need to print c padded with a phantom word of length m-n
            zerolen_string = '\makebox[0pt][l]{{{}}}'.format(str(c))  # Put in text in left-aligned box that does not consume width
            phantom_string = '\hphantom{{{}}}'.format(column_phantom[icol])  # Let the phantom string consume width instead
            colsep_string = ' ' if icol + 1 < num_cols else ''
            fmt_label += '{}{}{}'.format(zerolen_string, phantom_string, colsep_string)
        fmt_rows.append(fmt_label)
    return fmt_rows


def px_to_in(w, h):
    return 27. / 2560 * w, 27. / 1440 * h


def px_to_cm(w, h):
    w, h = px_to_in(w, h)
    return 2.54 * w, 2.54 * h


def labels_are_equal(lgnd_meta):
    # Gather all the labels for the different kinds of legend keys (e.g. "l", "m" type legends)
    labels = {}
    for idx, lgnd in enumerate(lgnd_meta):
        for lkey in lgnd['legend_keys']:
            label = lgnd[lkey]  # The legend description dict for this type of legend
            if not lkey in labels:
                labels[lkey] = []
            labels[lkey].append(label['label'])
    labels_equal = {}
    for idx, (lkey, label_list) in enumerate(labels.items()):
        l0 = label_list[0]
        labels_equal[lkey] = all(l == l0 for l in label_list)
    return labels_equal


def columns_are_equal(fmeta):
    # Gather all the labels for the different kinds of legend keys (e.g. "l", "m" type legends)
    # This will contain True, True, False, True... etc for each column,
    # where True/False tells whether this column is equal on all axes
    numaxes = len(fmeta['legends'].keys())
    numcols = len(fmeta['legends'][0].keys())
    columns_equal = []
    for icol in range(numcols):
        l0 = fmeta['legends'][0][icol]['label']
        columns_equal.append(
            all(fmeta['legends'][0][icol]['label'] == fmeta['legends'][iax][icol]['label'] for iax in range(numaxes)))

    return columns_equal

def add_legend5(fmeta):
    # meta should be a dict with keys 'ax' with the corresponding axes and 'legends' which names the extra legends
    # for (oidx,gso),(nrow,ncol) in zip(enumerate(fmeta['gso']), np.ndindex((fmeta['nrows'], fmeta['ncols']))):
    n = len(fmeta['ax'])
    if not all(len(x) == n for x in [fmeta['ax'], fmeta['lr'], fmeta['lb'], fmeta['legends']]):
        raise AssertionError("Not equal lengths")
    numaxes = len(fmeta['legends'].keys())
    numcols = len(fmeta['legends'][0].keys())
    columns_equal = columns_are_equal(fmeta)
    eqidx = [i for i, x in enumerate(columns_equal) if x is True]  # The indices of columns that are equal
    nqidx = [i for i, x in enumerate(columns_equal) if x is False]  # The indices of columns that are not equal
    # Let's first treat the common columns that can be factored out

    outside = fmeta.get('legendoutside')  # Put legend inside or outside each subplot
    collect = fmeta.get('legendcollect')  # Collect equal columns into a single legend outside
    location = fmeta.get('legendlocation')
    legend_ax = 'lr' if outside else 'ax'
    legend_eq = 'lc' if collect else legend_ax
    legend_nq = legend_ax
    legend_loc = rcParams['legend.loc'] if outside else (location if location is not None else 'best')
    loc_nq = 'upper left' if outside else (location if location is not None else 'best')
    loc_eq = 'center' if collect else loc_nq

    if legend_eq == legend_nq:
        # If legends all go together anyway, we may as well put them together again
        eqidx.extend(nqidx)
        nqidx = []

    iax_tgt = None
    if collect:
        # There may be an unused axis where we can put the common legend instead of lc
        for iax, ax in enumerate(fmeta['ax']):
            if not iax in fmeta['axes_used']:
                legend_eq = 'ax'
                iax_tgt = iax  # Put the common legend on the unused axis instead
                fmeta['axes_used'].append(iax)
                # Remove axis elements
                ax.set(yticklabels=[], xticklabels=[], ylabel=None, xlabel=None, xticks=[], yticks=[])
                ax.tick_params(labelbottom=False, labelleft=False)
                ax.axis('off')
                break

    for iax in range(numaxes):
        # Collect the columns
        columns = []
        handles = []
        for icol in eqidx:
            title = fmeta['legends'][iax][icol]['title']
            label = fmeta['legends'][iax][icol]['label']
            handle = fmeta['legends'][iax][icol]['handle']
            if not title or not label or not handle:
                continue
            columns.append([title] + label)
            handles.extend(handle)
            # Decide where to put it. ax, lr, lg or lc[0]/lc[1]
        if not columns or not handles:
            continue

        # pick a point in the middle, make it invisible with alpha = 0, and move to back with zorder so that
        # the legend entry appears in the beginning for instance with tikzplotlib
        xmin, xmax = fmeta['ax'][iax].get_xlim()
        ymin, ymax = fmeta['ax'][iax].get_ylim()
        titlepatch, = fmeta['ax'][iax].plot([0.5 * (xmin + xmax)], [0.5 * (ymin + ymax)], label=None, alpha=0.0, zorder=0)
        formatted_labels = get_formatted_columns(columns)
        legendtitle = fmeta.get('legendtitle')
        if iax_tgt is not None:
            iax = iax_tgt  # Put legends on common axis
        lg = fmeta[legend_eq][iax].legend(handles=[titlepatch] + handles, labels=formatted_labels,
                                          title=legendtitle,  # title=fmeta['legends'][iax][0]['header'],
                                          loc=loc_eq, prop=dict(stretch="ultra-condensed"))
        # lg._legend_box.align = "right"  # Align the legend title right (i.e. the title row)
        if collect:
            break

    # Now treat the unequal columns. These should go into each subplot
    for iax in range(numaxes):
        # Collect the columns
        columns = []
        handles = []
        for icol in nqidx:
            title = fmeta['legends'][iax][icol]['title']
            label = fmeta['legends'][iax][icol]['label']
            handle = fmeta['legends'][iax][icol]['handle']
            if not title or not label or not handle:
                continue
            columns.append([title] + label)
            handles.extend(handle)
        if not columns or not handles:
            continue
        # pick a point in the middle, make it invisible with alpha = 0, and move to back with zorder so that
        # the legend entry appears in the beginning for instance with tikzplotlib
        xmin, xmax = fmeta['ax'][iax].get_xlim()
        ymin, ymax = fmeta['ax'][iax].get_ylim()
        titlepatch, = fmeta['ax'][iax].plot([0.5 * (xmin + xmax)], [0.5 * (ymin + ymax)], label=None, alpha=0.0, zorder=0)
        formatted_labels = get_formatted_columns(columns)
        lg = fmeta[legend_nq][iax].legend(handles=[titlepatch] + handles, labels=formatted_labels, title=''.join(formatted_labels), loc=loc_nq,
                                          prop=dict(stretch="ultra-condensed"))
        lg._legend_box.align = "right"  # Align the legend title right (i.e. the title row)

    return fmeta


def prettify_plot5(fmeta):
    for idx, (ax, lr, lb) in enumerate(zip(fmeta['ax'], fmeta['lr'], fmeta['lb'])):
        if not idx in fmeta['axes_used']:
            print('Removing axes {}  (used {})'.format(idx, fmeta['axes_used']))
            fmeta['fig'].delaxes(ax)
            fmeta['fig'].delaxes(lr)
            fmeta['fig'].delaxes(lb)
            logger.info('Deleting unused subplot idx {}'.format(idx))
    # Add legends
    add_legend5(fmeta=fmeta)

    # Remove legend boxes that are not used
    for idx, (lr, lb) in enumerate(zip(fmeta['lr'], fmeta['lb'])):
        # Get all legend objects from lr
        objs_lr = [c for c in lr.get_children() if isinstance(c, legend.Legend)]
        objs_lb = [c for c in lb.get_children() if isinstance(c, legend.Legend)]
        if len(objs_lr) == 0 and lr in fmeta['fig'].axes:
            fmeta['fig'].delaxes(lr)
            logger.info('Deleting unused lr idx {}'.format(idx))

        if len(objs_lb) == 0 and lb in fmeta['fig'].axes:
            fmeta['fig'].delaxes(lb)
            logger.info('Deleting unused lb idx {}'.format(idx))

    for idx, lc in enumerate(fmeta['lc']):
        objs_lc = [c for c in lc.get_children() if isinstance(c, legend.Legend)]
        if len(objs_lc) == 0 and lc in fmeta['fig'].axes:
            fmeta['fig'].delaxes(lc)
            logger.info('Deleting unused lc idx {}'.format(idx))

    for ax in fmeta['ax']:
        if yformat := fmeta.get('yformat'):
            ax.yaxis.set_major_formatter(FormatStrFormatter(yformat))
        if xformat := fmeta.get('xformat'):
            ax.yaxis.set_major_formatter(FormatStrFormatter(xformat))
        if ymafmt := fmeta.get('ymafmt'):
            ax.yaxis.set_major_formatter(ymafmt)
        if xmafmt := fmeta.get('xmafmt'):
            ax.xaxis.set_major_formatter(xmafmt)
        if ymifmt := fmeta.get('ymifmt'):
            ax.yaxis.set_minor_formatter(ymifmt)
        if xmifmt := fmeta.get('xmifmt'):
            ax.xaxis.set_minor_formatter(xmifmt)
        if ymaloc := fmeta.get('ymaloc'):
            ax.yaxis.set_major_locator(ymaloc)
        if xmaloc := fmeta.get('xmaloc'):
            ax.xaxis.set_major_locator(xmaloc)
        if ymiloc := fmeta.get('ymiloc'):
            ax.yaxis.set_minor_locator(ymiloc)
        if xmiloc := fmeta.get('xmiloc'):
            ax.xaxis.set_minor_locator(xmiloc)
        # f['ax'][-1].ticklabel_format(axis='x', scilimits=[-3, 20])
        # f['ax'][-1].xaxis.grid(True, which='minor')
        # ax.locator_params(tight=True)
        # ax.autoscale()
        # print("xticks:", ax.get_xticks())


def save_figure(figs):
    if figs is None:
        return
    if isinstance(figs, list):
        for f in figs:
            if not isinstance(f, dict):
                raise TypeError("Expected type(f): dict")
            prettify_plot5(f)
            Path(f['filename']).parent.mkdir(parents=True, exist_ok=True)
            print('Saving figure: {}'.format(f['filename']))
            f['fig'].savefig('{}.pdf'.format(f['filename']), format='pdf')
            f['fig'].savefig('{}.png'.format(f['filename']), format='png')
            f['fig'].savefig('{}.svg'.format(f['filename']), format='svg')
            # f['fig'].savefig('{}.pgf'.format(f['filename']), format='pgf')

            with open('{}.git'.format(f['filename']), 'w') as gitfile:
                repo = Repo(search_parent_directories=True)
                commit_hash = repo.git.rev_parse("HEAD")
                gitfile.write(commit_hash)

            # tikzplotlib.save('{}.tex'.format(f['filename']), figure=f['fig'])
    else:
        raise TypeError("Unexpected type: ", type(figs))

def plt_scatter(x, y, xlabel, ylabel, ax, markerlist, label='', xscale='linear', yscale='linear', markersize=20):
    ax.scatter(x, y, marker=next(markerlist), label=label, alpha=.90, s=markersize)
    if xlabel != '':  ax.set_xlabel(xlabel)
    if ylabel != '':  ax.set_ylabel(ylabel)
    ax.set_xscale(xscale)
    ax.set_yscale(yscale)
    if label != '': ax.legend()

