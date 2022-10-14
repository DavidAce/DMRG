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


logger = logging.getLogger('tools')
import tikzplotlib

matplotlib_version = pkg_resources.get_distribution("matplotlib").version
if version.parse(matplotlib_version) == version.parse("3.5"):
    logger.warning("Matplotlib version {} may not work correctly with gridspec".format(matplotlib_version))


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
    print(sdata)
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
    w2 = wn[1]  # The width of distribution for pairwise interactions. The distribution is either U(J2_mean-w,J2_mean+w) or N(J2_mean,w)
    w3 = wn[2]  # The width of distribution for three-body interactions.

    if r == np.iinfo(np.uint64).max:
        r = L

    tmin1 = 1.0 / w1
    tmax1 = 1.0 / w1

    r2max = np.min([r, int((L - 1) / 2)])  # Number of sites from the center site to the edge site, max(|i-j|)/2
    Jmin2 = np.exp(-(r2max - 1) / x) * w2  # Order of magnitude of the smallest 2-body terms (furthest neighbor, up to L/2)
    Jmax2 = np.exp(-(1 - 1) / x) * w2  # Order of magnitude of the largest 2-body terms (nearest neighbor)
    tmax2 = 0.5 / Jmin2  # (0.5 to improve fits) Time that it takes for the most remote site to interact with the middle
    tmin2 = 1.0 / Jmax2  # Time that it takes for neighbors to interact

    tmin3 = 1.0 / w3
    tmax3 = 1.0 / w3

    tmin = np.max([tmin1, tmin2, tmin3])  # Should this be max? Also, can't have negatives in log-log.
    tmax = np.max([tmax1, tmax2, tmax3])

    tmin = np.min([tmin, tmax])  # Can't have negatives in log-log.
    tmax = np.max([tmin, tmax])  # Make sure the time points are ordered

    idx1 = np.where(tdata >= tmin)[0][0]
    idx2 = np.where(tdata <= tmax)[0][-1]
    if idx2 == len(tdata) - 1:
        idx2 = np.max([idx1, idx2 - 1])
    return idx1, idx2


@njit(parallel=True, cache=True)
def page_entropy(L):
    n = int(2 ** (L / 2))
    S = - float(n - 1) / (2 * n)
    for k in range(1 + n, 1 + n ** 2):  # Include last
        S = S + 1.0 / k
    return S


def get_prop(db, keyfmt, prop):
    if isinstance(keyfmt, list):
        keys = []
        for kf in keyfmt:
            keys.append(db[prop][kf.split(':')[0]] if ':' in kf else db[prop][kf])
        return keys
    else:
        return db[prop][keyfmt.split(':')[0]] if ':' in keyfmt else db[prop][keyfmt]


def get_vals(db, keyfmt):
    return get_prop(db, keyfmt, 'vals')


def get_keys(db, keyfmt):
    return get_prop(db, keyfmt, 'keys')


def get_tex(db, keyfmt):
    return get_prop(db, keyfmt, 'tex')


def get_title(db, keys):
    fmtvals = []
    for s in keys:
        key, fmt = s.split(':') if ':' in s else [s, '']
        if key in db['vals'] and key in db['tex']['keys']:
            if isinstance(db['vals'][key], list):
                fmtvals.append('${}=[{}]$'.format(
                    db['tex']['keys'][key].strip('$'),
                    ','.join('{0:{1}}'.format(val, fmt) for val in db['vals'][key])))
            else:
                fmtvals.append('${}={}$'.format(db['tex']['keys'][key].strip('$'), '{0:{1}}'.format(db['vals'][key], fmt)))
        elif key in db['tex']['eqs']:
            fmtvals.append(db['tex']['eqs'][key])
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
        raise LookupError('No field keys found in table:\n'
                          'Expected: {}\n'
                          'In table: {}'.format(keys, node.dtype.fields.keys()))
    return match


def get_table_data(node, keys=None, dtype='f8', transpose=None):
    if keys == None:
        return node[()].ravel().view((dtype, (1,)))[()], []

    elif node.dtype.fields == None:
        # This is a scalar put inside an ndarray, so that we have
        # a consistent return type
        return np.ndarray(node[()]), []
    else:
        # Find the column names that match keys
        cols = get_safe_table_field_name(node, keys)
        return node.fields(cols)[()].view((dtype, (len(cols),)))[()], cols


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
            print(dbval[key])
            raise
    return legendrow


def get_fig_meta(numplots: int, meta: dict):
    f = {
        'fig': None,
        'filename': meta.get('filename'),
        'constrained_layout': meta.get('constrained_layout'),
        'numplots': numplots,
        'axin': None,  # For insets
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
        'owr_pad': meta.get('owr_pad') if 'owr_pad' in meta else 1.25,
        'ohr_pad': meta.get('ohr_pad') if 'ohr_pad' in meta else 1.25,

        'go': None,  # Outer gridspec
        'gi': [],  # Inner gridspec
        'rc': [],  # List of coordinates "(row,col)" for each ax
        'ax': [],  # List of subplots with plots (i.e. [0,0] in gsi)
        'lr': [],  # List of subplots with legend right (i.e. [0,1] in gsi)
        'lb': [],  # List of subplots with legend below (i.e. [1,0] in gsi)
        'lc': [],  # List of subplots with legend common to all subplots
        'ymax': meta.get('ymax'),
        'ymin': meta.get('ymin'),
        'xmax': meta.get('xmax'),
        'xmin': meta.get('xmin'),
        'sharex': meta.get('sharex') if 'sharex' in meta else 'none',
        'sharey': meta.get('sharey') if 'sharey' in meta else 'none',
        'xscale': meta.get('xscale') if 'xscale' in meta else 'linear',
        'yscale': meta.get('yscale') if 'yscale' in meta else 'linear',
        'xnopos': meta.get('xnopos'),
        'ynopos': meta.get('ynopos'),
        'axes_used': [],
        'legends': [],
        'legendoutside': meta.get('legendoutside'),
        'legendcollect': meta.get('legendcollect'),
        'legendlocation': meta.get('legendlocation'),
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
        f['lr'].append(f['fig'].add_subplot(gi[0, 1]))
        f['lb'].append(f['fig'].add_subplot(gi[1, 0]))
        ax_matrix[ir, ic] = f['ax'][-1]

        # Turn off axis elements on the legend box
        f['lr'][-1].axis('off')
        f['lb'][-1].axis('off')

        is_last_row = ir + 1 == f['nrows']
        is_first_col = ic == 0
        print('sharex {} {}: {}'.format(ir, ic, sharex))
        if not is_last_row and f.get('sharex') in ['col', 'all', True]:
            f['ax'][-1].tick_params(bottom=False, labelbottom=False)
            f['ax'][-1].xaxis.label.set_visible(False)
            # f['ax'][-1].xaxis.set_visible(False)
        if not is_first_col and f.get('sharey') in ['row', 'all', True]:
            f['ax'][-1].tick_params(left=False, labelleft=False)
            f['ax'][-1].yaxis.label.set_visible(False)

        if box_aspect := f['box_aspect']:
            f['ax'][-1].set_box_aspect(box_aspect)

        if f['xnopos']:
            f['ax'][-1].set_xscale(f['xscale'], nonpositive=f['xnopos'])
        else:
            f['ax'][-1].set_xscale(f['xscale'])
        if f['ynopos']:
            f['ax'][-1].set_yscale(f['yscale'], nonpositive=f['ynopos'])
        else:
            f['ax'][-1].set_yscale(f['yscale'])

        if ymax := f['ymax']:
            f['ax'][-1].set_ylim(ymax=ymax)
        if ymin := f['ymin']:
            f['ax'][-1].set_ylim(ymin=ymin)
        if xmax := f['xmax']:
            f['ax'][-1].set_xlim(xmax=xmax)
        if xmin := f['xmin']:
            f['ax'][-1].set_xlim(xmin=xmin)

        # Apply things from meta
        if xticks := meta.get('xticks'):
            f['ax'][-1].set_xticks(xticks)
        if yticks := meta.get('yticks'):
            f['ax'][-1].set_yticks(yticks)
        if 'ylabel' in meta:
            f['ax'][-1].set_ylabel(meta['ylabel'])
        elif 'plotprefix' in meta:
            f['fig'].set_ylabel(meta['plotprefix'])
        else:
            f['fig'].set_ylabel('?')

        if 'xlabel' in meta:
            f['ax'][-1].set_xlabel(meta['xlabel'])

        if 'yformat' in meta:
            f['ax'][-1].yaxis.set_major_formatter(FormatStrFormatter(meta['yformat']))
        if 'xformat' in meta:
            f['ax'][-1].yaxis.set_major_formatter(FormatStrFormatter(meta['xformat']))


    # Set title
    if title := meta.get('suptitle'):
        f['fig'].suptitle(title)


    legend = {'handle': [], 'label': [], 'title': None}  # One such per legend column

    f['legends'] = {}
    for i, ax in enumerate(f['ax']):
        # We should probably never need more than 10 legend columns per ax
        f['legends'][i] = {}
        for j in range(10):
            f['legends'][i][j] = deepcopy(legend)  # i  for ax, j for legend column

    # m1_legend: {'handle': [], 'label': [], 'ncol': 1, 'unique': True, 'loc': 'upper left', 'insubfig': False,
    #           'title': ["$t_{}$".format('{\ln\ln}')]},

    return f


def add_legend(ax, handles_labels, title=None, loc=None, ncol=1):
    unique_handles = []
    unique_labels = []
    for handle, label in zip(handles_labels['line'], handles_labels['text']):
        if not label in unique_labels:
            unique_handles.append(handle)
            unique_labels.append(label)
    second_legend = plt.legend(handles=unique_handles, labels=unique_labels,
                               title=title, loc=loc,
                               framealpha=0.7, fontsize='x-small', labelspacing=0.25, ncol=ncol)
    ax.add_artist(second_legend)


def add_legend2(ax, legend_meta):
    loc = legend_meta['loc'] if 'loc' in legend_meta else rcParams['legend.loc']
    ncol = legend_meta['ncol'] if 'ncol' in legend_meta else 1
    title = legend_meta['title'] if 'title' in legend_meta else None
    fontsize = legend_meta['fontsize'] if 'fontsize' in legend_meta else rcParams['legend.fontsize']
    framealpha = legend_meta['framealpha'] if 'framealpha' in legend_meta else rcParams['legend.framealpha']
    labelspacing = legend_meta['labelspacing'] if 'labelspacing' in legend_meta else rcParams['legend.labelspacing']
    unique = legend_meta['unique'] if 'unique' in legend_meta else True
    insubfig = legend_meta['insubfig'] if 'insubfig' in legend_meta else False
    handles = []
    labels = []
    if insubfig:
        legend_objects = []
        for ax, handle, label in zip(legend_meta['ax'], legend_meta['handle'], legend_meta['label']):
            handle.set_label(label)
        ax.legend(title=title, loc=loc, framealpha=framealpha,
                  fontsize=fontsize, labelspacing=labelspacing, ncol=ncol)
    else:
        for handle, label in zip(legend_meta['handle'], legend_meta['label']):
            if unique:
                if not label in labels:
                    handles.append(handle)
                    labels.append(label)
            else:
                handles.append(handle)
                labels.append(label)

        legend_object = plt.legend(handles=handles, labels=labels,
                                   title=title, loc=loc, framealpha=framealpha,
                                   fontsize=fontsize, labelspacing=labelspacing, ncol=ncol)
        ax.add_artist(legend_object)


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
    print(columns)
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


def prettify_plot(fig, axes, cols, rows, axes_used, xmax=None, xmin=None, ymax=None, ymin=None, nlabel=None,
                  ntitle=None, Llabel=None, Ltitle=None, ncol=1, loc=None, legend_in_subfig=True):
    for idx, ax in enumerate(np.ravel(axes)):  # Last one is for the legend
        if not idx in axes_used:
            fig.delaxes(ax)
        else:
            ax.set_ylim(ymin=ymin, ymax=ymax)
            ax.set_xlim(xmin=xmin, xmax=xmax)

    handles_labels = [ax.get_legend_handles_labels() for ax in fig.axes]
    unique_handles = []
    unique_labels = []
    for handles, labels in handles_labels:
        for handle, label in zip(handles, labels):
            if not label in unique_labels:
                unique_handles.append(handle)
                unique_labels.append(label)

    if legend_in_subfig:
        if not loc:
            loc = 'best'
        for idx, ax in enumerate(np.ravel(axes)):
            print('loc:', loc)
            ax.legend(loc=loc, ncol=ncol)
        return
    else:
        # Create a legend for the first line.
        if not loc:
            if nlabel or Llabel:
                loc = 'upper center'
            else:
                loc = 'center right'

        ax = fig.add_subplot(rows, cols + 1, rows * cols + 1)

        first_legend = plt.legend(handles=unique_handles, labels=unique_labels, loc=loc, ncol=ncol)

        # Add the legend manually to the current Axes.
        plt.gca().add_artist(first_legend)

        if nlabel:
            add_legend(ax, nlabel, ntitle, loc='best', ncol=2)

        if Llabel:
            add_legend(ax, Llabel, Ltitle, loc='best', ncol=2)
        ax.axis('off')


def prettify_plot2(fig, axes, cols, rows, axes_used, xmax=None, xmin=None, ymax=None, ymin=None, extra_legend=None):
    for idx, ax in enumerate(np.ravel(axes)):  # Last one is for the legend
        if not idx in axes_used:
            fig.delaxes(ax)
        else:
            ax.set_ylim(ymin=ymin, ymax=ymax)
            ax.set_xlim(xmin=xmin, xmax=xmax)

    ax = fig.add_subplot(rows, cols, rows * cols)
    ax.axis('off')

    # Create a legend for the first line.
    if extra_legend:
        if isinstance(extra_legend, list):
            for l in extra_legend:
                add_legend2(ax, l)
        elif isinstance(extra_legend, dict):
            add_legend2(ax, extra_legend)


def add_legend3(fig, lgnd_meta, lgnd_list=None):
    # meta should be a dict with keys 'ax' with the corresponding axes and 'legends' which names the extra legends

    for lkey in lgnd_meta['legend_keys']:
        ax = lgnd_meta['ax']
        l = lgnd_meta[lkey]
        handles = l['handle']
        labels = l['label']
        loc = l['loc'] if 'loc' in l else rcParams['legend.loc']
        ncol = l['ncol'] if 'ncol' in l else 1
        title = l['title'] if 'title' in l else None
        fontsize = l['fontsize'] if 'fontsize' in l else rcParams['legend.fontsize']
        framealpha = l['framealpha'] if 'framealpha' in l else rcParams['legend.framealpha']
        labelspacing = l['labelspacing'] if 'labelspacing' in l else rcParams['legend.labelspacing']
        unique = l['unique'] if 'unique' in l else True
        insubfig = l['insubfig'] if 'insubfig' in l else False

        if not insubfig:
            gs = ax.get_gridspec()
            rows = gs.nrows
            cols = gs.ncols
            # ax = fig.add_subplot(rows, cols, rows * cols)
            # ax.axis('off')
            raise AssertionError("insubfig==false has not been implemented in add_legend3")

        if len(handles) == 0:
            continue
        print(lgnd_meta)
        # Detects wether this is an 'l' or an 'm' plot
        # 'l' plots are for the main data characteristics, such as 'L', 'num', bond dimensions, etc.
        # 'm' is to show the saturation window
        if 'l' in lkey:
            anchor = (1.01, 1.01)  # anchor the upper left corner of the legend to this coordinate on the subplot
        else:
            anchor = (1.01, 0.5)  # anchor the upper left corner of the legend to this coordinate on the subplot

        # Shrink current axis by 20%
        # box = ax.get_position()
        # ax.set_position([box.x0, box.y0, box.width * 0.5, box.height])
        dummy = mpatches.Patch(alpha=0.0, label=None)

        formatted_labels = get_formatted_columns([title] + labels)

        lg = ax.legend(handles=[dummy] + handles, labels=formatted_labels,
                       title=None, framealpha=framealpha,
                       fontsize=fontsize, labelspacing=labelspacing, ncol=ncol,
                       bbox_to_anchor=anchor,
                       bbox_transform=ax.transAxes,
                       borderaxespad=0.
                       )
        lg._legend_box.align = "right"  # Align the legend title right (i.e. the dummy row)
        ax.add_artist(lg)
        if isinstance(lgnd_list, list):
            lgnd_list.append(lg)

    # else:
    #     for handle, label in zip(legend_meta['handle'], legend_meta['label']):
    #         if unique:
    #             if not label in labels:
    #                 handles.append(handle)
    #                 labels.append(label)
    #         else:
    #             handles.append(handle)
    #             labels.append(label)
    #
    #     legend_object = plt.legend(handles=handles, labels=labels,
    #                                title=title, loc=loc, framealpha=framealpha,
    #                                fontsize=fontsize, labelspacing=labelspacing, ncol=ncol)
    #     ax.add_artist(legend_object)


def prettify_plot3(fig, axes, cols, rows, axes_used, xmax=None, xmin=None, ymax=None, ymin=None, lgnd_meta=None, lgnd_list=None):
    for idx, ax in enumerate(np.ravel(axes)):  # Last one is for the legend
        if not idx in axes_used:
            fig.delaxes(ax)
        else:
            ax.set_ylim(ymin=ymin, ymax=ymax)
            ax.set_xlim(xmin=xmin, xmax=xmax)

    # Create a legend for the first line.
    if lgnd_meta:
        if isinstance(lgnd_meta, list):
            for l in lgnd_meta:
                add_legend3(fig=fig, lgnd_meta=l, lgnd_list=lgnd_list)
        elif isinstance(lgnd_meta, dict):
            add_legend3(fig=fig, lgnd_meta=lgnd_meta, lgnd_list=lgnd_list)

        # cols = 1
        # if(not isinstance(axes,plt.Axes)): # axes can be Axes or array
        #     cols = np.shape(axes)[-1]

        # Where is the right edge of the last subfigure, excluding legend, as a fraction of the whole figure?
        # It depends on cols. Example: Let cols == 3, and rect=[0,0,0.75,1] then we should fit 3 subfigures,
        # where 2 of them have a legend included (due to tigth_layout(rect=[...])) but the third legend would
        # not fit inside the figure. So then, interpreting 100% as the figure including the legend, the right
        # edge must be at (3-1)+ 0.75 /3
        # frac = 0.75
        # redg = 0.95*((cols - 1) + frac)/cols
        # fig.tight_layout(rect=[0, 0, frac, 1]) # Reshapes each subplot to be 75% of its original width, to fit a legend
        # fig.subplots_adjust(wspace=frac, right=redg) # Reshape the figure so that the right-most legend will fit. 0.8 is a fraction of the total figure size

        if lgnd_list:
            # Figure out w_pad for tight layout
            fig.canvas.draw()
            # figsize = fig.get_size_inches() * fig.dpi  # figure (width,height) in pixels
            font_size_in_pixels = rcParams['font.size']
            # Get the widest legend
            w_max_in_pixels = 0.
            h_max_in_pixels = 0.
            # print('lgnd_list', lgnd_list)

            for lg in lgnd_list:
                # print('getting lg frame', lg)
                f = lg.get_frame()
                w_max_in_pixels = max(w_max_in_pixels, f.get_width())
                h_max_in_pixels = max(h_max_in_pixels, f.get_height())

            # We have the widest legend box in units of pixels. Convert to units of font size
            w_max_in_centimeter = w_max_in_pixels / 2560 * 27 * 2.54
            w_max_in_font_units = w_max_in_pixels / font_size_in_pixels
            h_max_in_font_units = h_max_in_pixels / font_size_in_pixels
            print("font_size_in_pixels: ", font_size_in_pixels)
            print("w_max_in_centimeter: ", w_max_in_centimeter)
            print("w_max_in_pixels    : ", w_max_in_pixels)
            print("w_max_in_font_units: ", w_max_in_font_units)

            ax_w_max_in_pixels = 0.
            ax_h_max_in_pixels = 0.
            for lmeta in lgnd_meta:
                ax = lmeta['ax']
                ax_tbox = ax.get_tightbbox(fig.canvas.get_renderer())
                # ax_tbox = ax.get_tightbbox(fig.canvas.get_renderer()).transformed(fig.dpi_scale_trans.inverted())
                # ax_bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
                ax_bbox = ax.get_window_extent()
                w_bbox, h_bbox = ax_bbox.width, ax_bbox.height
                w_tbox, h_tbox = ax_tbox.width, ax_tbox.height
                # ax_w_max_in_pixels = max(ax_w_max_in_pixels, width)
                # ax_h_max_in_pixels = max(ax_h_max_in_pixels, height)
                print('bbox: {} px | {} in | {} cm'.format((w_bbox, h_bbox), px_to_in(w_bbox, h_bbox), px_to_cm(w_bbox, h_bbox)))
                print('tbox: {} px | {} in | {} cm'.format((w_tbox, h_tbox), px_to_in(w_tbox, h_tbox), px_to_cm(w_tbox, h_tbox)))

            # We have the widest axes box in units of pixels. Convert to units of font size
            # w_max_in_centimeter = w_max_in_pixels / 2560*27*2.54
            # w_max_in_font_units = w_max_in_pixels / font_size_in_pixels
            # h_max_in_font_units = h_max_in_pixels / font_size_in_pixels
            # print("font_size_in_pixels: ", font_size_in_pixels)
            # print("w_max_in_centimeter: ", w_max_in_centimeter)
            # print("w_max_in_pixels    : ", w_max_in_pixels)
            # print("w_max_in_font_units: ", w_max_in_font_units)

            figcols = 1
            if (not isinstance(axes, plt.Axes)):  # axes can be Axes or array
                figcols = np.shape(axes)[-1]

            # Where is the right edge of the last subfigure, excluding legend, as a fraction of the whole figure?
            # It depends on figcols. Example: Let cols == 3, and rect=[0,0,0.75,1] then we should fit 3 subfigures,
            # where 2 of them have a legend included (due to tigth_layout(rect=[...])) but the third legend would
            # not fit inside the figure. So then, interpreting 100% as the figure including the legend, the right
            # edge must be at (3-1)+ 0.75 /3
            frac = 1
            redg = 0.95 * ((figcols - 1) + frac) / figcols
            # fig.tight_layout(rect=[0, 0, frac, 1]) # Reshapes each subplot to be 75% of its original width, to fit a legend
            # fig.tight_layout(w_pad = w_max_in_font_units,
            # rect=[0, 0, frac, 1]
            # ) # h_pad/w_pad: padding (height/width) between edges of adjacent subplots, as a fraction of the font size.
            # Reshape the figure so that the right-most legend will fit. 0.8 is a fraction of the total figure
            # fig.subplots_adjust(
            # wspace=frac,
            # right=redg
            # )


    else:
        fig.tight_layout()


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


def add_legend4(fmeta, lgnd_meta):
    # meta should be a dict with keys 'ax' with the corresponding axes and 'legends' which names the extra legends
    # for (oidx,gso),(nrow,ncol) in zip(enumerate(fmeta['gso']), np.ndindex((fmeta['nrows'], fmeta['ncols']))):
    n = len(fmeta['ax'])
    if not all(len(x) == n for x in [fmeta['ax'], fmeta['lr'], fmeta['lb'], lgnd_meta]):
        raise AssertionError("Not equal lengths")
    labels_equal = labels_are_equal(lgnd_meta)
    for idx, (ax, lr, lb, lgnd) in enumerate(zip(fmeta['ax'], fmeta['lr'], fmeta['lb'], lgnd_meta)):
        for lkey in lgnd['legend_keys']:
            l = lgnd[lkey]  # The legend description dict
            handles = l['handle']
            labels = l['label']
            loc = l['loc'] if 'loc' in l else rcParams['legend.loc']
            ncol = l['ncol'] if 'ncol' in l else 1
            title = l['title'] if 'title' in l else None
            fontsize = l['fontsize'] if 'fontsize' in l else rcParams['legend.fontsize']
            framealpha = l['framealpha'] if 'framealpha' in l else rcParams['legend.framealpha']
            labelspacing = l['labelspacing'] if 'labelspacing' in l else rcParams['legend.labelspacing']
            unique = l['unique'] if 'unique' in l else True
            insubfig = l['insubfig'] if 'insubfig' in l else False

            if len(handles) == 0:
                continue
            if labels_equal[lkey] and idx > 0:
                continue

            titlepatch = mpatches.Patch(alpha=0.0, label=None)
            formatted_labels = get_formatted_columns([title] + labels)

            # Decide where to put it. ax, lr, lg or lc[0]/lc[1]
            if insubfig:
                lg = ax.legend(handles=[titlepatch] + handles, labels=formatted_labels, title=None, ncol=ncol)
            elif labels_equal[lkey] and 'l' in lkey:
                lg = fmeta['lc'][0].legend(handles=[titlepatch] + handles, labels=formatted_labels, title=None, ncol=ncol)
            elif labels_equal[lkey] and 'm' in lkey:
                lg = fmeta['lc'][1].legend(handles=[titlepatch] + handles, labels=formatted_labels, title=None, ncol=ncol)
            else:
                if 'l' in lkey:
                    lg = lr.legend(handles=[titlepatch] + handles, labels=formatted_labels, title=None, ncol=ncol)
                else:
                    lg = lb.legend(handles=[titlepatch] + handles, labels=formatted_labels, title=None, ncol=ncol)

            lg._legend_box.align = "right"  # Align the legend title right (i.e. the dummy row)
            # ax.add_artist(lg)


def prettify_plot4(fmeta, lgnd_meta=None):
    for idx, (ax, lr, lb) in enumerate(zip(fmeta['ax'], fmeta['lr'], fmeta['lb'])):
        if not idx in fmeta['axes_used']:
            fmeta['fig'].delaxes(ax)
            fmeta['fig'].delaxes(lr)
            fmeta['fig'].delaxes(lb)

    # Add legends
    if lgnd_meta:
        add_legend4(fmeta=fmeta, lgnd_meta=lgnd_meta)

        # Remove legend boxes that are not used
        for lr, lb in zip(fmeta['lr'], fmeta['lb']):
            # Get all legend objects from lr
            objs_lr = [c for c in lr.get_children() if isinstance(c, legend.Legend)]
            objs_lb = [c for c in lb.get_children() if isinstance(c, legend.Legend)]
            if (len(objs_lr) == 0):
                fmeta['fig'].delaxes(lr)
            if (len(objs_lb) == 0):
                fmeta['fig'].delaxes(lb)

        for lc in fmeta['lc']:
            objs_lc = [c for c in lc.get_children() if isinstance(c, legend.Legend)]
            if (len(objs_lc) == 0):
                fmeta['fig'].delaxes(lc)


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
        if iax_tgt is not None:
            iax = iax_tgt  # Put legends on common axis
        lg = fmeta[legend_eq][iax].legend(handles=[titlepatch] + handles, labels=formatted_labels, title=None, loc=loc_eq)
        lg._legend_box.align = "right"  # Align the legend title right (i.e. the title row)
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
        lg = fmeta[legend_nq][iax].legend(handles=[titlepatch] + handles, labels=formatted_labels, title=''.join(formatted_labels), loc=loc_nq)
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
    return

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


def save_figure(f):
    prettify_plot5(f)
    f['fig'].savefig('{}.pdf'.format(f['filename']), format='pdf')
    f['fig'].savefig('{}.png'.format(f['filename']), format='png')
    f['fig'].savefig('{}.svg'.format(f['filename']), format='svg')
    f['fig'].savefig('{}.pgf'.format(f['filename']), format='pgf')
    tikzplotlib.save('{}.tex'.format(f['filename']), figure=f['fig'])


def get_markerlist():
    return itertools.cycle(('^', 'v', '<', '>', 'o'))


def plt_scatter(x, y, xlabel, ylabel, ax, markerlist, label='', xscale='linear', yscale='linear', markersize=20):
    ax.scatter(x, y, marker=next(markerlist), label=label, alpha=.90, s=markersize)
    if xlabel != '':  ax.set_xlabel(xlabel)
    if ylabel != '':  ax.set_ylabel(ylabel)
    ax.set_xscale(xscale)
    ax.set_yscale(yscale)
    if label != '': ax.legend()


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
