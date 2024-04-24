import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import rcParams, legend
from matplotlib.ticker import NullFormatter, FormatStrFormatter, ScalarFormatter, LogLocator, LinearLocator
from copy import deepcopy
import itertools


def write_attributes(*args, **kwargs):
    for a in args:
        print(a)
    for k, v in kwargs.items():
        print("%s = %s" % (k, v))


def get_data(node, keys=None, dtype='f8'):
    if node.dtype.fields == None or keys == None:
        if np.isscalar(node[()]):
            # This is a scalar. We put it in a ndarray, so that we have a consistent return type
            return np.ndarray(node[()])
        else:
            # This is an array dataset. Use ravel to make sure we get a 1d-array
            return np.ravel(node[()])
    else:
        # Find the column names that match key
        cols = []
        if isinstance(keys, str):
            cols = [x for x in node.dtype.fields.keys() if keys in x]
            if not cols:
                raise ValueError("Could not find fields matching [{}] in node [{}] with fields: {}".format(keys, node.name, node.dtype.fields.keys()))

            # Get the dtype of the first column. All the others should match or this won't work
            # dtype = node.dtype.fields[cols[0]][0]
            return node.fields(cols)[()].view(dtype)[()]
        elif isinstance(keys, list):
            cols = [x for x in node.dtype.fields.keys() if any([k in x for k in keys])]
            return node.fields(cols)[()].view((dtype, (len(cols),)))[()]


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


def get_matching_path(db, keys):
    matches = []
    for path, node in db['dsets'].items():
        ismatch = True
        for key in keys:
            if isinstance(key, list):
                if not any('/{}'.format(k) in path for k in key):
                    ismatch = False
                    break
            elif not '/{}'.format(key) in path:
                ismatch = False
                break
        if ismatch:
            matches.append(node['path']['data'])
    return matches


def get_matching_data(h5, grps, dset, keys=None, dtype='f8'):
    data = []
    for grp in grps:
        data.extend(get_data(h5[grp][dset], keys=keys, dtype=dtype))
    return data


def get_matching_vals(h5, grps, db, keys):
    data = []
    for grp in grps:
        dbval = db['dsets'][grp]
        vals = get_vals(dbval, keys)
        if isinstance(vals, np.ndarray):
            data.extend(vals)
        elif np.isscalar(vals):
            data.append(vals)
    return data


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


def get_legend_rows(db, h5, grps, legend_col_keys):
    legendrows = []
    for grp in grps:
        legendrows.append(get_legend_row(db, h5[grp], legend_col_keys=legend_col_keys))
    return legendrows


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


def get_fig_meta(numplots: int, meta: dict):
    f = {
        'fig': None,
        'numplots': numplots,
        'nrows': None,
        'ncols': None,
        'crows': 1,  # Common row put at the bottom of all subplots
        'ccols': 1,  # Common col put at the right of all subplots
        'irows': 2,  # We make a 2x2 grid where 0,0 is the plot, and 0,1 and 1,0 are legends
        'icols': 2,  # We make a 2x2 grid where 0,0 is the plot, and 0,1 and 1,0 are legends

        'owr': None,
        'ohr': None,
        'iwr': [1000, 1],  # Width ratio between plot and right legend
        'ihr': [1000, 1],  # Height ratio between plot and bottom legend
        'owr_pad': meta.get('owr_pad') if 'owr_pad' in meta else 1.20,  # The left-most subplot needs more space for the ylabel (NEEDS FINE TUNING)
        'owl_pad': meta.get('owl_pad') if 'owl_pad' in meta else 1.20,  # The bottom subplot needs more space for the xlabel (NEEDS FINE TUNING)

        'go': None,  # Outer gridspec
        'gi': [],  # Inner gridspec
        'ax': [],  # List of subplots with plots (i.e. [0,0] in gsi)
        'lr': [],  # List of subplots with legend right (i.e. [0,1] in gsi)
        'lb': [],  # List of subplots with legend below (i.e. [1,0] in gsi)
        'lc': [],  # List of subplots with legend common to all subplots
        # 'ymax': meta.get('ymax'),
        # 'ymin': meta.get('ymin'),
        # 'xmax': meta.get('xmax'),
        # 'xmin': meta.get('xmin'),
        'sharex': meta.get('sharex'),  # if 'sharex' in meta else 'none',
        'sharey': meta.get('sharey'),  # if 'sharey' in meta else 'none',
        'xscale': meta.get('xscale') if 'xscale' in meta else 'linear',
        'yscale': meta.get('yscale') if 'yscale' in meta else 'linear',
        'xlabel': meta.get('xlabel'),
        'ylabel': meta.get('ylabel'),
        'xticks': meta.get('xticks'),
        'yticks': meta.get('yticks'),
        'xformat': meta.get('xformat'),
        'yformat': meta.get('yformat'),
        'xmajor_formatter': meta.get('xmajor_formatter'),
        'ymajor_formatter': meta.get('ymajor_formatter'),
        'xminor_formatter': meta.get('xminor_formatter'),
        'yminor_formatter': meta.get('yminor_formatter'),
        'axes_used': [],
        'legends': [],
        'legendoutside': meta.get('legendoutside'),
        'legendcollect': meta.get('legendcollect'),
    }
    # Initialize a figure. The size is taken from the stylesheet
    # constrained_layout will make sure to add objects to fill the area
    f['fig'] = plt.figure(constrained_layout=True)
    # Set padding between subplots
    # w_pad: Width padding in inches.
    #        This is the pad around Axes and is meant to make sure there is enough room for fonts to look good.
    #        Defaults to 3 pts = 0.04167 inches
    # h_pad: Height padding in inches. Defaults to 3 pts. Helps distance title from suptitle
    # wspace: Width padding between subplots, expressed as a fraction of the subplot width.
    #         The total padding ends up being w_pad + wspace.
    # hspace: Height padding between subplots, expressed as a fraction of the subplot width.
    #         The total padding ends up being h_pad + hspace.
    # , sharex = True, sharey = True
    # f['fig'].set_constrained_layout_pads(w_pad=1 / 72, h_pad=1 / 72, hspace=0, wspace=0)
    f['fig'].set_constrained_layout_pads(w_pad=0, h_pad=0.015, hspace=0, wspace=0)

    # Get the outer grid rows and columns
    f['nrows'], f['ncols'] = get_optimal_subplot_num(numplots=numplots)

    # Set width ratios for outer and inner grid rows
    f['owr'] = [1.0] * (f['ncols'] + f['ccols'])
    f['ohr'] = [1.0] * (f['nrows'] + f['crows'])
    if f['sharey'] == 'all' or f['sharey'] == 'row':
        f['owr'][0] *= f['owr_pad']  # The left-most subplot needs more space for the ylabel (NEEDS FINE TUNING)
    if f['sharex'] == 'all' or f['sharex'] == 'col':
        f['ohr'][-2] *= f['owl_pad']  # The bottom subplot needs more space for the xlabel (NEEDS FINE TUNING)
    f['owr'][-1] = 0.001  # The right common area can be small to begin with
    f['ohr'][-1] = 0.001  # The bottom common area can be small to begin with

    # Create the outer subplot grid  ([g]rid[o]uter)
    f['go'] = f['fig'].add_gridspec(nrows=f['nrows'] + f['crows'], ncols=f['ncols'] + f['ccols'], width_ratios=f['owr'],
                                    height_ratios=f['ohr'], wspace=0.00, hspace=5.00)

    # Generate the outer grid common areas
    f['lc'].append(f['fig'].add_subplot(f['go'][:, -1]))  # The right common area
    f['lc'].append(f['fig'].add_subplot(f['go'][-1, :]))  # The bottom common area
    f['lc'][0].set(yticklabels=[], xticklabels=[], ylabel=None, xlabel=None, xticks=[], yticks=[])
    f['lc'][0].tick_params(labelbottom=False, labelleft=False)
    f['lc'][1].set(yticklabels=[], xticklabels=[], ylabel=None, xlabel=None, xticks=[], yticks=[])
    f['lc'][1].tick_params(labelbottom=False, labelleft=False)

    # Generate the inner grids
    for ir, ic in np.ndindex(f['nrows'], f['ncols']):
        go = f['go'][ir, ic]
        f['gi'].append(go.subgridspec(nrows=f['irows'], ncols=f['icols'], width_ratios=f['iwr'], height_ratios=f['ihr'],
                                      wspace=0.01, hspace=0.01))

        gi = f['gi'][-1]
        f['ax'].append(f['fig'].add_subplot(gi[0, 0]))
        f['lr'].append(f['fig'].add_subplot(gi[0, 1]))
        f['lb'].append(f['fig'].add_subplot(gi[1, 0]))

        # Turn off axis elements on the legend box
        f['lr'][-1].set(yticklabels=[], xticklabels=[], ylabel=None, xlabel=None, xticks=[], yticks=[])
        f['lb'][-1].set(yticklabels=[], xticklabels=[], ylabel=None, xlabel=None, xticks=[], yticks=[])
        f['lr'][-1].tick_params(labelbottom=False, labelleft=False)
        f['lb'][-1].tick_params(labelbottom=False, labelleft=False)

        f['ax'][-1].set_xscale(f['xscale'])
        f['ax'][-1].set_yscale(f['yscale'])

        # Apply things from meta
        if xticks := f.get('xticks'):
            f['ax'][-1].set_xticks(xticks)
        if yticks := f.get('yticks'):
            f['ax'][-1].set_yticks(yticks)
        f['ax'][-1].set_xlabel(f['xlabel'])
        f['ax'][-1].set_ylabel(f['ylabel'])
        if yformat := f.get('yformat'):
            f['ax'][-1].yaxis.set_major_formatter(FormatStrFormatter(yformat))
        if xformat := f.get('xformat'):
            f['ax'][-1].xaxis.set_major_formatter(FormatStrFormatter(xformat))

        if formatter := f.get('ymajor_formatter'):
            f['ax'][-1].yaxis.set_major_formatter(formatter)
        if formatter := f.get('yminor_formatter'):
            f['ax'][-1].yaxis.set_minor_formatter(formatter)

        if formatter := f.get('xmajor_formatter'):
            f['ax'][-1].xaxis.set_major_formatter(formatter)
        if formatter := f.get('xminor_formatter'):
            f['ax'][-1].xaxis.set_minor_formatter(formatter)

        # locmaj = LogLocator(base=10.0, subs=(1.0,), numticks=None)
        # f['ax'][-1].yaxis.set_major_locator(locmaj)
        #
        # locmin = LogLocator(base=10.0, subs=np.arange(2, 10) * .1,numticks=None)
        # f['ax'][-1].yaxis.set_minor_locator(locmin)
        # f['ax'][-1].yaxis.set_minor_formatter(NullFormatter())

        # f['ax'][-1].xaxis.set_major_formatter(ScalarFormatter())
        # locmaj = LinearLocator(base=10.0, subs=(4.0,), numticks=None)
        # f['ax'][-1].xaxis.set_major_locator(locmaj)
        #
        # locmin = LinearLocator(base=10.0, subs=(1.0,),numticks=None)
        # f['ax'][-1].xaxis.set_minor_locator(locmin)
        # f['ax'][-1].xaxis.set_minor_formatter(NullFormatter())

        # Set aspect ratios
        f['ax'][-1].set(aspect='auto', anchor='C')

        # Hide y labels in the bulk
        if sharey := f.get('sharey'):
            if ic > 0 and (sharey == 'row' or f['sharey'] == 'all'):
                f['ax'][-1].set(yticklabels=[], ylabel=None)
        # Hide x labels in the bulk
        if sharex := f.get('sharex'):
            if ir + 1 < f['nrows'] and (sharex == 'col' or f['sharex'] == 'all'):
                f['ax'][-1].set(xticklabels=[], xlabel=None)

    # Set title
    if title := meta.get('suptitle'):
        f['fig'].suptitle(title)

    # Let all subplots share axes
    if 'sharey' in f and f['sharey'] != 'none':
        f['ax'][0].get_shared_y_axes().join(*f['ax'])
    if 'sharex' in f and f['sharex'] != 'none':
        f['ax'][0].get_shared_x_axes().join(*f['ax'])

    # Hide axes for legend areas
    for lr, lb in zip(f['lr'], f['lb']):
        lr.axis('off')
        lb.axis('off')
    for lc in f['lc']:
        lc.axis('off')

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


def prettify_plot(fig, axes, cols, rows, axes_used, xmax=None, xmin=None, ymax=None, ymin=None, nlabel=None):
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
    ax = fig.add_subplot(rows, cols, rows * cols)
    # Create a legend for the first line.
    if nlabel:
        loc = 'upper center'
    else:
        loc = 'center'

    first_legend = plt.legend(handles=unique_handles, labels=unique_labels,
                              loc=loc, labelspacing=0.25, fontsize='small')

    # Add the legend manually to the current Axes.
    plt.gca().add_artist(first_legend)

    if nlabel:
        second_legend = plt.legend(handles=nlabel['line'], labels=nlabel['text'],
                                   title='Realizations', loc='lower center',
                                   framealpha=0.7, fontsize='x-small', labelspacing=0.25, ncol=1)
        plt.gca().add_artist(second_legend)
    ax.axis('off')


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


def columns_are_equal(fmeta):
    # Gather all the labels for the different kinds of legend keys (e.g. "l", "m" type legends)
    # This will contain True, True, False, True... etc for each column,
    # where True/False tells whether this column is equal on all axes
    numaxes = len(fmeta['legends'].keys())
    numcols = len(fmeta['legends'][0].keys())
    columns_equal = []
    for icol in range(numcols):
        l0 = fmeta['legends'][0][icol]['label']
        equal = []
        for iax in range(numaxes):
            li = fmeta['legends'][iax][icol]['label']
            if len(l0) == len(li):  # Skip empty li unless l0 is also empty
                equal.append(l0 == li)
        columns_equal.append(all(equal))
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

    legend_ax = 'lr' if outside else 'ax'
    legend_eq = 'lc' if collect else legend_ax
    legend_nq = legend_ax
    legend_loc = rcParams['legend.loc'] if outside else 'best'
    loc_nq = 'upper left' if outside else 'best'
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
        titlepatch = mpatches.Patch(alpha=0.0, label=None)
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
        # print("columns:",columns)
        # print("handles:",handles)
        titlepatch = mpatches.Patch(alpha=0.0, label=None)
        formatted_labels = get_formatted_columns(columns)
        # print("formatted labels:",formatted_labels)
        lg = fmeta[legend_nq][iax].legend(handles=[titlepatch] + handles, labels=formatted_labels, title=None, loc=loc_nq)
        lg._legend_box.align = "right"  # Align the legend title right (i.e. the title row)

    return fmeta


def prettify_plot5(fmeta):
    # Add legends
    fmeta = add_legend5(fmeta=fmeta)  # May modify fmeta

    for idx, (ax, lr, lb) in enumerate(zip(fmeta['ax'], fmeta['lr'], fmeta['lb'])):
        if not idx in fmeta['axes_used']:
            fmeta['fig'].delaxes(ax)
            fmeta['fig'].delaxes(lr)
            fmeta['fig'].delaxes(lb)

    # Remove legend boxes that are not used
    for lr, lb in zip(fmeta['lr'], fmeta['lb']):
        # Get all legend objects from lr
        objs_lr = [c for c in lr.get_children() if isinstance(c, legend.Legend)]
        objs_lb = [c for c in lb.get_children() if isinstance(c, legend.Legend)]
        if len(objs_lr) == 0 and lr in fmeta['fig'].axes:
            fmeta['fig'].delaxes(lr)
        if len(objs_lb) == 0 and lb in fmeta['fig'].axes:
            fmeta['fig'].delaxes(lb)

    for lc in fmeta['lc']:
        objs_lc = [c for c in lc.get_children() if isinstance(c, legend.Legend)]
        if len(objs_lc) == 0 and lc in fmeta['fig'].axes:
            fmeta['fig'].delaxes(lc)


def get_markerlist():
    return itertools.cycle(('^', 'v', '<', '>', 'o'))


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
