from matplotlib.ticker import LogLocator, \
    LogFormatter, LogFormatterExponent, LogFormatterSciNotation, \
    LogFormatterMathtext, NullFormatter, MultipleLocator, MaxNLocator
import numpy as np

# mplstyle = '../common/stylesheets/prb.mplstyle'
mplstyle = '../common/stylesheets/slack.mplstyle'
legendoutside = False
legendcollect = False


def get_meta(plotdir):
    meta = {
        'common': {
            'include': {
                # 'L': ['L_8', 'L_12', 'L_16', 'L_20', 'L_24'],
                # 'L': ['L_24', 'L_32'],
                # 'x': ['x_0.5000', 'x_1.0000'],
                # 'f': ['f_0.0125', 'f_0.0250', 'f_0.0500','f_0.0750'],
                # 'f': ['f_0.1000', 'f_0.1250', 'f_0.1500','f_0.1750','f_0.2000'],
                # 'f': ['f_0.1500','f_0.2000'],
                # 'f': ['f_0.1500', 'f_0.2500', 'f_0.3000'],
                # 'u': ['u_2', 'u_3'],
                # 'u': ['u_4'],
                # 'w': ['w[+1.0000_+0.2500_+0.1000]', 'w[+1.0000_+0.5000_+0.1000]']
            },
        },

        'ent': {
            'groupname': 'measurements',
            'colname': 'entanglement_entropy_midchain',
            'normpage': False,
            'titlename': 'Entanglement Entropy',
            'ylabel': '$\langle \langle S_\mathrm{E}(L/2)\\rangle \\rangle$',
            'xlabel': '$t$',
            # 'xmaloc': LogLocator(base=10, numticks=10, numdecs=32),
            # 'xmiloc': LogLocator(base=10, numticks=10, subs=(.1, .2, .3, .4, .5, .6, .7, .8, .9)),
            'xmafmt': LogFormatterMathtext(),
            'plotprefix': 'SE',
            'plotdir': plotdir,
            'findsaturation': True,  # Instead of taking the last value, take the average of the plateau
            'findloglogwindow': True,
            'markloglogwindow': True,
            'timeloglevel': 1,
            'mplstyle': mplstyle,
            'legendcols': ['f', 'x', 'num', 'bmax:.0f', 'bavg:.0f', 'tsim'],  # Choose 'num', 'bmax','tsim'
            'legendoutside': legendoutside,
            'legendcollect': legendcollect,
            'legendlocation': 'lower right',
        },
        'num1': {
            'groupname': 'measurements',
            'colname': 'number_entropy_midchain',
            'normpage': False,
            'titlename': 'Number Entropy',
            'ylabel': '$\langle S_N(L/2)\\rangle$',
            'yformat': '%.2f',
            'plotprefix': 'SN',
            'plotdir': plotdir,
            # 'ymin': 0.41,
            # 'ymax': 0.45,
            # 'xmin': 1,
            'findsaturation': True,  # Instead of taking the last value, take the average of the plateau
            'findloglogwindow': True,
            'markloglogwindow': True,
            'fitloglogwindow': True,
            'timeloglevel': 1,
            'mplstyle': mplstyle,
            'legendcols': ['f', 'x', 'num', 'bmax:.0f', 'bavg:.0f', 'tsim'],  # Choose 'num', 'bmax','tsim'
            'legendoutside': legendoutside,
            'legendcollect': legendcollect,
            'legendlocation': 'lower right',
        },
        'num2': {
            'groupname': 'measurements',
            'colname': 'number_entropy_midchain',
            'normpage': False,
            'titlename': 'Number Entropy',
            'ylabel': '$\langle S_N(L/2)\\rangle$',
            'yformat': '%.2f',
            'plotprefix': 'SN',
            'plotdir': plotdir,
            # 'ymin': 0.40,
            # 'ymax': 0.45,
            'findsaturation': True,  # Instead of taking the last value, take the average of the plateau
            'findloglogwindow': True,
            'markloglogwindow': True,
            'fitloglogwindow': True,
            'timeloglevel': 2,
            'mplstyle': mplstyle,
            'legendcols': ['f', 'x', 'num', 'bmax:.0f', 'bavg:.0f', 'tsim'],  # Choose 'num', 'bmax','tsim'
            'legendoutside': legendoutside,
            'legendcollect': legendcollect,
            'legendlocation': 'lower right',
        },
        'numH1': {  # Hartley number entropy
            'groupname': 'measurements',
            'colname': 'hartley_number_entropy_midchain',
            'normpage': False,
            'titlename': ' HartleyNumber Entropy',
            'ylabel': '$\langle S_H(L/2)\\rangle$',
            'yformat': '%.2f',
            # 'ymin': 1.26,
            'plotprefix': 'SH',
            'plotdir': plotdir,
            'findsaturation': True,  # Instead of taking the last value, take the average of the plateau
            'findloglogwindow': True,
            'markloglogwindow': True,
            'fitloglogwindow': True,
            'timeloglevel': 1,
            'mplstyle': mplstyle,
            'legendcols': ['f', 'x', 'num', 'bmax:.0f', 'bavg:.0f', 'tsim'],  # Choose 'num', 'bmax','tsim'
            'legendoutside': legendoutside,
            'legendcollect': legendcollect,
            'legendlocation': 'lower right',
        },
        'numH2': {  # Hartley number entropy
            'groupname': 'measurements',
            'colname': 'hartley_number_entropy_midchain',
            'normpage': False,
            'titlename': ' HartleyNumber Entropy',
            'ylabel': '$\langle S_H(L/2)\\rangle$',
            'yformat': '%.2f',
            # 'ymin': 1.27,
            'plotprefix': 'SH',
            'plotdir': plotdir,
            'findsaturation': True,  # Instead of taking the last value, take the average of the plateau
            'findloglogwindow': True,
            'markloglogwindow': True,
            'fitloglogwindow': True,
            'timeloglevel': 2,
            'mplstyle': mplstyle,
            'legendcols': ['f', 'x', 'num', 'bmax:.0f', 'bavg:.0f', 'tsim'],  # Choose 'num', 'bmax','tsim'
            'legendoutside': legendoutside,
            'legendcollect': legendcollect,
            'legendlocation': 'lower right',
        },
        'numa': {
            'groupname': 'measurements',
            'colname': 'number_entropy_midchain',
            'normpage': False,
            'titlename': 'Number Entropy Approach to saturation',
            'ylabel': '$|S_\mathrm{N}(L/2, t) - \\tilde{S}_\mathrm{N}(L/2)|$',
            # 'yformat': '%.2f',
            'yscale': 'log',
            'plotprefix': 'SNA',
            'plotdir': plotdir,
            'findsaturation': True,  # Instead of taking the last value, take the average of the plateau
            'findloglogwindow': True,
            'markloglogwindow': True,
            'fitloglogwindow': False,
            'plotsatapproach': True,
            'timeloglevel': 1,
            'mplstyle': mplstyle,
            'legendcols': ['f', 'x', 'num', 'bmax:.0f', 'bavg:.0f', 'tsim'],  # Choose 'num', 'bmax','tsim'
            'legendoutside': legendoutside,
            'legendcollect': legendcollect,
            'legendlocation': 'lower left',
        },
        'numHa': {
            'groupname': 'measurements',
            'colname': 'hartley_number_entropy_midchain',
            'normpage': False,
            'titlename': 'Hartley Number Entropy Approach to saturation',
            'ylabel': '$|S_\mathrm{H}(L/2, t) - \\tilde{S}_\mathrm{H}(L/2)|$',
            # 'yformat': '%.2f',
            'yscale': 'log',
            'plotprefix': 'SHA',
            'plotdir': plotdir,
            'findsaturation': True,  # Instead of taking the last value, take the average of the plateau
            'findloglogwindow': True,
            'markloglogwindow': True,
            'fitloglogwindow': False,
            'plotsatapproach': True,
            'timeloglevel': 1,
            'mplstyle': mplstyle,
            'legendcols': ['f', 'x', 'num', 'bmax:.0f', 'bavg:.0f', 'tsim'],  # Choose 'num', 'bmax','tsim'
            'legendoutside': legendoutside,
            'legendcollect': legendcollect,
            'legendlocation': 'lower left',
        },
        'nument1': {
            'groupname': 'measurements',
            'colname': ['number_entropy_midchain', 'entanglement_entropy_midchain'],
            'normpage': False,
            'titlename': 'Entanglement and Number Entropies',
            'ylabel': '$\langle S(L/2)\\rangle$',
            'plotprefix': 'SESN',
            'plotdir': plotdir,
            'findsaturation': True,  # Instead of taking the last value, take the average of the plateau
            'findloglogwindow': True,
            'markloglogwindow': True,
            'fitloglogwindow': True,
            'zoomloglogwindow': {
                'colnum': [0],  # Which columns (colname) to zoom in on
                # 'pos': [0.03, 0.6, 0.40, 0.40], # Positon of the inset, x0 y0 width height
                'pos': [0.5, 0.35, 0.35, 0.35],  # Positon of the inset, x0 y0 width height
                'coords': [None, None, None, None],
                # These zoom limits x1,x2,y1,y2, must be set by finding the maximum log log window
                'legendtitle': '$S_N$'
            },
            'timeloglevel': 1,
            'linestyle': ['solid', 'dashed'],
            'mplstyle': mplstyle,
            'legendcols': ['f', 'x', 'num', 'bmax:.0f', 'bavg:.0f', 'tsim'],  # Choose 'num', 'bmax','tsim'
            'legendoutside': legendoutside,
            'legendcollect': legendcollect,
            'legendlocation': 'upper left',
        },
        'nument2': {
            'groupname': 'measurements',
            'colname': ['number_entropy_midchain', 'entanglement_entropy_midchain'],
            'normpage': False,
            'titlename': 'Entanglement and Number Entropies',
            'ylabel': '$\langle S(L/2)\\rangle$',
            'plotprefix': 'SESN',
            'plotdir': plotdir,
            'findsaturation': True,  # Instead of taking the last value, take the average of the plateau
            'findloglogwindow': True,
            'markloglogwindow': True,
            'fitloglogwindow': True,
            'zoomloglogwindow': {
                'colnum': [0],  # Which columns (colname) to zoom in on
                'pos': [0.03, 0.6, 0.32, 0.32],  # Positon of the inset, x0 y0 width height
                'coords': [None, None, None, None],
                # These zoom limits x1,x2,y1,y2, must be set by finding the maximum log log window
                'legendtitle': '$S_N$'
            },
            'timeloglevel': 1,
            'linestyle': ['solid', 'dashed'],
            'mplstyle': mplstyle,
            'legendcols': ['f', 'x', 'num', 'bmax:.0f', 'bavg:.0f', 'tsim'],  # Choose 'num', 'bmax','tsim'
            'legendoutside': legendoutside,
            'legendcollect': legendcollect,
            'legendlocation': 'lower right',
        },

        'chi': {
            'groupname': 'measurements',
            'colname': ['bond_mid', 'bond_dimension_midchain'],
            'titlename': 'Bond Dimension',
            'ylabel': '$\langle\chi\\rangle$',
            'plotprefix': 'chi',
            'plotdir': plotdir,
            'findsaturation': True,  # Instead of taking the last value, take the average of the plateau
            'findloglogwindow': False,
            'timeloglevel': 1,
            'mplstyle': mplstyle,
            'legendcols': ['x', 'num', 'bavg:.0f', 'tsim'],  # Choose 'num', 'bmax','tsim'
            'legendoutside': legendoutside,
            'legendcollect': legendcollect,
            'legendlocation': 'best',
        },
        'tsim': {
            'groupname': 'measurements',
            'colname': 'algorithm_time',
            'titlename': 'Simulation Time',
            'ylabel': 'Avg. Simulation Time [min]',
            'yscale': 'log',
            'plotprefix': 'tsim',
            'plotdir': plotdir,
            'normalize': 60,
            'realizations': 600,
            'findsaturation': False,  # Instead of taking the last value, take the average of the plateau
            'findloglogwindow': False,
            'timeloglevel': 1,
            'mplstyle': mplstyle,
            'legendcols': ['num'],  # Choose 'num', 'bmax','tsim'
            'legendoutside': legendoutside,
            'legendcollect': legendcollect,
        },
        'titr': {
            'groupname': 'measurements',
            'colname': 'algorithm_time',
            'titlename': 'Iteration Time',
            'ylabel': 'Avg Iteration Time [min]',
            'ydiff': True,  # Plot the difference in y instead.
            'yscale': 'log',
            'plotprefix': 'titr',
            'plotdir': plotdir,
            'normalize': 60,
            'realizations': 600,
            'findsaturation': False,  # Instead of taking the last value, take the average of the plateau
            'findloglogwindow': False,
            'timeloglevel': 1,
            'mplstyle': mplstyle,
            'legendcols': ['num'],  # Choose 'num', 'bmax','tsim'
            'legendoutside': legendoutside,
            'legendcollect': legendcollect,
        },
        'trn': {
            'groupname': 'measurements',
            'colname': 'truncation_error',
            'normpage': False,
            'titlename': 'Truncation Error',
            'ylabel': '$\epsilon$',
            'yscale': 'log',
            'sharex': 'all',
            'sharey': 'all',
            # MultipleLocator(tick_spacing)
            # 'xmaloc': LogLocator(base=10, numticks=99, numdecs=16),
            # 'xmiloc': LogLocator(base=10, numticks=99, subs=(.1, .2, .3, .4, .5, .6, .7, .8, .9)),
            # 'xmafmt': LogFormatter(base=10.0, labelOnlyBase=True, minor_thresholds=(1.0, 0.4), linthresh=None),
            # 'xmaloc': LogLocator(numticks=9, numdecs=32),
            # 'xmiloc': LogLocator(subs=(.1, .2, .3, .4, .5, .6, .7, .8, .9)),
            'xmafmt': LogFormatterMathtext(labelOnlyBase=True),

            'plotprefix': 'eps',
            'plotdir': plotdir,
            'findsaturation': False,  # Instead of taking the last value, take the average of the plateau
            'findloglogwindow': False,
            'markloglogwindow': False,
            'timeloglevel': 1,
            'mplstyle': mplstyle,
            'legendcols': ['x', 'w', 'num', 'bmax:.0f', 'bavg:.0f', 'tsim'],  # Choose 'num', 'bmax','tsim'
            'legendoutside': legendoutside,
            'legendcollect': legendcollect,
            'legendlocation': 'best',
        },
        'lbit': {
            'groupname': 'lbits',
            'dsetname': 'decay',
            'titlename': 'l-bit',
            'ylabel': '$\langle \langle \mathrm{l-bit} \\rangle\\rangle$ ',
            'xlabel': "$|i-j|$",
            'yscale': 'log',
            'ynopos': 'clip',
            'plotprefix': 'lbit',
            'plotdir': plotdir,
            'mplstyle': mplstyle,
            'xmax': 12,
            'ymin': 1e-15,
            'legendcols': ['num'],  # Choose 'num', 'bmax','tsim'
            'legendoutside': legendoutside,
            'legendcollect': legendcollect,
            'legendlocation': 'best',
        },
        'dist-num': {
            'groupname': 'number_entropies',
            'dsetname': 'data',
            'normpage': False,
            'titlename': 'Number Entropy',
            # 'figsize': (3.375, 3.00),
            'ylabel': '$p(S_N(L/2))$',
            'yscale': 'log',
            'yformat': '%.2f',
            'sharex': 'all',
            'sharey': 'all',
            'tidx': 'window',  # Time indices for which to plot the distribution
            'bins': 40,
            'xlabel': '$S_N$',
            'density': True,
            'plotprefix': 'SN',
            'plotdir': plotdir,
            # 'ymin': 0.41,
            # 'ymax': 0.45,
            # 'xmin': 1,
            'findsaturation': True,  # Instead of taking the last value, take the average of the plateau
            'findloglogwindow': True,
            'markloglogwindow': True,
            'fitloglogwindow': True,
            # 'timeloglevel': 1,
            'mplstyle': mplstyle,
            # 'legendcols': ['f', 'x', 'num', 'bmax:.0f', 'bavg:.0f', 'tsim', 't:.1e'],  # Choose 'num', 'bmax','tsim'
            'legendcols': ['L', 'f', 'w', 't:<8.1e'],  # Choose 'num', 'bmax','tsim'
            # 'legendcols': ['f', 'x', 'w', 'num', 't:<8.1e'],  # Choose 'num', 'bmax','tsim'
            'legendoutside': False,
            'legendcollect': False,
            'legendlocation': 'lower left',
        },
        'divg-ent': {  # Distribution of infinite time averaged entropy
            'groupname': 'entanglement_entropies',
            'dsetname': 'data',
            'normpage': False,
            'titlename': 'Entanglement Entropy',
            # 'figsize': (3.375, 3.00),
            'ylabel': '$p(S_E^\infty(L/2))$',
            'yscale': 'log',
            'yformat': '%.2f',
            'sharex': 'all',
            'sharey': 'all',
            'xmin': 0.0,
            'xmax': 1.5,
            'tidx': 'window',  # Time indices for which to plot the distribution
            'bins': 40,
            'xlabel': '$S_E^\infty$',
            'density': True,
            'plotprefix': 'SE',
            'plotdir': plotdir,
            # 'ymin': 0.41,
            # 'ymax': 0.45,
            # 'xmin': 1,
            'findsaturation': True,  # Instead of taking the last value, take the average of the plateau
            'findloglogwindow': True,
            'markloglogwindow': True,
            'fitloglogwindow': True,
            # 'timeloglevel': 1,
            'mplstyle': mplstyle,
            # 'legendcols': ['f', 'x', 'num', 'bmax:.0f', 'bavg:.0f', 'tsim', 't:.1e'],  # Choose 'num', 'bmax','tsim'
            # 'legendcols': ['f', 'w'],  # Choose 'num', 'bmax','tsim'
            'legendcols': ['L', 'f', 'x', 'w', 'num'],  # Choose 'num', 'bmax','tsim'
            'legendoutside': False,
            'legendcollect': False,
            'legendlocation': 'lower left',
        },
        'divg-num': {  # Distribution of infinite time averaged entropy
            'groupname': 'number_entropies',
            'dsetname': 'data',
            'normpage': False,
            'titlename': 'Number Entropy',
            # 'figsize': (3.375, 3.00),
            'ylabel': '$p(S_N^\infty(L/2))$',
            'yscale': 'log',
            'yformat': '%.2f',
            'sharex': 'all',
            'sharey': 'all',
            'xmin': 0.0,
            'xmax': 1.5,
            'tidx': 'window',  # Time indices for which to plot the distribution
            'bins': 40,
            'xlabel': '$S_N^\infty$',
            'density': True,
            'plotprefix': 'SN',
            'plotdir': plotdir,
            # 'ymin': 0.41,
            # 'ymax': 0.45,
            # 'xmin': 1,
            'findsaturation': True,  # Instead of taking the last value, take the average of the plateau
            'findloglogwindow': True,
            'markloglogwindow': True,
            'fitloglogwindow': True,
            # 'timeloglevel': 1,
            'mplstyle': mplstyle,
            # 'legendcols': ['f', 'x', 'num', 'bmax:.0f', 'bavg:.0f', 'tsim', 't:.1e'],  # Choose 'num', 'bmax','tsim'
            # 'legendcols': ['f', 'w'],  # Choose 'num', 'bmax','tsim'
            'legendcols': ['L', 'f', 'x', 'w', 'num'],  # Choose 'num', 'bmax','tsim'
            'legendoutside': False,
            'legendcollect': False,
            'legendlocation': 'lower left',
        },
        'tavg-ent': {
            'groupname': 'entanglement_entropies',
            'dsetname': 'data',
            'normpage': False,
            'titlename': 'Infinite time average of Entanglement Entropy',
            # 'figsize': (3.375, 3.00),
            'ylabel': '$\langle \langle S_\mathrm{E}^\infty(L/2) \\rangle \\rangle$',
            'yformat': '%.2f',
            'xmaloc': MaxNLocator(integer=True),
            'sharex': 'all',
            'sharey': 'all',
            'xlabel': '$L$',
            'density': True,
            'plotprefix': 'SE',
            'plotdir': plotdir,
            # 'ymin': 0.41,
            # 'ymax': 0.45,
            # 'xmin': 1,
            'findsaturation': True,  # Instead of taking the last value, take the average of the plateau
            'findloglogwindow': True,
            'markloglogwindow': True,
            'fitloglogwindow': True,
            # 'timeloglevel': 1,
            'mplstyle': mplstyle,
            # 'legendcols': ['f', 'x', 'num', 'bmax:.0f', 'bavg:.0f', 'tsim', 't:.1e'],  # Choose 'num', 'bmax','tsim'
            'legendcols': ['L', 'f', 'x', 'num'],  # Choose 'num', 'bmax','tsim'
            'legendoutside': False,
            'legendcollect': False,
            'legendlocation': 'lower left',
        },
        'tavg-num': {
            'groupname': 'number_entropies',
            'dsetname': 'data',
            'normpage': False,
            'titlename': 'Infinite time average of Number Entropies',
            # 'figsize': (3.375, 3.00),
            'ylabel': '$\langle \langle S_\mathrm{N}^\infty(L/2) \\rangle \\rangle$',
            'yformat': '%.2f',
            'xmaloc': MaxNLocator(integer=True),
            'sharex': 'all',
            'sharey': 'all',
            'xlabel': '$L$',
            'density': True,
            'plotprefix': 'SN',
            'plotdir': plotdir,
            # 'ymin': 0.41,
            # 'ymax': 0.45,
            # 'xmin': 1,
            'findsaturation': True,  # Instead of taking the last value, take the average of the plateau
            'findloglogwindow': True,
            'markloglogwindow': True,
            'fitloglogwindow': True,
            # 'timeloglevel': 1,
            'mplstyle': mplstyle,
            # 'legendcols': ['f', 'x', 'num', 'bmax:.0f', 'bavg:.0f', 'tsim', 't:.1e'],  # Choose 'num', 'bmax','tsim'
            'legendcols': ['L', 'f', 'x', 'num'],  # Choose 'num', 'bmax','tsim'
            'legendoutside': False,
            'legendcollect': False,
            'legendlocation': 'lower left',
        },
    }

    return meta