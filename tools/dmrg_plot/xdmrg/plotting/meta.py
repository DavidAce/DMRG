import numpy as np
from matplotlib.ticker import NullFormatter, FormatStrFormatter, ScalarFormatter, LogLocator, LinearLocator

# mplstyle = './src/plotting/stylesheets/prb.mplstyle'
mplstyle = './src/plotting/stylesheets/slack.mplstyle'
legendoutside = False
legendcollect = False


def get_meta(plotdir):
    meta = {
        'common': {
            'include': {
                # 'L': ['L_16', 'L_18'],
                # 'g': ['g_0.0100','g_0.0150'],
                # 'd': ('d_+0.0000'),
                # 'l': ['g_0.0100', 'g_0.0200', 'g_0.0300','g_0.0400'],
                # 'algo'  : ['xDMRG', 'ed-e0.00'],
                # 'state' : ['state_0', 'state_1', 'states'],
                # 'point' : ['fes', 'bondpoint', 'tables', 'measurements'], # measurements correspond to ed data
            },
        },

        'mid-ent': {
            'titlename': 'Midchain von Neumann Entropy',
            'load': 'measurements',
            'ygrp': 'measurements',
            'ycol': 'entanglement_entropy_midchain',
            # 'yformat'    : '%.2f',
            'ylabel': '$\langle\langle S_\mathrm{vN}(L/2) \\rangle\\rangle /L$',
            'yscale': 'linear',
            'ydiv': 'L',  # Divide y-data by L
            # 'xscale'     : 'log',
            'sharey': 'all',
            # 'yticks'     : [0.01, 0.1,1,10],
            'xlim': {'d': (-10, 10)},
            # 'ylim'       : (0.01,1.5),

            # 'xlabel'     : '$\Delta = \log{\\bar J} - \log{\\bar h}$',
            'plotprefix': 'SE',
            'plotdir': plotdir,
            'mplstyle': mplstyle,
            'legendcols': ['L', 'l:.2f', 'd', 'num', 'time:.2f'],  # Choose 'num', 'time', 'algo', 'state'
            'legendoutside': False,
            'legendcollect': True,
            'inc': {
                # 'algo' :  ['xDMRG'],
                'algo': ['xDMRG', 'ed-e0.00'],
                # 'state': ['state_0', 'state_1'],
                'state': ['state_0', 'state_1', 'states'],
                'point': ['tables', 'measurements'],  # measurements correspond to ed data
            },
        },

        'chain-ent': {
            'titlename': 'von Neumann Entropy vs Site',
            'load': 'entanglement_entropies',
            'ygrp': 'entanglement_entropies',
            'ycol': 'L_',  # Matches columns L_0, L_1 ...
            'ytag': 'SvN',  # Used in plot filename
            'ytitle': 'von Neumann Entropy',
            'ylabel': '$\langle\langle S_\mathrm{vN} \\rangle\\rangle$',
            'yformat': '%.2f',
            'yfit': 'S_vs_l_effective_central_charge',
            # 'yscale': 'log',
            'sharey': 'none',
            # 'yticks': [1, 2, 4, 6, 8, 10, 12],
            'ymajor_formatter': FormatStrFormatter('%.2f'),
            'yminor_formatter': NullFormatter(),

            'xgrp': None,
            'xcol': None,
            'xtag': 'site',
            'xtitle': 'Site',
            'xlabel': '$l$',
            # 'xformat': '%.0f',
            'xfold': True,  # Mirror the x-axis around the center value, e.g. for S vs L plots

            'xscale': 'log',
            'sharex': 'none',
            # 'xticks': [1, 2, 4, 6, 8, 10, 12],
            'xmajor_formatter': FormatStrFormatter('%.0f'),
            'xminor_formatter': FormatStrFormatter('%.0f'),
            # 'xticks': [32, 64, 128, 256],

            # 'yscale': 'log',
            # 'xscale': 'log',

            # 'yticks': [0.01, 0.1, 1, 10],
            # 'xmin': 0,
            # 'xmax': 64,

            # 'xlabel'     : '$\Delta = \log{\\bar J} - \log{\\bar h}$',
            'plotprefix': 'SE',
            'plotdir': plotdir,
            'mplstyle': mplstyle,
            'legendcols': ['L', 'l:.2f', 'd', 'num', 'time:.2f'],  # Choose 'num', 'b', 'x','time', 'algo', 'state'
            'legendoutside': False,  # Put legend inside or outside each subplot
            'legendcollect': True,  # Collect equal columns into a single legend outside
            'inc': {
                # 'algo' :  ['xDMRG'],
                'algo': ['xDMRG', 'ed-e0.00'],
                # 'state': ['state_0', 'state_1'],
                'state': ['state_0', 'state_1', 'states'],
                'point': ['tables', 'measurements'],  # measurements correspond to ed data
            },
        },

        'chain-bond': {
            'titlename': 'Average Bond Dimension vs Site',
            'load': 'bond_dims',
            'ygrp': 'bond_dims',
            'ycol': 'L_',  # Matches columns L_0, L_1 ...
            'ytag': 'bond',  # Used in plot filename
            'ytitle': 'Bond Dimension',
            'ylabel': '$\langle\langle \chi \\rangle\\rangle$',
            'yformat': '%.2f',

            'xgrp': None,
            'xcol': None,
            'xtag': 'site',
            'xtitle': 'Site',
            'xlabel': '$l$',
            'xformat': '%.0f',
            'xticks': [1, 2, 4, 6, 8, 10, 12],

            # 'yscale': 'log',
            # 'xscale': 'log',

            # 'yticks': [0.01, 0.1, 1, 10],
            # 'xmin': 0,
            # 'xmax': 64,

            # 'xlabel'     : '$\Delta = \log{\\bar J} - \log{\\bar h}$',
            'plotprefix': 'Chi',
            'plotdir': plotdir,
            'mplstyle': mplstyle,
            'legendcols': ['L', 'l:.2f', 'd', 'num', 'time:.2f'],  # Choose 'num', 'b', 'x','time', 'algo', 'state'
            'legendoutside': False,  # Put legend inside or outside each subplot
            'legendcollect': False,  # Collect equal columns into a single legend outside
            'inc': {
                # 'algo' :  ['xDMRG'],
                'algo': ['xDMRG', 'ed-e0.00'],
                # 'state': ['state_0', 'state_1'],
                'state': ['state_0', 'state_1', 'states'],
                'point': ['tables', 'measurements'],  # measurements correspond to ed data
            },
        },

        'fes-ent': {
            'titlename': 'von Neumann Entropy vs Bond dimension limit',
            'load': ['measurements', 'status'],

            'ygrp': 'measurements',
            'ycol': 'entanglement_entropy_midchain',
            'ytag': 'SvN',
            'ytitle': 'von Neumann Entropy',
            'ylabel': '$\langle\langle S_\mathrm{vN} \\rangle\\rangle$',
            'yformat': '%.2f',
            # 'ymin': 0,
            # 'ymax': 64,
            # 'yscale': 'log',
            # 'yticks': [0.01, 0.1, 1, 10],

            'xgrp': 'status',
            'xcol': 'bond_lim',
            'xtag': 'bond',
            'xtitle': 'Bond dimension limit',
            'xlabel': '$\chi_\\mathrm{lim}$',
            'xformat': '%.0f',
            'xscale': 'log',
            # 'xticks': [32, 64, 128, 256],
            # 'xmin': 0,
            # 'xmax': 64,

            # 'xlabel'     : '$\Delta = \log{\\bar J} - \log{\\bar h}$',
            'plotprefix': 'FES_SE',
            'plotdir': plotdir,
            'mplstyle': mplstyle,
            'legendcols': ['L', 'l:.2f', 'd', 'num', 'time:.2f'],  # Choose 'num', 'b', 'x','time', 'algo', 'state'
            'legendoutside': False,  # Put legend inside or outside each subplot
            'legendcollect': False,  # Collect equal columns into a single legend outside
            'inc': {
                # 'algo' :  ['xDMRG'],
                'algo': ['xDMRG', 'ed-e0.00'],
                # 'state': ['state_0', 'state_1'],
                'state': ['state_0', 'state_1', 'states'],
                'point': ['fes', 'bondpoint', 'measurements'],  # measurements correspond to ed data
            },
        },

        'hist-var': {
            'load': 'measurements',
            'ygrp': 'measurements',
            'ycol': 'energy_variance',
            'titlename': 'Energy Variance',
            # 'yformat': '%.2f',
            'ylabel': '$\mathrm{Var}H$',
            'yscale': 'log',
            'xscale': 'log',
            # 'yticks': [0.01, 0.1, 1, 10],
            # 'xticks': [32, 64, 128, 256],
            'bins': np.logspace(start=-18, stop=0, num=64, base=10, endpoint=True),
            'xmin': 1e-18,
            'xmax': 1,

            # 'xlabel'     : '$\Delta = \log{\\bar J} - \log{\\bar h}$',
            'plotprefix': 'Var',
            'plotdir': plotdir,
            'mplstyle': mplstyle,
            'legendcols': ['L', 'l:.2f', 'd:.2f', 'num', 'time:.2f'],  # ['num', 'time', 'b', 'x'],  # Choose 'num', 'b', 'x','time', 'algo', 'state'
            'legendoutside': False,
            'legendcollect': True,
            'inc': {
                # 'algo' :  ['xDMRG'],
                'algo': ['xDMRG', 'ed-e0.00'],
                # 'state': ['state_0', 'state_1'],
                'state': ['state_0', 'state_1', 'states'],
                # 'point': ['fes', 'bondpoint', 'measurements'],  # measurements correspond to ed data
                'point': ['tables', 'measurements'],  # measurements correspond to ed data
            },
        },
        'time': {
            'load': 'measurements',
            'colname': 'algorithm_time',
            'titlename': 'Simulation Time',
            # 'yformat': '%.2f',
            'ylabel': '$\mathrm{Time } [s]$',
            'yscale': 'log',
            # 'yticks': [0.01, 0.1, 1, 10],
            # 'xmin': -10,
            # 'xmax': 10,

            # 'xlabel'     : '$\Delta = \log{\\bar J} - \log{\\bar h}$',
            'plotprefix': 'Time',
            'plotdir': plotdir,
            'mplstyle': mplstyle,
            'legendcols': ['b', 'num'],  # Choose 'num', 'b', 'x','time', 'algo', 'state'
            'legendoutside': False,
            'inc': {
                # 'algo' :  ['xDMRG'],
                'algo': ['xDMRG', 'ed-e0.00'],
                # 'state': ['state_0', 'state_1'],
                'state': ['state_0', 'state_1', 'states'],
                'point': ['tables', 'measurements'],  # measurements correspond to ed data
            },
        },
    }
    return meta
