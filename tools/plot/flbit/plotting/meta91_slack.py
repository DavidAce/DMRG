from matplotlib.ticker import LogLocator, \
    LogFormatter, LogFormatterExponent, LogFormatterSciNotation, \
    LogFormatterMathtext, NullFormatter, MultipleLocator, MaxNLocator
import numpy as np
from pathlib import Path

# mplstyle = '../common/stylesheets/prb.mplstyle'
mplstyle = '../common/stylesheets/slack.mplstyle'
prb = 'prb' in mplstyle



def get_meta(plotdir):

    default = {
        'box_aspect': 1,
        'subspec_title': False,
        'figspec_title': False,
        'legendoutside' : False,
        'legendcollect' : False
    }

    meta = {
        'common': {
            'include': {
                # 'L': ['L_32'],
                # 'L': ['L_8', 'L_12', 'L_16', 'L_20'],
                # 'L': ['L_12'],
                # 'x': ['x_0.5000', 'x_1.0000'],
                # 'f': ['f_0.0125', 'f_0.0250', 'f_0.0500','f_0.0750'],
                # 'f': ['f_0.1000', 'f_0.1250', 'f_0.1500','f_0.1750','f_0.2000'],
                # 'f': ['f_0.5000'],
                # 'f': ['f_0.1500', 'f_0.2500', 'f_0.3000'],
                # 'f': ['f_0.4000', 'f_0.5000', 'f_0.6000'],
                # 'f': ['f_0.5000', 'f_1.0000', 'f_2.0000', 'f_3.0000', 'f_4.0000', 'f_5.0000'],
                # 'u': ['u_2', 'u_3'],
                # 'f': ['f_0.4500'],
                # 'u': ['u_5', 'u_6'],
                # 'u': ['u_5'],
                # 'w': ['w[+1.0000_+0.2500_+0.1000]', 'w[+1.0000_+0.5000_+0.1000]']
            },
            'include_v3': {
                # 'L': [36],
                # 'L': [12],
                # 'L': ['L_8', 'L_12', 'L_16', 'L_20'],
                # 'L': ['L_12'],
                # 'x': ['x_0.5000', 'x_1.0000'],
                # 'f': ['f_0.0125', 'f_0.0250', 'f_0.0500','f_0.0750'],
                # 'f': ['f_0.1000', 'f_0.1250', 'f_0.1500','f_0.1750','f_0.2000'],
                # 'f': [0.5000],
                # 'f': ['f_0.1500', 'f_0.2500', 'f_0.3000'],
                # 'f': ['f_0.4000', 'f_0.5000', 'f_0.6000'],
                # 'f': ['f_0.5000', 'f_1.0000', 'f_2.0000', 'f_3.0000', 'f_4.0000', 'f_5.0000'],
                # 'u': ['u_2', 'u_3'],
                # 'f': ['f_0.4500'],
                # 'u': ['u_5', 'u_6'],
                # 'u': [32],
                # 'w': ['w[+1.0000_+0.2500_+0.1000]', 'w[+1.0000_+0.5000_+0.1000]']
                # 'tgw8': ['ID'],
                # 'cgw8': ['EX'],
                # 'u': [8,16,24,32,40,48,56,64],
                # 'u': [8, 16, 32, 64, 80],
                # 'ubond': [256],
            },
        },

        'ent': {
            'default' : default,
            'axtitle': False,
            'groupname': 'measurements',
            'colname': 'entanglement_entropy',
            'normpage': False,
            'titlename': 'Entanglement Entropy',
            'ylabel': '$\langle \langle S(L/2)\\rangle \\rangle$',
            'xlabel': '$t$',
            'yscale': 'log',
            # 'xmaloc': LogLocator(base=10, numticks=10, numdecs=32),
            # 'xmiloc': LogLocator(base=10, numticks=10, subs=(.1, .2, .3, .4, .5, .6, .7, .8, .9)),
            # 'ymax': 1e10,
            'xmafmt': LogFormatterMathtext(),
            'plotprefix': 'SE',
            'plotdir': Path(plotdir, Path(mplstyle).stem),
            'findsaturation': True,  # Instead of taking the last value, take the average of the plateau
            'findloglogwindow': True,
            'markloglogwindow': True,
            'timeselection': 'lnt',
            'mplstyle': mplstyle,
            'legendcols': ['f', 'bavg:.0f', 'tsim'],  # Choose 'num', 'bmax','tsim'
            # 'legendcols': ['f', 'x', 'num', 'bmax:.0f', 'bavg:.0f', 'tsim'],  # Choose 'num', 'bmax','tsim'
            
            
            'legendlocation': 'lower right',
        },
        'num1': {
            'default' : default,
            'axtitle': False,
            'groupname': 'measurements',
            'colname': 'number_entropy',
            'normpage': False,
            'titlename': 'Number Entropy',
            'ylabel': '$\langle \langle S_N(L/2)\\rangle \\rangle$',
            'yformat': '%.3f',
            # 'yscale': 'log',
            'plotprefix': 'SN',
            'plotdir': Path(plotdir, Path(mplstyle).stem),
            # 'ymin': 0.120,
            # 'ymax': 0.325,
            # 'xmin': 5e-1,
            # 'xmax': 5e+3,
            'findsaturation': True,  # Instead of taking the last value, take the average of the plateau
            'findloglogwindow': True,
            'markloglogwindow': True,
            'fitloglogwindow': True,
            'shadederror': False,
            'timeselection': 'lnt',
            'mplstyle': mplstyle,
            # 'legendcols': ['f', 'x', 'num', 'bmax:.0f', 'bavg:.0f', 'tsim'],  # Choose 'num', 'bmax','tsim'
            'legendcols': ['L', 'f', 'bavg:.0f','tsim'],  # Choose 'num', 'bmax','tsim'
            
            
            'legendlocation': 'lower right',
            # 'legendlocation': (0.01, 0.65),
        },
        'num2': {
            'default' : default,
            'axtitle': False,
            'groupname': 'measurements',
            'colname': 'number_entropy',
            'normpage': False,
            'titlename': 'Number Entropy',
            'filter':{
                # 'f': [0.3,0.5],
                # 'u': [16],
            },
            'ylabel': '$\langle\langle S_N(L/2)\\rangle \\rangle$',
            'yformat': '%.3f',
            'yscale': 'log',
            'plotprefix': 'SN',
            'plotdir': Path(plotdir, Path(mplstyle).stem),
            'sharex': 'all',
            'sharey': 'none',
            'ymin': 0.120,
            'ymax': 0.325,
            'xmin': 0.0,
            'xmax': 2.5,
            'findsaturation': True,  # Instead of taking the last value, take the average of the plateau
            'findloglogwindow': True,
            'markloglogwindow': True,
            'fitloglogwindow': True,
            'fillerror': False,
            'timeselection': 'lnlnt',
            'mplstyle': mplstyle,
            # 'legendcols': ['f', 'x', 'num', 'bmax:.0f', 'bavg:.0f', 'tsim'],  # Choose 'num', 'bmax','tsim'
            'legendcols': ['L', 'f', 'bavg:.0f','tsim'],  # Choose 'num', 'bmax','tsim'
            
            
            'legendlocation': 'lower right',
            # 'legendlocation': (0.01, 0.65),
        },
        'numH1': {  # Hartley number entropy
            'default' : default,
            'axtitle': False,
            'groupname': 'measurements',
            'colname': 'hartley_number_entropy',
            'normpage': False,
            'titlename': ' HartleyNumber Entropy',
            'ylabel': '$\langle S_H(L/2)\\rangle$',
            'yformat': '%.2f',
            # 'ymin': 1.26,
            'plotprefix': 'SH',
            'plotdir': Path(plotdir, Path(mplstyle).stem),
            'findsaturation': True,  # Instead of taking the last value, take the average of the plateau
            'findloglogwindow': True,
            'markloglogwindow': True,
            'fitloglogwindow': True,
            'timeselection': 'lnt',
            'mplstyle': mplstyle,
            'legendcols': ['f', 'bmax:.0f', 'bavg:.0f', 'tsim'],  # Choose 'num', 'bmax','tsim'
            
            
            'legendlocation': 'lower right',
        },
        'numH2': {  # Hartley number entropy
            'default' : default,
            'groupname': 'measurements',
            'colname': 'hartley_number_entropy',
            'normpage': False,
            'titlename': ' HartleyNumber Entropy',
            'ylabel': '$\langle S_H(L/2)\\rangle$',
            'yformat': '%.2f',
            # 'ymin': 1.27,
            'plotprefix': 'SH',
            'plotdir': Path(plotdir, Path(mplstyle).stem),
            'findsaturation': True,  # Instead of taking the last value, take the average of the plateau
            'findloglogwindow': True,
            'markloglogwindow': True,
            'fitloglogwindow': True,
            'timeselection': 'lnlnt',
            'mplstyle': mplstyle,
            'legendcols': ['f', 'bmax:.0f', 'bavg:.0f', 'tsim'],  # Choose 'num', 'bmax','tsim'
            
            
            'legendlocation': 'lower right',
        },
        'numa': {
            'default' : default,
            'groupname': 'measurements',
            'colname': 'number_entropy',
            'normpage': False,
            'titlename': 'Number Entropy Approach to saturation',
            'ylabel': '$|S_\mathrm{N}(L/2, t) - \\tilde{S}_\mathrm{N}(L/2)|$',
            # 'yformat': '%.2f',
            'yscale': 'log',
            'plotprefix': 'SNA',
            'plotdir': Path(plotdir, Path(mplstyle).stem),
            'findsaturation': True,  # Instead of taking the last value, take the average of the plateau
            'findloglogwindow': True,
            'markloglogwindow': True,
            'fitloglogwindow': False,
            'plotsatapproach': True,
            'timeselection': 'lnt',
            'mplstyle': mplstyle,
            'legendcols': ['f', 'bmax:.0f', 'bavg:.0f', 'tsim'],  # Choose 'num', 'bmax','tsim'
            
            
            'legendlocation': 'lower left',
        },
        'numHa': {
            'default' : default,
            'groupname': 'measurements',
            'colname': 'hartley_number_entropy',
            'normpage': False,
            'titlename': 'Hartley Number Entropy Approach to saturation',
            'ylabel': '$|S_\mathrm{H}(L/2, t) - \\tilde{S}_\mathrm{H}(L/2)|$',
            # 'yformat': '%.2f',
            'yscale': 'log',
            'plotprefix': 'SHA',
            'plotdir': Path(plotdir, Path(mplstyle).stem),
            'findsaturation': True,  # Instead of taking the last value, take the average of the plateau
            'findloglogwindow': True,
            'markloglogwindow': True,
            'fitloglogwindow': False,
            'plotsatapproach': True,
            'timeselection': 'lnt',
            'mplstyle': mplstyle,
            'legendcols': ['f', 'bmax:.0f', 'bavg:.0f', 'tsim'],  # Choose 'num', 'bmax','tsim'
            
            
            'legendlocation': 'lower left',
        },
        'nument1': {
            'default' : default,
            'groupname': 'measurements',
            'colname': ['number_entropy', 'entanglement_entropy'],
            'normpage': False,
            'titlename': 'Entanglement and Number Entropies',
            'ylabel': '$\langle S(L/2)\\rangle$',
            'plotprefix': 'SESN',
            'plotdir': Path(plotdir, Path(mplstyle).stem),
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
            'timeselection': 'lnt',
            'linestyle': ['solid', 'dashed'],
            'mplstyle': mplstyle,
            # 'legendcols': ['f', 'x', 'num', 'bmax:.0f', 'bavg:.0f', 'tsim'],  # Choose 'num', 'bmax','tsim'
            'legendcols': ['f'],  # Choose 'num', 'bmax','tsim'
            
            
            'legendlocation': 'upper left',
        },
        'nument2': {
            'default' : default,
            'groupname': 'measurements',
            'colname': ['number_entropy', 'entanglement_entropy'],
            'normpage': False,
            'titlename': 'Entanglement and Number Entropies',
            'ylabel': '$\langle S(L/2)\\rangle$',
            'plotprefix': 'SESN',
            'plotdir': Path(plotdir, Path(mplstyle).stem),
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
            'timeselection': 'lnt',
            'linestyle': ['solid', 'dashed'],
            'mplstyle': mplstyle,
            'legendcols': ['f', 'x', 'num', 'bmax:.0f', 'bavg:.0f', 'tsim'],  # Choose 'num', 'bmax','tsim'
            
            
            'legendlocation': 'lower right',
        },

        'chi': {
            'default' : default,
            'groupname': 'measurements',
            'colname': ['bond_mid'],
            'titlename': 'Midchain Bond Dimension',
            'ylabel': '$\langle\chi\\rangle$',
            'plotprefix': 'chi',
            'plotdir': Path(plotdir, Path(mplstyle).stem),
            'findsaturation': True,  # Instead of taking the last value, take the average of the plateau
            'findloglogwindow': False,
            'timeselection': 'lnt',
            'mplstyle': mplstyle,
            'legendcols': ['bavg:.0f', 'tsim'],  # Choose 'num', 'bmax','tsim'
            
            
            'legendlocation': 'best',
        },
        'tsim': {
            'default' : default,
            'groupname': 'measurements',
            'colname': 'algorithm_time',
            'titlename': 'Simulation Time',
            'ylabel': 'Avg. Simulation Time [min]',
            'yscale': 'log',
            'plotprefix': 'tsim',
            'plotdir': Path(plotdir, Path(mplstyle).stem),
            'normalize': 60,
            'realizations': 600,
            'findsaturation': False,  # Instead of taking the last value, take the average of the plateau
            'findloglogwindow': False,
            'timeselection': 'lnt',
            'mplstyle': mplstyle,
            'legendcols': [],  # Choose 'num', 'bmax','tsim'
            
            
        },
        'titr': {
            'default' : default,
            'groupname': 'measurements',
            'colname': 'algorithm_time',
            'titlename': 'Iteration Time',
            'ylabel': 'Avg Iteration Time [min]',
            'ydiff': True,  # Plot the difference in y instead.
            'yscale': 'log',
            'plotprefix': 'titr',
            'plotdir': Path(plotdir, Path(mplstyle).stem),
            'normalize': 60,
            'realizations': 600,
            'findsaturation': False,  # Instead of taking the last value, take the average of the plateau
            'findloglogwindow': False,
            'timeselection': 'lnt',
            'mplstyle': mplstyle,
            'legendcols': [],  # Choose 'num', 'bmax','tsim'
        },
        'trn': {
            'default' : default,
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
            'plotdir': Path(plotdir, Path(mplstyle).stem),
            'findsaturation': False,  # Instead of taking the last value, take the average of the plateau
            'findloglogwindow': False,
            'markloglogwindow': False,
            'timeselection': 'lnt',
            'mplstyle': mplstyle,
            'legendcols': ['w', 'bmax:.0f', 'bavg:.0f', 'tsim'],  # Choose 'num', 'bmax','tsim'
            'legendlocation': 'best',
        },
        'lbit-avg': {
            'default' : default,
            'groupname': 'lbits',
            'dsetname': 'corrmat',
            'titlename': 'l-bit decay fit $C e^{-(|i-j|/\\xi)^\\beta}$',
            'filter': {
                'L': [12,20],
                'f': [0.2,0.5],
                # 'u': [16],
            },
            # 'titlename': 'l-bit decay fit $C e^{-|i-j|/\\xi}$',
            'ylabel': '$\log_{10} \langle \langle O_{ij}\\rangle\\rangle_\mathrm{arithmetic}$ ',
            # 'ylabel': '$\log_{10} \\bar O(|i-j|)$ ',
            'xlabel': "$j$",
            'xticks': [0, 6, 12, 18] if prb else None,
            # 'yticks': [1e0, 1e-4, 1e-9],
            'yticks': [0, -3, -6, -9, -12, -15],
            # 'yscale': 'log',
            # 'ynopos': 'mask',
            'plotprefix': 'lbit-arithmetic',
            'plotdir': Path(plotdir, Path(mplstyle).stem),
            'mplstyle': mplstyle,
            'ymin': -16,
            'xnormalize': True,
            # 'xmax': 32 if prb else None,
            # 'ymin': 1e-14,
            # 'legendcols': ['f', 'tstd', 'tgw8', 'cstd', 'cgw8', 'ubond'],  # Choose 'num', 'bmax','tsim'
            'legendcols': [] if prb else [],  # Choose 'num', 'bmax','tsim'
            'legendfits': ['xi', 'beta'] if prb else ['C', 'xi', 'beta', 'pos'],
            'legendoutside': True,
            'legendcollect': False,
            'legendlocation': (0.48, 0.58) if prb else 'center left',
            # 'legendtitle': '$y = C e^{-|i-j|/\\xi_\\tau}$',
            # 'legendtitle': '$\log \\bar O(x) = a - x \\xi_\\tau^{-1}$',
            'lbit-site': [0, 'mid', 'last'],
            'lbit-mean': 'arithmetic',
            'fit-beta': True,
            'fit-ymin': None,
            'fit-skip': 0 if prb else 0,
            'fit-mark': False,
            'fit-plot': True,
            'inset-cls': {
                # 'pos': [0.03, 0.6, 0.40, 0.40], # Positon of the inset, x0 y0 width height
                'pos': [0.17, 0.15, 0.25, 0.25],  # Positon of the inset, x0 y0 width height
                'coords': [None, None, None, None],
                # These zoom limits x1,x2,y1,y2, must be set by finding the maximum log log window
                'legendtitle': '$\\xi_\\tau$',
            },
        },
        'lbit-typ': {
            'default' : default,
            'groupname': 'lbits',
            'dsetname': 'corrmat',
            'titlename': 'l-bit decay fit $C e^{-(|i-j|/\\xi)^\\beta}$ (geometric avg)',
            'filter': {
                'L': [12,20],
                'f': [0.2,0.5],
                # 'u': [16],
            },
            # 'titlename': 'l-bit decay fit $C e^{-|i-j|/\\xi}$',
            # 'ylabel': '$\langle \langle O(|i-j|) \\rangle\\rangle$ ',
            'ylabel': '$\log_{10} \langle \langle O_{ij}\\rangle \\rangle_\mathrm{geometric}$ ',
            'xlabel': "$j$",
            'xticks': [0, 6, 12, 18] if prb else None,
            # 'yticks': [1e0, 1e-4, 1e-9],
            'yticks': [0, -3, -6, -9, -12, -15],
            # 'yscale': 'log',
            # 'ynopos': 'mask',
            'plotprefix': 'lbit-geometric',
            'plotdir': Path(plotdir, Path(mplstyle).stem),
            'mplstyle': mplstyle,
            'ymin': -16,
            'xnormalize' : True,
            # 'xmax': 32 if prb else None,
            # 'ymin': 1e-14,
            # 'legendcols': ['f', 'tstd', 'tgw8', 'cstd', 'cgw8', 'ubond'],  # Choose 'num', 'bmax','tsim'
            'legendcols': [] if prb else [],  # Choose 'num', 'bmax','tsim'
            'legendfits': ['xi', 'beta'] if prb else ['C', 'xi', 'beta', 'pos'],
            'legendlocation': (0.48, 0.58) if prb else 'center left',
            'legendoutside': True,
            'legendcollect': False,
            'lbit-site': [0, 'mid', 'last'],
            'lbit-mean': 'geometric',
            'fit-beta': True,
            'fit-ymin': None,
            'fit-skip': 0 if prb else 0,
            'fit-mark': False,
            'fit-plot': True,
            'inset-cls': {
                # 'pos': [0.03, 0.6, 0.40, 0.40], # Positon of the inset, x0 y0 width height
                'pos': [0.17, 0.15, 0.25, 0.25],  # Positon of the inset, x0 y0 width height
                'coords': [None, None, None, None],
                # These zoom limits x1,x2,y1,y2, must be set by finding the maximum log log window
                'legendtitle': '$\\xi_\\tau$',
            },
        },

        'cls-avg': {
            'default' : default,
            'groupname': 'lbits',
            'dsetname': 'corrmat',
            'titlename': '$\ell$-bit localization length (arithmetic avg.)',
            'filter': {
                # 'L': [8,12,16,20],
                # 'f': [0.2, 0.3, 0.4, 0.6, 0.8, 1.0],
                # 'u': [4, 8, 12, 16],
            },
            'ylabel': '$\langle\langle \\xi_\\tau \\rangle\\rangle$',
            # 'yticks': [0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4] if prb else None,
            # 'xticks': [16,32,48,64,80],
            'plotprefix': 'cls-avg',
            'plotdir': Path(plotdir, Path(mplstyle).stem),
            'mplstyle': mplstyle,
            'legendcols': [] if prb else ['f', 'x'],  # Choose 'num', 'bmax','tsim'
            'legendoutside': True,
            'legendcollect': True,
            'legendfits': ['xi', 'beta'] if prb else ['C', 'xi', 'beta', 'pos'],
            'lbit-site': ['mid'],
            'lbit-mean': 'arithmetic',
            'fit-beta': True,
            'fit-ymin': None,
            'fit-skip': 0 if prb else 0,
            # 'legendlocation': 'best',
            'legendlocation': 'best',  # "(0.52, 0.05),
        },
        'cls-typ': {
            'default' : default,
            'groupname': 'lbits',
            'dsetname': 'corrmat',
            'titlename': '$\ell$-bit localization length (geometric avg.)',
            'filter': {
                # 'L': [8,12,16,20],
                # 'f': [0.2, 0.3, 0.4, 0.6, 0.8, 1.0],
                #'u': [4, 8, 12, 16],
            },
            'ylabel': '$\langle\langle \\xi_\\tau \\rangle\\rangle$',
            # 'yticks': [0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4] if prb else None,
            # 'xticks': [16,32,48,64,80],
            'plotprefix': 'cls-typ',
            'plotdir': Path(plotdir, Path(mplstyle).stem),
            'mplstyle': mplstyle,
            'legendcols': [] if prb else ['f'],  # Choose 'num', 'bmax','tsim'
            'legendoutside': True,
            'legendcollect': True,
            'legendfits': ['xi', 'beta'] if prb else ['C', 'xi', 'beta', 'pos'],
            'lbit-site': ['mid'],
            'lbit-mean': 'geometric',
            'fit-beta': True,
            'fit-ymin': None,
            'fit-skip': 0 if prb else 0,
            # 'legendlocation': 'best',
            'legendlocation': 'best',  # "(0.52, 0.05),
        },
        'dist-num': {
            'default' : default,
            'groupname': 'measurements',
            'dsetname': 'number_entropy',
            'normpage': False,
            'titlename': 'Number Entropy',
            # 'figsize': (3.375, 3.00),
            'ylabel': '$p(S_N(L/2))$',
            'yscale': 'log',
            # 'yformat': '%.2f',
            'sharex': 'all',
            'sharey': 'all',
            'tidx': 'window',  # Time indices for which to plot the distribution
            'bins': 40,
            'xlabel': '$S_N$',
            'density': True,
            'plotprefix': 'SN',
            'plotdir': Path(plotdir, Path(mplstyle).stem),
            # 'ymin': 0.41,
            # 'ymax': 0.45,
            # 'xmin': 1,
            'findsaturation': True,  # Instead of taking the last value, take the average of the plateau
            'findloglogwindow': True,
            'markloglogwindow': True,
            'fitloglogwindow': True,
            # 'timeselection': 'lnt',
            'mplstyle': mplstyle,
            # 'legendcols': ['f', 'x', 'num', 'bmax:.0f', 'bavg:.0f', 'tsim', 't:.1e'],  # Choose 'num', 'bmax','tsim'
            # 'legendcols': ['L', 'f', 'w', 't:<8.1e'],  # Choose 'num', 'bmax','tsim'
            'legendcols': ['L', 'f'],  # Choose 'num', 'bmax','tsim'
            # 'legendcols': ['f', 'x', 'w', 'num', 't:<8.1e'],  # Choose 'num', 'bmax','tsim'
            'legendoutside': False,
            'legendcollect': False,
            'legendlocation': 'lower left',
        },
        'divg-ent': {  # Distribution of infinite time averaged entropy
            'default' : default,
            'groupname': 'measurements',
            'dsetname': 'entanglement_entropy',
            'normpage': False,
            'titlename': 'Entanglement Entropy',
            # 'figsize': (3.375, 3.00),
            'ylabel': '$p(S_E^\infty(L/2))$',
            'yscale': 'log',
            # 'yformat': '%.2f',
            'sharex': 'all',
            'sharey': 'all',
            # 'xmin': 0.0,
            # 'xmax': 1.5,
            'tidx': 'window',  # Time indices for which to plot the distribution
            'bins': 60,
            'xlabel': '$S_E^\infty$',
            'density': True,
            'plotprefix': 'SE',
            'plotdir': Path(plotdir, Path(mplstyle).stem),
            'ymin': 1e-3,
            # 'ymax': 0.45,
            # 'xmin': 1,
            'mplstyle': mplstyle,
            # 'legendcols': ['f', 'x', 'num', 'bmax:.0f', 'bavg:.0f', 'tsim', 't:.1e'],  # Choose 'num', 'bmax','tsim'
            # 'legendcols': ['f', 'w'],  # Choose 'num', 'bmax','tsim'
            'legendcols': ['L', 'f'],  # Choose 'num', 'bmax','tsim'
            'legendoutside': False,
            'legendcollect': False,
            'legendlocation': 'lower left',
        },
        'divg-num': {  # Distribution of infinite time averaged entropy
            'default' : default,
            'groupname': 'measurements',
            'dsetname': 'number_entropy',
            'normpage': False,
            'titlename': 'Number Entropy',
            'filter': {
                #'L': [12, 16, 20],
                #'f': [0.2,0.3,0.4,0.5],
                #'u': [16],
            },
            # 'figsize': (3.375, 3.00),
            'ylabel': '$p(S_N^\infty(L/2))$',
            'yscale': 'log',
            # 'yformat': '%.2f',
            'sharex': 'all',
            'sharey': 'all',
            'xmin': 0.0,
            'xmax': 1.2,
            'ymin': 1e-3,
            'ymax': 1e+1,
            'tidx': 'window',  # Time indices for which to plot the distribution
            'bins': 60,
            'xlabel': '$S_N^\infty$',
            'density': True,
            'marklog2':True,
            'plotprefix': 'SN',
            'plotdir': Path(plotdir, Path(mplstyle).stem),
            # 'ymin': 0.41,
            # 'ymax': 0.45,
            # 'xmin': 1,
            'mplstyle': mplstyle,
            # 'legendcols': ['f', 'x', 'num', 'bmax:.0f', 'bavg:.0f', 'tsim', 't:.1e'],  # Choose 'num', 'bmax','tsim'
            # 'legendcols': ['f', 'w'],  # Choose 'num', 'bmax','tsim'
            'legendcols': ['L', 'f'],  # Choose 'num', 'bmax','tsim'
            # 'legendcols': ['L'],  # Choose 'num', 'bmax','tsim'
            'legendoutside': False,
            'legendcollect': False,
            'legendlocation': 'lower left',
        },
        'dist-mem': {  # Distribution of infinite time averaged entropy
            'default' : default,
            'groupname': 'mem_usage',
            'dsetname': 'data',
            'colname' : 'hwm',
            'titlename': 'Max Memory Usage',
            'filter': {
                # 'L': [12, 16, 20],
                # 'f': [0.2,0.3,0.4,0.5],
                # 'u': [16],
            },
            # 'figsize': (3.375, 3.00),
            'ylabel': 'Histogram',
            'yscale': 'log',
            # 'yformat': '%.2f',
            'sharex': 'all',
            'sharey': 'all',
            # 'xmin': 0.0,
            # 'xmax': 1.2,
            # 'ymin': 1e-5,
            # 'tidx': 'window',  # Time indices for which to plot the distribution
            'bins': 60,
            'xlabel': 'RAM [MB]',
            'density': False,
            'plotprefix': 'mem',
            'plotdir': Path(plotdir, Path(mplstyle).stem),
            # 'ymin': 0.41,
            # 'ymax': 0.45,
            # 'xmin': 1,
            'findsaturation': True,  # Instead of taking the last value, take the average of the plateau
            'findloglogwindow': True,
            'markloglogwindow': True,
            'fitloglogwindow': True,
            # 'timeselection': 'lnt',
            'mplstyle': mplstyle,
            # 'legendcols': ['f', 'x', 'num', 'bmax:.0f', 'bavg:.0f', 'tsim', 't:.1e'],  # Choose 'num', 'bmax','tsim'
            # 'legendcols': ['f', 'w'],  # Choose 'num', 'bmax','tsim'
            'legendcols': ['L', 'f', 'bavg:.0f', 'tsim'],  # Choose 'num', 'bmax','tsim'
            # 'legendcols': ['L'],  # Choose 'num', 'bmax','tsim'
            'legendoutside': False,
            'legendcollect': False,
            'legendlocation': 'upper right',
        },
        'tavg-ent': {
            'default' : default,
            'groupname': 'measurements',
            'dsetname': 'entanglement_entropy',
            'normpage': False,
            'titlename': 'Average Entanglement entropy as $t\\rightarrow \infty$',
            # 'figsize': (3.375, 3.00),
            'ylabel': '$\langle \langle S_\mathrm{E}^\infty(L/2) \\rangle \\rangle$',
            'yformat': '%.2f',
            'xmaloc': MaxNLocator(integer=True),
            'sharex': 'all',
            'sharey': 'all',
            'xlabel': '$L$',
            'density': True,
            'plotprefix': 'SE',
            'plotdir': Path(plotdir, Path(mplstyle).stem),
            # 'ymin': 0.41,
            # 'ymax': 0.45,
            # 'xmin': 1,
            'findsaturation': True,  # Instead of taking the last value, take the average of the plateau
            'findloglogwindow': True,
            'markloglogwindow': True,
            'fitloglogwindow': True,
            # 'timeselection': 'lnt',
            'mplstyle': mplstyle,
            # 'legendcols': ['f', 'x', 'num', 'bmax:.0f', 'bavg:.0f', 'tsim', 't:.1e'],  # Choose 'num', 'bmax','tsim'
            'legendcols': ['L', 'f'],  # Choose 'num', 'bmax','tsim'
            'legendoutside': False,
            'legendcollect': False,
            'legendlocation': 'lower left',
        },
        'tavg-num': {
            'default' : default,
            'groupname': 'measurements',
            'dsetname': 'number_entropy',
            'normpage': False,
            'titlename': 'Average Number Entropy as $t\\rightarrow \infty$',
            # 'figsize': (3.375, 3.00),
            'ylabel': '$\langle \langle S_\mathrm{N}^\infty(L/2) \\rangle \\rangle$',
            'yformat': '%.2f',
            'xmaloc': MaxNLocator(integer=True, nbins=4),
            'sharex': 'all',
            'sharey': 'all',
            'xlabel': '$L$',
            'density': True,
            'plotprefix': 'SN',
            'plotdir': Path(plotdir, Path(mplstyle).stem),
            # 'ymin': 0.41,
            # 'ymax': 0.45,
            # 'xmin': 1,
            'findsaturation': True,  # Instead of taking the last value, take the average of the plateau
            'findloglogwindow': True,
            'markloglogwindow': True,
            'fitloglogwindow': True,
            # 'timeselection': 'lnt',
            'mplstyle': mplstyle,
            # 'legendcols': ['f', 'x', 'num', 'bmax:.0f', 'bavg:.0f', 'tsim', 't:.1e'],  # Choose 'num', 'bmax','tsim'
            # 'legendcols': ['L', 'u', 'f', 'x', 'num'],  # Choose 'num', 'bmax','tsim'
            # 'legendcols': ['f', 'x', 'num'],  # Choose 'num', 'bmax','tsim'
            'legendcols': ['f'],  # Choose 'num', 'bmax','tsim'
            'legendoutside': False,
            'legendcollect': False,
            'legendlocation': 'lower right',
        },
        'dist-tsim': {
            'default' : default,
            'groupname': 'measurements',
            'dsetname': 'algorithm_time',
            'normpage': False,
            'titlename': 'Simulation time',
            # 'figsize': (3.375, 3.00),
            'ylabel': '$p(t_\mathrm{sim})$',
            'yscale': 'log',
            'xscale': 'log',
            # 'yformat': '%.2f',
            'sharex': 'all',
            'sharey': 'all',
            'tidx': [-1],  # Time indices for which to plot the distribution
            'bins': 60,
            'normalize': 60,
            'xlabel': '$t_\mathrm{sim}$',
            'density': True,
            'plotprefix': 'tsim',
            'plotdir': Path(plotdir, Path(mplstyle).stem),
            'ymin': 1e-4,
            'xmin': 5e-2,
            'xmax': 1e+3,
            'findsaturation': False,  # Instead of taking the last value, take the average of the plateau
            'findloglogwindow': False,
            'markloglogwindow': False,
            'fitloglogwindow': False,
            # 'timeselection': 'lnt',
            'mplstyle': mplstyle,
            # 'legendcols': ['f', 'x', 'num', 'bmax:.0f', 'bavg:.0f', 'tsim', 't:.1e'],  # Choose 'num', 'bmax','tsim'
            # 'legendcols': ['L', 'f', 'w', 't:<8.1e'],  # Choose 'num', 'bmax','tsim'
            'legendcols': ['L', 'tsim'],  # Choose 'num', 'bmax','tsim'
            # 'legendcols': ['f', 'x', 'w', 'num', 't:<8.1e'],  # Choose 'num', 'bmax','tsim'
            'legendoutside': False,
            'legendcollect': False,
            'legendlocation': 'upper right',
        },
        'dist-chi': {
            'default' : default,
            'groupname': 'measurements',
            'dsetname': 'bond_mid',
            'normpage': False,
            'titlename': 'Midchain Bond dimension ($t\\rightarrow\infty)$',
            # 'figsize': (3.375, 3.00),
            'ylabel': '$p(\chi)$',
            # 'yscale': 'log',
            'xscale': 'log',
            # 'yformat': '%.2f',
            'sharex': 'all',
            'sharey': 'all',
            'tidx': [-1],  # Time indices for which to plot the distribution
            'bins': 20,
            'xlabel': '$\chi$',
            'density': True,
            'plotprefix': 'bond',
            'plotdir': Path(plotdir, Path(mplstyle).stem),
            # 'ymin': 0.41,
            # 'ymax': 0.45,
            # 'xmin': 1,
            'findsaturation': False,  # Instead of taking the last value, take the average of the plateau
            'findloglogwindow': False,
            'markloglogwindow': False,
            'fitloglogwindow': False,
            # 'timeselection': 'lnt',
            'mplstyle': mplstyle,
            # 'legendcols': ['f', 'x', 'num', 'bmax:.0f', 'bavg:.0f', 'tsim', 't:.1e'],  # Choose 'num', 'bmax','tsim'
            # 'legendcols': ['L', 'f', 'w', 't:<8.1e'],  # Choose 'num', 'bmax','tsim'
            'legendcols': ['L', 'tsim'],  # Choose 'num', 'bmax','tsim'
            # 'legendcols': ['f', 'x', 'w', 'num', 't:<8.1e'],  # Choose 'num', 'bmax','tsim'
            'legendoutside': False,
            'legendcollect': False,
            'legendlocation': 'upper right',
        },
        'linearFit-num2': {
            'default' : default,
            'groupname': 'measurements',
            'colname': 'number_entropy',
            'normpage': False,
            'titlename': 'Number Entropy fit to $\\alpha + \\beta \log \log t$',
            'ylabel': '$\\beta$',
            'yformat': '%.3f',
            'plotprefix': 'SN-linearFit',
            'plotdir': Path(plotdir, Path(mplstyle).stem),
            'sharex': 'all',
            'sharey': 'none',
            # 'xmin': 0.0,
            # 'xmax': 2.1,
            'findsaturation': True,  # Instead of taking the last value, take the average of the plateau
            'findloglogwindow': True,
            'markloglogwindow': True,
            'fitloglogwindow': True,
            'fillerror': False,
            'timeselection': 'lnlnt',
            'mplstyle': mplstyle,
            # 'legendcols': ['f', 'x', 'num', 'bmax:.0f', 'bavg:.0f', 'tsim'],  # Choose 'num', 'bmax','tsim'
            'legendcols': ['L', 'f', 'tstd', 'tgw8', 'cstd', 'cgw8'],  # Choose 'num', 'bmax','tsim'
            
            
            'legendlocation': 'center left',
            # 'legendlocation': (0.01, 0.65),
        },
        'rise-num2': {
            'default' : default,
            'groupname': 'measurements',
            'colname': 'number_entropy',
            'normpage': False,
            'titlename': 'Number Entropy rise',
            'ylabel': '$S_N(t_2) - S_N(t_1)$',
            'yformat': '%.3f',
            'plotprefix': 'SN-rise',
            'plotdir': Path(plotdir, Path(mplstyle).stem),
            'sharex': 'all',
            'sharey': 'none',
            # 'xmin': 0.0,
            # 'xmax': 2.75,
            'relative': False,
            'ylabel_relative': '$\\frac{S_N(t_2) - S_N(t_1)}{S_N(t_2)}$',
            'findsaturation': True,  # Instead of taking the last value, take the average of the plateau
            'findloglogwindow': True,
            'markloglogwindow': True,
            'fitloglogwindow': True,
            'fillerror': False,
            'timeselection': 'lnlnt',
            'mplstyle': mplstyle,
            # 'legendcols': ['f', 'x', 'num', 'bmax:.0f', 'bavg:.0f', 'tsim'],  # Choose 'num', 'bmax','tsim'
            'legendcols': ['L', 'f'],  # Choose 'num', 'bmax','tsim'
            
            
            'legendlocation': 'center right',
            # 'legendlocation': (0.01, 0.65),
        },
        'num-svnt': {
            'default' : default,
            'groupname': 'measurements',
            'colname': 'number_entropy',
            'normpage': True,
            'titlename': 'Number Entropy',
            'filter': {
                # 'L': [8,12,16,20],
                # 'f': [0.3,0.5],
                # 'u': [16],
            },
            'ylabel': '$\langle \langle S_N(L/2)\\rangle \\rangle$',
            'use_configurational_entropy': False,
            'yformat': '%.3f',
            'plotprefix': 'SN',
            'plotdir': Path(plotdir, Path(mplstyle).stem),
            'ymin': 0.1250,
            'ymax': 0.3150,
            'xmin': 0.02,
            'xmax': 0.16,
            'findsaturation': True,  # Instead of taking the last value, take the average of the plateau
            'findloglogwindow': True,
            'markloglogwindow': True,
            'fitloglogwindow': True,
            'shadederror': False,
            'timeselection': 't',,
            'mplstyle': mplstyle,
            # 'legendcols': ['f', 'x', 'num', 'bmax:.0f', 'bavg:.0f', 'tsim'],  # Choose 'num', 'bmax','tsim'
            'legendcols': ['L', 'f', 'bavg:.0f','tsim'],  # Choose 'num', 'bmax','tsim'
            
            
            'legendlocation': 'lower right',
            # 'legendlocation': (0.01, 0.65),
        },
        'crossup': {
            'default' : default,
            'groupname': 'lbits',
            'dsetname': 'corrmat',
            'titlename': '$\ell$-bit half-chain crossover',
            'filter': {
                # 'L': [8,12,16,20],
                # 'f': [0.2,0.3,0.4,0.5,1.0],
                # 'u': [16],
            },
            # 'titlename': 'l-bit decay fit $C e^{-|i-j|/\\xi}$',
            'ylabel': '$X_\ell$',
            # 'ylabel': '$\log_{10} \\bar O(|i-j|)$ ',
            # 'xlabel': "$j$",
            # 'xticks': [0, 6, 12, 18] if prb else None,
            # 'yticks': [1e0, 1e-4, 1e-9],
            # 'yticks': [0, -3, -6, -9, -12, -15],
            # 'yscale': 'log',
            # 'ynopos': 'mask',
            'plotprefix': 'lbit',
            'plotdir': Path(plotdir, Path(mplstyle).stem),
            'mplstyle': mplstyle,
            # 'ymin': -16,
            # 'xmax': 32 if prb else None,
            # 'ymin': 1e-14,
            # 'legendcols': ['f', 'tstd', 'tgw8', 'cstd', 'cgw8', 'ubond'],  # Choose 'num', 'bmax','tsim'
            'legendcols': [] if prb else [],  # Choose 'num', 'bmax','tsim'
            'legendfits': ['xi', 'beta'] if prb else ['C', 'xi', 'beta', 'pos'],
            
            
            'legendlocation': 'lower right',
            # 'legendtitle': '$y = C e^{-|i-j|/\\xi_\\tau}$',
            # 'legendtitle': '$\log \\bar O(x) = a - x \\xi_\\tau^{-1}$',
            'lbit-mean': 'arithmetic',
            'lbit-site': [0, 'mid', 'last'],

            'fit-beta': True,
            'fit-ymin': None,
            'fit-skip': 0 if prb else 0,
            'inset-cls': {
                # 'pos': [0.03, 0.6, 0.40, 0.40], # Positon of the inset, x0 y0 width height
                'pos': [0.17, 0.15, 0.25, 0.25],  # Positon of the inset, x0 y0 width height
                'coords': [None, None, None, None],
                # These zoom limits x1,x2,y1,y2, must be set by finding the maximum log log window
                'legendtitle': '$\\xi_\\tau$',
            },
        }

    }



    return meta
