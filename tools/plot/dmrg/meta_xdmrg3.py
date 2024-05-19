from pathlib import Path
import numpy as np
from matplotlib.ticker import LogLocator, \
    LogFormatterMathtext, FixedLocator, ScalarFormatter

mplstyle = 'src/plots/stylesheets/prl.mplstyle'
# mplstyle = '../common/stylesheets/slack.mplstyle'
prl = 'prl' in mplstyle


def get_meta(plotdir, cachedir):
    # figsize1x1_halfcol = 3.404 * 0.55, 3.404 * 0.50, # Half-size in one column of a two-column document

    # figsize1x1_halfcol = 3.404 * 0.55, 3.404 * 0.50,  # Half-size in one column of a two-column document
    # subplots1x1 = {
    #     'top': 1.0,
    #     'bottom': 0.115,
    #     'left': 0.25,
    #     'right': 0.99,
    # }
    # default = {
    #     'box_aspect': 1,
    #     'subspec_title': False,
    #     'figspec_title': False,
    #     'legendoutside': False,
    #     'legendcollect': False,
    #     'constrained_layout': False,
    #     'axtitle': False,
    #     'figsize': figsize1x1_halfcol,
    #     'subplots': subplots1x1,
    # }

    # In revtex4-2 with a two-column article we get (in inches)
    # \usepackage{xprintlen}
    # \printlen[10][in]{\columnwidth} : 3.4039025050 in = 86.45958 mm
    # \printlen[10][in]{\linewidth}   : 3.4039025050 in = 86.45958 mm
    # \printlen[10][in]{\textwidth}   : 7.0568710470 in = 179.24548 mm

    # For a 2x3 grid, the easiest is to divide the width 7.05687/3=2.35229 for each figure
    # 2.26926: is 66% of a column width
    # All the sizes below refer to the size of ONE single plot inside a nxn grid in either spanning 1 or 2 columns
    #
    # figsize1x1_halfcol = 2.35229, 3.404 * 0.450, # Third-size for spanning a two-column document  with 3 figures

    box_aspect = 2. / 3
    figsize1x1_1col = 3.4039, 3.4039 * box_aspect,  # Third-size for spanning a two-column document  with 3 figures
    figsize1x2_1col = 0.5 * 3.4039, 0.7 * 3.4039,  # For two plots side by side in 1 column of a two-column doc
    figsize2x3_2col = 7.0568 / 3, 7.0568 / 3 * box_aspect,  # Third-size for spanning a two-column document  with 3 figures
    figsize2x2_2col = 7.0568 / 2, 7.0568 / 2 * box_aspect,  # Third-size for spanning a two-column document  with 3 figures
    figsize3x3_2col = 7.0568 / 1, 7.0568 / 1 * box_aspect,  # Third-size for spanning a two-column document  with 3 figures

    # figsize1x1_halfcol = 2.26926, 3.4039 * 0.440, # Third-size for spanning a two-column document  with 3 figures
    # figsize1x1_halfcol = 2.26926, 3.4039 * 0.440, # Third-size for spanning a two-column document  with 3 figures

    subplots1x1_1col = {
        'top': 0.997,
        'bottom': 0.11,
        'left': 0.110,
        'right': 0.997,
        'wspace': 0,
        'hspace': 0,
    }
    subplots1x2_1col = {
        'top': 0.997,
        'bottom': 0.10,
        'left': 0.21,
        'right': 0.997,
        'wspace': 0,
        'hspace': 0,
    }
    subplots2x3_2col = {
        'top': 0.997,
        'bottom': 0.17,
        'left': 0.180,
        'right': 0.997,
        'wspace': 0,
        'hspace': 0,
    }
    subplots2x2_2col = {
        'top': 0.997,
        'bottom': 0.12,
        'left': 0.100,
        'right': 0.997,
        'wspace': 0,
        'hspace': 0,
    }
    subplots3x3_2col = {
        'top': 0.997,
        'bottom': 0.05,
        'left': 0.050,
        'right': 0.997,
        'wspace': 0,
        'hspace': 0,
    }
    default = {
        'box_aspect': box_aspect,
        'subspec_title': False,
        'figspec_title': False,
        'legendoutside': False,
        'legendcollect': False,
        'constrained_layout': False,
        'axtitle': False,
        'figsize': figsize2x3_2col,
        'subplots': subplots2x3_2col,
        'cachedir': Path(cachedir),
        'plotdir': Path(plotdir, Path(mplstyle).stem),
    }
    default2x2 = {
        'box_aspect': box_aspect,
        'subspec_title': False,
        'figspec_title': False,
        'legendoutside': False,
        'legendcollect': False,
        'constrained_layout': False,
        'axtitle': False,
        'figsize': figsize2x2_2col,
        'subplots': subplots2x2_2col,
        'cachedir': Path(cachedir),
        'plotdir': Path(plotdir, Path(mplstyle).stem),
    }
    default3x3 = {
        'box_aspect': box_aspect,
        'subspec_title': False,
        'figspec_title': False,
        'legendoutside': False,
        'legendcollect': False,
        'constrained_layout': False,
        'axtitle': False,
        'figsize': figsize3x3_2col,
        'subplots': subplots3x3_2col,
        'cachedir': Path(cachedir),
        'plotdir': Path(plotdir, Path(mplstyle).stem),
    }
    default1x1_1col = {
        'box_aspect': box_aspect,
        'subspec_title': False,
        'figspec_title': False,
        'legendoutside': False,
        'legendcollect': False,
        'constrained_layout': False,
        'axtitle': False,
        'figsize': figsize1x1_1col,
        'subplots': subplots1x1_1col,
        'cachedir': Path(cachedir),
        'plotdir': Path(plotdir, Path(mplstyle).stem),
    }
    default1x2_1col = {
        'box_aspect': 1.6175,
        'subspec_title': False,
        'figspec_title': False,
        'legendoutside': False,
        'legendcollect': False,
        'constrained_layout': False,
        'axtitle': False,
        'figsize': figsize1x2_1col,
        'subplots': subplots1x2_1col,
        'cachedir': Path(cachedir),
        'plotdir': Path(plotdir, Path(mplstyle).stem),
    }
    subplots_halfheight = {
        'top': 1.0,
        'bottom': 0.15,
        'left': 0.02,  # 0.17,#0.02,
        'right': 0.72,  # 1.0,
        'wspace': 0,
        'hspace': 0,
    }
    default_halfheight = {
        'box_aspect': 0.95,  # 0.86 0#0.46,
        'figsize': (3.404 * 5.0 / 9, 3.404 * 0.25),
        'subspec_title': False,
        'figspec_title': False,
        'legendoutside': False,
        'legendcollect': False,
        'subplots': subplots_halfheight,
        'cachedir': Path(cachedir),
        'plotdir': Path(plotdir, Path(mplstyle).stem),
    }
    figsize1x1_inset = 2.26926 * 0.34, 3.404 * 0.440 * 0.42,  # Tiny size for insets
    subplots_inset = {
        'top': 0.94,
        'bottom': 0.18,
        'left': 0.36,  # 0.17,#0.02,
        'right': 0.970,  # 1.0,
        'wspace': 0,
        'hspace': 0,
    }
    default_inset = {
        'box_aspect': 1.0,  # 0.46,
        'figsize': figsize1x1_inset,
        'subspec_title': False,
        'figspec_title': False,
        'legendoutside': False,
        'legendcollect': False,
        'subplots': subplots_inset,
        'font.size': 7.0,
    }

    meta = {
        'common': {
            'cachedir': Path(cachedir),
            'plotdir': Path(plotdir, Path(mplstyle).stem),
        },

        'opdm-spectrum-g': {
            'default': default,
            'include': {
                'L': [12,14],
                'g': [0.0, 0.01, 0.02, 0.03],
                'd': [0.0],
            },
            'groupbase': 'tables',
            'groupname': 'opdm_spectrum',
            'colname': 'eigval',
            'figspec': ['L'],
            'subspec': ['d'],
            'linspec': ['g'],
            'plotdir': Path(plotdir, Path(mplstyle).stem),
            'filename': 'opdm-spectrum-d',
            'palettes': ["viridis_r"],
            'mplstyle': mplstyle,
            # 'titlename': 'l-bit decay fit $C e^{-(|i-j|/\\xi)^\\beta}$',
            # 'filter': {
            # 'L': [24],
            # 'f': [0.2,0.4],
            # 'u': [16],
            # },
            # 'titlename': 'l-bit decay fit $C e^{-|i-j|/\\xi}$',
            # 'ylabel': '$\log_{10} \langle \langle O(|i-j|) \\rangle\\rangle$ ',
            'ylabel': 'opdm spectrum',
            # 'ylabel': '$\log_{10} \\bar O(|i-j|)$ ',
            'xlabel': "eigv $i$",
            # 'xlabelpad': -8,
            # 'xcoords': (0.5, -0.04),
            # 'xticks': [0, 0.25, 0.5, 0.75, 1.0] if prl else None,
            # 'xticks': [0, 15],
            # 'xticklabels': ["0","27"],
            # 'yticks': [1e-2, 1e-12],
            # 'ylabelpad': -16,

            # 'yticklabels': ['$-1$','$-15$'],
            # 'ycoords': (-0.34, 0.34),
            # 'xticklabels': ['$0$', '$L/2$', '$L$'],
            # 'yticks': [0, -3, -6, -9, -12, -15],
            # 'yscale': 'log',
            # 'ynopos': 'mask',

            # 'ymin': -16,
            # 'ymin': 0,
            # 'ymax': 2,
            # 'xmin': 0,
            # 'xmax': 15,
            # 'xnormalize': False,
            # 'xshift2mid': False,
            # 'xmax': 32 if prl else None,
            # 'ymin': 1e-14,
            # 'legendcols': ['f', 'tstd', 'tgw8', 'cstd', 'cgw8', 'ubond'],  # Choose 'num', 'bmax','tsim'
            'legendcols': ['L', 'g', 'd'],  # Choose 'num', 'bmax','tsim'
            # 'legendfits': ['xi', 'beta'] if prl else ['C', 'xi', 'beta', 'pos'],
            'legendoutside': False,
            'legendcollect': False,
            # 'legendlocation': (0.18, 0.0),
            'legendlocation': ['upper right', 'upper right', 'upper right', 'upper right', ],
            # 'bbox_to_anchor': (1.0, 0.70),  # Use with loc 'upper right'
            'bbox_to_anchor': (1.00, 1.01),
            'frameon': False,
            # 'legendtitle': 'Arithmetic average',
            'legendtitle': None,  # '$\\overline O \propto e^{(\\frac{|i-j|}{\\xi_\\tau})^\\beta}$',
            # 'legendtitle': '$y = C e^{-|i-j|/\\xi_\\tau}$',
            # 'legendtitle': '$\log \\bar O(x) = a - x \\xi_\\tau^{-1}$',
            # 'lbit-site': [0, 'mid', 'last'],
            # 'lbit-site': ['mid'],
            # 'lbit-mean': 'arithmetic',
            # 'lbit-axis': '',
            # 'fit-beta': True,
            # 'fit-ymin': 1e-16,
            # 'fit-skip': 0 if prl else 0,
            # 'fit-mark': False,
            # 'fit-plot': False,
            # 'inset-cls': {
            #     # 'pos': [0.03, 0.6, 0.40, 0.40], # Positon of the inset, x0 y0 width height
            #     'pos': [0.17, 0.15, 0.25, 0.25],  # Positon of the inset, x0 y0 width height
            #     'coords': [None, None, None, None],
            #     # These zoom limits x1,x2,y1,y2, must be set by finding the maximum log log window
            #     'legendtitle': '$\\xi_\\tau$',
            # },
        },
        'opdm-spectrum-d': {
            'default': default,
            'include': {
                'L': [12,14],
                'g': [0.03],
                'd': range(-4, 5, 2),
            },
            'groupbase': 'tables',
            'groupname': 'opdm_spectrum',
            'colname': 'eigval',
            'figspec': ['L'],
            'subspec': ['g'],
            'linspec': ['d'],
            'plotdir': Path(plotdir, Path(mplstyle).stem),
            'filename': 'opdm-spectrum-d',
            'palettes': ["viridis_r"],
            'mplstyle': mplstyle,
            # 'titlename': 'l-bit decay fit $C e^{-(|i-j|/\\xi)^\\beta}$',
            # 'filter': {
            # 'L': [24],
            # 'f': [0.2,0.4],
            # 'u': [16],
            # },
            # 'titlename': 'l-bit decay fit $C e^{-|i-j|/\\xi}$',
            # 'ylabel': '$\log_{10} \langle \langle O(|i-j|) \\rangle\\rangle$ ',
            'ylabel': 'opdm spectrum',
            # 'ylabel': '$\log_{10} \\bar O(|i-j|)$ ',
            'xlabel': "eigv $i$",
            # 'xlabelpad': -8,
            # 'xcoords': (0.5, -0.04),
            # 'xticks': [0, 0.25, 0.5, 0.75, 1.0] if prl else None,
            # 'xticks': [0, 15],
            # 'xticklabels': ["0","27"],
            # 'yticks': [1e-2, 1e-12],
            # 'ylabelpad': -16,

            # 'yticklabels': ['$-1$','$-15$'],
            # 'ycoords': (-0.34, 0.34),
            # 'xticklabels': ['$0$', '$L/2$', '$L$'],
            # 'yticks': [0, -3, -6, -9, -12, -15],
            # 'yscale': 'log',
            # 'ynopos': 'mask',

            # 'ymin': -16,
            # 'ymin': 0,
            # 'ymax': 2,
            # 'xmin': 0,
            # 'xmax': 15,
            # 'xnormalize': False,
            # 'xshift2mid': False,
            # 'xmax': 32 if prl else None,
            # 'ymin': 1e-14,
            # 'legendcols': ['f', 'tstd', 'tgw8', 'cstd', 'cgw8', 'ubond'],  # Choose 'num', 'bmax','tsim'
            'legendcols': ['L', 'g', 'd'],  # Choose 'num', 'bmax','tsim'
            # 'legendfits': ['xi', 'beta'] if prl else ['C', 'xi', 'beta', 'pos'],
            'legendoutside': False,
            'legendcollect': False,
            # 'legendlocation': (0.18, 0.0),
            'legendlocation': ['upper right', 'upper right', 'upper right', 'upper right', ],
            # 'bbox_to_anchor': (1.0, 0.70),  # Use with loc 'upper right'
            'bbox_to_anchor': (1.00, 1.01),
            'frameon': False,
            # 'legendtitle': 'Arithmetic average',
            'legendtitle': None,  # '$\\overline O \propto e^{(\\frac{|i-j|}{\\xi_\\tau})^\\beta}$',
            # 'legendtitle': '$y = C e^{-|i-j|/\\xi_\\tau}$',
            # 'legendtitle': '$\log \\bar O(x) = a - x \\xi_\\tau^{-1}$',
            # 'lbit-site': [0, 'mid', 'last'],
            # 'lbit-site': ['mid'],
            # 'lbit-mean': 'arithmetic',
            # 'lbit-axis': '',
            # 'fit-beta': True,
            # 'fit-ymin': 1e-16,
            # 'fit-skip': 0 if prl else 0,
            # 'fit-mark': False,
            # 'fit-plot': False,
            # 'inset-cls': {
            #     # 'pos': [0.03, 0.6, 0.40, 0.40], # Positon of the inset, x0 y0 width height
            #     'pos': [0.17, 0.15, 0.25, 0.25],  # Positon of the inset, x0 y0 width height
            #     'coords': [None, None, None, None],
            #     # These zoom limits x1,x2,y1,y2, must be set by finding the maximum log log window
            #     'legendtitle': '$\\xi_\\tau$',
            # },
        },
        'opdm-gapsize-g': {
            'default': default,
            'include': {
                'L': [12, 14],
                'g': [0.000, 0.005, 0.010, 0.015, 0.020, 0.025, 0.030],
                'd': [0.0],
            },
            'groupbase': 'tables',
            'groupname': 'opdm_spectrum',
            'colname': 'eigval',
            'figspec': ['d'],
            'subspec': ['d'],
            'linspec': ['L'],
            'xaxspec': ['g'],
            'plotdir': Path(plotdir, Path(mplstyle).stem),
            'filename': 'opdm-gapsize-g',
            'palettes': ["viridis_r"],
            'mplstyle': mplstyle,
            # 'titlename': 'l-bit decay fit $C e^{-(|i-j|/\\xi)^\\beta}$',
            # 'filter': {
            # 'L': [24],
            # 'f': [0.2,0.4],
            # 'u': [16],
            # },
            # 'titlename': 'l-bit decay fit $C e^{-|i-j|/\\xi}$',
            # 'ylabel': '$\log_{10} \langle \langle O(|i-j|) \\rangle\\rangle$ ',
            'ylabel': '$1-$ opdm gap',
            # 'ylabel': '$\log_{10} \\bar O(|i-j|)$ ',
            'xlabel': "$g$",
            # 'xlabelpad': -8,
            # 'xcoords': (0.5, -0.04),
            # 'xticks': [0, 0.25, 0.5, 0.75, 1.0] if prl else None,
            # 'xticks': [0, 15],
            # 'xticklabels': ["0","27"],
            # 'yticks': [1e-2, 1e-12],
            # 'ylabelpad': -16,

            # 'yticklabels': ['$-1$','$-15$'],
            # 'ycoords': (-0.34, 0.34),
            # 'xticklabels': ['$0$', '$L/2$', '$L$'],
            # 'yticks': [0, -3, -6, -9, -12, -15],
            # 'yscale': 'log',
            # 'ynopos': 'mask',

            # 'ymin': -16,
            # 'ymin': 0,
            # 'ymax': 2,
            # 'xmin': 0,
            # 'xmax': 15,
            # 'xnormalize': False,
            # 'xshift2mid': False,
            # 'xmax': 32 if prl else None,
            # 'ymin': 1e-14,
            # 'legendcols': ['f', 'tstd', 'tgw8', 'cstd', 'cgw8', 'ubond'],  # Choose 'num', 'bmax','tsim'
            'legendcols': ['L', 'g', 'd'],  # Choose 'num', 'bmax','tsim'
            # 'legendfits': ['xi', 'beta'] if prl else ['C', 'xi', 'beta', 'pos'],
            'legendoutside': False,
            'legendcollect': False,
            # 'legendlocation': (0.18, 0.0),
            'legendlocation': ['upper right', 'upper right', 'upper right', 'upper right', ],
            # 'bbox_to_anchor': (1.0, 0.70),  # Use with loc 'upper right'
            'bbox_to_anchor': (1.00, 0.7),
            'frameon': False,
            # 'legendtitle': 'Arithmetic average',
            'legendtitle': None,  # '$\\overline O \propto e^{(\\frac{|i-j|}{\\xi_\\tau})^\\beta}$',
            # 'legendtitle': '$y = C e^{-|i-j|/\\xi_\\tau}$',
            # 'legendtitle': '$\log \\bar O(x) = a - x \\xi_\\tau^{-1}$',
            # 'lbit-site': [0, 'mid', 'last'],
            # 'lbit-site': ['mid'],
            # 'lbit-mean': 'arithmetic',
            # 'lbit-axis': '',
            # 'fit-beta': True,
            # 'fit-ymin': 1e-16,
            # 'fit-skip': 0 if prl else 0,
            # 'fit-mark': False,
            # 'fit-plot': False,
            # 'inset-cls': {
            #     # 'pos': [0.03, 0.6, 0.40, 0.40], # Positon of the inset, x0 y0 width height
            #     'pos': [0.17, 0.15, 0.25, 0.25],  # Positon of the inset, x0 y0 width height
            #     'coords': [None, None, None, None],
            #     # These zoom limits x1,x2,y1,y2, must be set by finding the maximum log log window
            #     'legendtitle': '$\\xi_\\tau$',
            # },
        },
        'opdm-gapsize-d-L': {
            'default': default,
            'include': {
                'L': [12, 14],
                'g': [0.000, 0.005, 0.010, 0.015, 0.020, 0.025, 0.030],
                'd': [-6, -5, -4, -3, -2, -1, 0, 1, 2, 3 ,4 ,5 ,6],
                'ø': [None], # Dummy variable for subspec
            },
            'groupbase': 'tables',
            'groupname': 'opdm_spectrum',
            'colname': 'eigval',
            'figspec': ['g'],
            'subspec': ['ø'],
            'linspec': ['L'],
            'xaxspec': ['d'],
            'plotdir': Path(plotdir, Path(mplstyle).stem),
            'filename': 'opdm-gapsize-d',
            'palettes': ["viridis_r"],
            'mplstyle': mplstyle,
            # 'titlename': 'l-bit decay fit $C e^{-(|i-j|/\\xi)^\\beta}$',
            # 'filter': {
            # 'L': [24],
            # 'f': [0.2,0.4],
            # 'u': [16],
            # },
            # 'titlename': 'l-bit decay fit $C e^{-|i-j|/\\xi}$',
            # 'ylabel': '$\log_{10} \langle \langle O(|i-j|) \\rangle\\rangle$ ',
            'ylabel': '$1-$ opdm gap',
            # 'ylabel': '$\log_{10} \\bar O(|i-j|)$ ',
            'xlabel': "$\Delta$",
            # 'xlabelpad': -8,
            # 'xcoords': (0.5, -0.04),
            # 'xticks': [0, 0.25, 0.5, 0.75, 1.0] if prl else None,
            # 'xticks': [0, 15],
            # 'xticklabels': ["0","27"],
            # 'yticks': [1e-2, 1e-12],
            # 'ylabelpad': -16,

            # 'yticklabels': ['$-1$','$-15$'],
            # 'ycoords': (-0.34, 0.34),
            # 'xticklabels': ['$0$', '$L/2$', '$L$'],
            # 'yticks': [0, -3, -6, -9, -12, -15],
            # 'yscale': 'log',
            # 'ynopos': 'mask',

            # 'ymin': -16,
            # 'ymin': 0,
            # 'ymax': 2,
            # 'xmin': 0,
            # 'xmax': 15,
            # 'xnormalize': False,
            # 'xshift2mid': False,
            # 'xmax': 32 if prl else None,
            # 'ymin': 1e-14,
            # 'legendcols': ['f', 'tstd', 'tgw8', 'cstd', 'cgw8', 'ubond'],  # Choose 'num', 'bmax','tsim'
            'legendcols': ['L', 'g'],  # Choose 'num', 'bmax','tsim'
            # 'legendfits': ['xi', 'beta'] if prl else ['C', 'xi', 'beta', 'pos'],
            'legendoutside': False,
            'legendcollect': False,
            # 'legendlocation': (0.18, 0.0),
            'legendlocation': ['upper right', 'upper right', 'upper right', 'upper right', ],
            # 'bbox_to_anchor': (1.0, 0.70),  # Use with loc 'upper right'
            'bbox_to_anchor': (1.00, 0.9),
            'frameon': False,
            # 'legendtitle': 'Arithmetic average',
            'legendtitle': None,  # '$\\overline O \propto e^{(\\frac{|i-j|}{\\xi_\\tau})^\\beta}$',
            # 'legendtitle': '$y = C e^{-|i-j|/\\xi_\\tau}$',
            # 'legendtitle': '$\log \\bar O(x) = a - x \\xi_\\tau^{-1}$',
            # 'lbit-site': [0, 'mid', 'last'],
            # 'lbit-site': ['mid'],
            # 'lbit-mean': 'arithmetic',
            # 'lbit-axis': '',
            # 'fit-beta': True,
            # 'fit-ymin': 1e-16,
            # 'fit-skip': 0 if prl else 0,
            # 'fit-mark': False,
            # 'fit-plot': False,
            # 'inset-cls': {
            #     # 'pos': [0.03, 0.6, 0.40, 0.40], # Positon of the inset, x0 y0 width height
            #     'pos': [0.17, 0.15, 0.25, 0.25],  # Positon of the inset, x0 y0 width height
            #     'coords': [None, None, None, None],
            #     # These zoom limits x1,x2,y1,y2, must be set by finding the maximum log log window
            #     'legendtitle': '$\\xi_\\tau$',
            # },
        },
        'opdm-gapsize-d-g': {
            'default': default,
            'include': {
                'L': [12, 14],
                'g': [0.000, 0.005, 0.010, 0.015, 0.020, 0.025, 0.030],
                'd': [-6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6],
                'ø': [None],  # Dummy variable for subspec
            },
            'groupbase': 'tables',
            'groupname': 'opdm_spectrum',
            'colname': 'eigval',
            'figspec': ['L'],
            'subspec': ['ø'],
            'linspec': ['g'],
            'xaxspec': ['d'],
            'plotdir': Path(plotdir, Path(mplstyle).stem),
            'filename': 'opdm-gapsize-d',
            'palettes': ["viridis_r"],
            'mplstyle': mplstyle,
            # 'titlename': 'l-bit decay fit $C e^{-(|i-j|/\\xi)^\\beta}$',
            # 'filter': {
            # 'L': [24],
            # 'f': [0.2,0.4],
            # 'u': [16],
            # },
            # 'titlename': 'l-bit decay fit $C e^{-|i-j|/\\xi}$',
            # 'ylabel': '$\log_{10} \langle \langle O(|i-j|) \\rangle\\rangle$ ',
            'ylabel': '$1-$ opdm gap',
            # 'ylabel': '$\log_{10} \\bar O(|i-j|)$ ',
            'xlabel': "$\Delta$",
            # 'xlabelpad': -8,
            # 'xcoords': (0.5, -0.04),
            # 'xticks': [0, 0.25, 0.5, 0.75, 1.0] if prl else None,
            # 'xticks': [0, 15],
            # 'xticklabels': ["0","27"],
            # 'yticks': [1e-2, 1e-12],
            # 'ylabelpad': -16,

            # 'yticklabels': ['$-1$','$-15$'],
            # 'ycoords': (-0.34, 0.34),
            # 'xticklabels': ['$0$', '$L/2$', '$L$'],
            # 'yticks': [0, -3, -6, -9, -12, -15],
            # 'yscale': 'log',
            # 'ynopos': 'mask',

            # 'ymin': -16,
            # 'ymin': 0,
            # 'ymax': 2,
            # 'xmin': 0,
            # 'xmax': 15,
            # 'xnormalize': False,
            # 'xshift2mid': False,
            # 'xmax': 32 if prl else None,
            # 'ymin': 1e-14,
            # 'legendcols': ['f', 'tstd', 'tgw8', 'cstd', 'cgw8', 'ubond'],  # Choose 'num', 'bmax','tsim'
            'legendcols': ['L', 'g'],  # Choose 'num', 'bmax','tsim'
            # 'legendfits': ['xi', 'beta'] if prl else ['C', 'xi', 'beta', 'pos'],
            'legendoutside': False,
            'legendcollect': False,
            # 'legendlocation': (0.18, 0.0),
            'legendlocation': ['upper right', 'upper right', 'upper right', 'upper right', ],
            # 'bbox_to_anchor': (1.0, 0.70),  # Use with loc 'upper right'
            'bbox_to_anchor': (1.00, 0.9),
            'frameon': False,
            # 'legendtitle': 'Arithmetic average',
            'legendtitle': None,  # '$\\overline O \propto e^{(\\frac{|i-j|}{\\xi_\\tau})^\\beta}$',
            # 'legendtitle': '$y = C e^{-|i-j|/\\xi_\\tau}$',
            # 'legendtitle': '$\log \\bar O(x) = a - x \\xi_\\tau^{-1}$',
            # 'lbit-site': [0, 'mid', 'last'],
            # 'lbit-site': ['mid'],
            # 'lbit-mean': 'arithmetic',
            # 'lbit-axis': '',
            # 'fit-beta': True,
            # 'fit-ymin': 1e-16,
            # 'fit-skip': 0 if prl else 0,
            # 'fit-mark': False,
            # 'fit-plot': False,
            # 'inset-cls': {
            #     # 'pos': [0.03, 0.6, 0.40, 0.40], # Positon of the inset, x0 y0 width height
            #     'pos': [0.17, 0.15, 0.25, 0.25],  # Positon of the inset, x0 y0 width height
            #     'coords': [None, None, None, None],
            #     # These zoom limits x1,x2,y1,y2, must be set by finding the maximum log log window
            #     'legendtitle': '$\\xi_\\tau$',
            # },
        },

        'var-dist-g': {
            'default': default2x2,
            'include': {
                'L': [12, 14],
                'g': [0.0, 0.01, 0.02, 0.03],
                'd': [0],
            },
            'groupbase': 'tables',
            'groupname': 'measurements',
            'dsetname': 'data',
            'colname': 'energy_variance',
            'exactKeys': True,  # Match this column name exactly
            'figspec': ['L'],
            'subspec': ['g'],
            'linspec': ['d'],
            'plotdir': Path(plotdir, Path(mplstyle).stem),
            'filename': 'energy-variance-distribution-g',
            'palettes': ["Dark2"],
            'mplstyle': mplstyle,
            # 'titlename': 'l-bit decay fit $C e^{-(|i-j|/\\xi)^\\beta}$',
            # 'filter': {
            # 'L': [24],
            # 'f': [0.2,0.4],
            # 'u': [16],
            # },
            # 'titlename': 'l-bit decay fit $C e^{-|i-j|/\\xi}$',
            # 'ylabel': '$\log_{10} \langle \langle O(|i-j|) \\rangle\\rangle$ ',
            'ylabel': 'Histogram',
            # 'ylabel': '$\log_{10} \\bar O(|i-j|)$ ',
            'xlabel': "$\mathrm{Var}(H)$",
            # 'xlabelpad': -8,
            # 'xcoords': (0.5, -0.04),
            # 'xticks': [0, 0.25, 0.5, 0.75, 1.0] if prl else None,
            # 'xticks': [0, 15],
            # 'xticklabels': ["0","27"],
            # 'yticks': [1e-2, 1e-12],
            # 'ylabelpad': -16,

            # 'yticklabels': ['$-1$','$-15$'],
            # 'ycoords': (-0.34, 0.34),
            # 'xticklabels': ['$0$', '$L/2$', '$L$'],
            # 'yticks': [0, -3, -6, -9, -12, -15],
            'yscale': 'log',
            'xscale': 'log',
            # 'ynopos': 'mask',
            'xticks': [1e-15, 1e-10, 1e-5, 1e-0],
            'bins': np.logspace(start=-18, stop=0, num=19, endpoint=True),
            'density': False,
            # 'ymin': -16,
            # 'ymin': 0,
            # 'ymax': 2,
            # 'xmin': 0,
            # 'xmax': 15,
            # 'xnormalize': False,
            # 'xshift2mid': False,
            # 'xmax': 32 if prl else None,
            # 'ymin': 1e-14,
            # 'legendcols': ['f', 'tstd', 'tgw8', 'cstd', 'cgw8', 'ubond'],  # Choose 'num', 'bmax','tsim'
            'legendcols': ['L', 'g', 'd'],  # Choose 'num', 'bmax','tsim'
            # 'legendfits': ['xi', 'beta'] if prl else ['C', 'xi', 'beta', 'pos'],
            'legendoutside': False,
            'legendcollect': False,
            # 'legendlocation': (0.18, 0.0),
            'legendlocation': ['upper right', 'upper right', 'upper right', 'upper right', ],
            # 'bbox_to_anchor': (1.0, 0.70),  # Use with loc 'upper right'
            'bbox_to_anchor': (1.00, 1.01),
            'frameon': False,
            # 'legendtitle': 'Arithmetic average',
            'legendtitle': None,  # '$\\overline O \propto e^{(\\frac{|i-j|}{\\xi_\\tau})^\\beta}$',
            # 'legendtitle': '$y = C e^{-|i-j|/\\xi_\\tau}$',
            # 'legendtitle': '$\log \\bar O(x) = a - x \\xi_\\tau^{-1}$',
            # 'lbit-site': [0, 'mid', 'last'],
            # 'lbit-site': ['mid'],
            # 'lbit-mean': 'arithmetic',
            # 'lbit-axis': '',
            # 'fit-beta': True,
            # 'fit-ymin': 1e-16,
            # 'fit-skip': 0 if prl else 0,
            # 'fit-mark': False,
            # 'fit-plot': False,
            # 'inset-cls': {
            #     # 'pos': [0.03, 0.6, 0.40, 0.40], # Positon of the inset, x0 y0 width height
            #     'pos': [0.17, 0.15, 0.25, 0.25],  # Positon of the inset, x0 y0 width height
            #     'coords': [None, None, None, None],
            #     # These zoom limits x1,x2,y1,y2, must be set by finding the maximum log log window
            #     'legendtitle': '$\\xi_\\tau$',
            # },
        },
        'var-dist-d': {
            'default': default3x3,
            'include': {
                'L': [12, 14],
                'g': [0.01],
                'd': [-3, -1, 0, -1, 3, 6],
            },
            'groupbase': 'tables',
            'groupname': 'measurements',
            'dsetname': 'data',
            'colname': 'energy_variance',
            'exactKeys': True,  # Match this column name exactly
            'figspec': ['L'],
            'subspec': ['d'],
            'linspec': ['g'],
            'plotdir': Path(plotdir, Path(mplstyle).stem),
            'filename': 'energy-variance-distribution-d',
            'palettes': ["Dark2"],
            'mplstyle': mplstyle,
            # 'titlename': 'l-bit decay fit $C e^{-(|i-j|/\\xi)^\\beta}$',
            # 'filter': {
            # 'L': [24],
            # 'f': [0.2,0.4],
            # 'u': [16],
            # },
            # 'titlename': 'l-bit decay fit $C e^{-|i-j|/\\xi}$',
            # 'ylabel': '$\log_{10} \langle \langle O(|i-j|) \\rangle\\rangle$ ',
            'ylabel': 'Histogram',
            # 'ylabel': '$\log_{10} \\bar O(|i-j|)$ ',
            'xlabel': "$\mathrm{Var}(H)$",
            # 'xlabelpad': -8,
            # 'xcoords': (0.5, -0.04),
            # 'xticks': [0, 0.25, 0.5, 0.75, 1.0] if prl else None,
            # 'xticks': [0, 15],
            # 'xticklabels': ["0","27"],
            # 'yticks': [1e-2, 1e-12],
            # 'ylabelpad': -16,

            # 'yticklabels': ['$-1$','$-15$'],
            # 'ycoords': (-0.34, 0.34),
            # 'xticklabels': ['$0$', '$L/2$', '$L$'],
            # 'yticks': [0, -3, -6, -9, -12, -15],
            'xticks': [1e-15, 1e-10, 1e-5, 1e-0],
            'yscale': 'log',
            'xscale': 'log',
            # 'ynopos': 'mask',
            'bins': np.logspace(start=-18, stop=0, num=19, endpoint=True),
            'density': False,
            # 'ymin': -16,
            # 'ymin': 0,
            # 'ymax': 2,
            # 'xmin': 0,
            # 'xmax': 15,
            # 'xnormalize': False,
            # 'xshift2mid': False,
            # 'xmax': 32 if prl else None,
            # 'ymin': 1e-14,
            # 'legendcols': ['f', 'tstd', 'tgw8', 'cstd', 'cgw8', 'ubond'],  # Choose 'num', 'bmax','tsim'
            'legendcols': ['L', 'g', 'd'],  # Choose 'num', 'bmax','tsim'
            # 'legendfits': ['xi', 'beta'] if prl else ['C', 'xi', 'beta', 'pos'],
            'legendoutside': False,
            'legendcollect': False,
            # 'legendlocation': (0.18, 0.0),
            'legendlocation': ['upper right', 'upper right', 'upper right', 'upper right', ],
            # 'bbox_to_anchor': (1.0, 0.70),  # Use with loc 'upper right'
            'bbox_to_anchor': (1.00, 1.01),
            'frameon': False,
            # 'legendtitle': 'Arithmetic average',
            'legendtitle': None,  # '$\\overline O \propto e^{(\\frac{|i-j|}{\\xi_\\tau})^\\beta}$',
            # 'legendtitle': '$y = C e^{-|i-j|/\\xi_\\tau}$',
            # 'legendtitle': '$\log \\bar O(x) = a - x \\xi_\\tau^{-1}$',
            # 'lbit-site': [0, 'mid', 'last'],
            # 'lbit-site': ['mid'],
            # 'lbit-mean': 'arithmetic',
            # 'lbit-axis': '',
            # 'fit-beta': True,
            # 'fit-ymin': 1e-16,
            # 'fit-skip': 0 if prl else 0,
            # 'fit-mark': False,
            # 'fit-plot': False,
            # 'inset-cls': {
            #     # 'pos': [0.03, 0.6, 0.40, 0.40], # Positon of the inset, x0 y0 width height
            #     'pos': [0.17, 0.15, 0.25, 0.25],  # Positon of the inset, x0 y0 width height
            #     'coords': [None, None, None, None],
            #     # These zoom limits x1,x2,y1,y2, must be set by finding the maximum log log window
            #     'legendtitle': '$\\xi_\\tau$',
            # },
        },

        'ent': {
            'default': default,
            'groupbase': 'tables',
            'groupname': 'measurements',
            'colname': 'entanglement_entropy',
            'normpage': False,
            'normsinf': False,
            'ylabel': '$\overline S_\mathrm{E}$',
            'xlabel': '$t$',
            'xticks': [1e1, 1e4, 1e7, 1e10, 1e13],
            'xscale': 'log',
            # 'xticks': [],
            'xmin': 1e-1,
            'xmax': 1e+15,
            'ymin': 0,
            'ymax': 1.73,
            # 'yscale': 'log',
            'yformat': '%.2f',
            'plotprefix': 'SE',
            'plotdir': Path(plotdir, Path(mplstyle).stem),
            'filename': 'SE',
            'markerlist': ['growth-begin', 'num-lnlnt-cease', 'num-saturated', 'ent-lnt-cease', 'ent-saturated'],
            'findsaturation': True,  # Instead of taking the last value, take the average of the plateau
            'marklogwindow': True,
            'marksaturation': True,
            'findloglogwindow': True,
            'markloglogwindow': True,
            'markloglogoffset': 0.1,
            'fillerror': True,
            'markerror': False,
            'timeselection': 'lnt',
            'timenormalization': '',
            'mplstyle': mplstyle,
            'legendcols': ['L', 'l', 'num'],  # Choose 'num', 'bmax','tsim'
            'legendlocation': ['upper right', 'upper right', 'upper right', 'upper right', ],
            # 'bbox_to_anchor': (1.0, 0.70),  # Use with loc 'upper right'
            'bbox_to_anchor': (1.00, 1.00),

            # 'legendcols': ['L', 'f','num'],  # Choose 'num', 'bmax','tsim'
            # 'legendcols': ['f', 'x', 'num', 'bmax:.0f', 'bavg:.0f', 'tsim'],  # Choose 'num', 'bmax','tsim'
            # 'legendlocation': [ (0.75, 0.42), (0.75, 0.42), (0.75, 0.42), (0.75, 0.42), ],
            # 'legendlocation': ['upper right', 'upper right', 'upper right', 'upper right', ],
            # 'bbox_to_anchor': (1.05, 1.05),  # Use with loc 'upper right'
        },

    }

    return meta
