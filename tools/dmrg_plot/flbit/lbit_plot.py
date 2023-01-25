import glob
import os.path

import matplotlib.pyplot as plt
from dmrg_plot.common.io.h5ops import *
from dmrg_plot.common.io.parse import parse

from database.database import *
from lbit_avg import lbit_avg
from plotting.meta import *
from plotting.multiplot import *


def lbit_plot(args):
    projdir = args.basedir
    algo_filter = args.algos
    state_filter = args.states
    model_filter = args.models
    point_filter = args.points
    batches = args.batches
    h5avgs = []
    dbs = []
    plotdir = None
    meta = None
    for batch in batches:
        batchglob = glob.glob('{}/{}*'.format(projdir, batch))
        if batchglob:
            print(batchglob)
            batchdir = batchglob[0]
            avgfile = '{}/analysis/data/averaged.h5'.format(batchdir)
            plotdir = '{}/analysis/plots'.format(batchdir)
            # if meta is None:
            meta = get_meta(plotdir)
            h5avgs.append(h5py.File(avgfile, 'r'))
            dbs.append(load_time_database2(h5avgs[-1], meta, algo_filter=algo_filter, model_filter=model_filter,
                                           state_filter=state_filter, debug=False))

    if not os.path.exists(plotdir):
        os.makedirs(plotdir)

    palettes = [  # Distinguishable colors
        "tab10",
        "husl",
        "YlOrBr",
        "Set2",
        "magma",
    ]

    palettes = [  # Sequential colors
        "summer_r",
        "magma_r",
        "winter_r",
    ]

    # fig3 = ['J', 'x:.2f', 'r']
    # sub3 = ['L', 'u', 'w']
    # l1   = ['f']

    fig3 = ['J', 'x:.2f', 'r']
    sub3 = ['f', 'u', 'w']
    sub2 = ['f', 'u']
    l1 = ['L']

    # #
    # fig3 = ['J', 'x:.2f', 'r']
    # sub3 = ['L', 'f', 'w']
    # l1   = ['u']
    # MAKE PLOTS OF FIGURE 4 IN https://journals.aps.org/prb/abstract/10.1103/PhysRevB.103.024203

    # f = None
    # for db, palette in zip(dbs, palettes):
    #     f = plot_v2_cls_fig3_sub3_line1(db=db, meta=meta['cls'], figspec=['J', 'r', 'w'], subspec=['L', 'u', 'x:.2f'], linspec=['f'], figs=f, palette_name=palette)
    # save_figure(f)
    #

    # plt.show()
    # exit(0)

    # f = None
    # for db, palette in zip(dbs, palettes):
    #     f = plot_v2_time_fig3_sub3_line1(db=db, meta=meta['ent'], figspec=fig3, subspec=sub3, linspec=l1, figs=f, palette_name=palette)
    # save_figure(f)
    #
    # f = None
    # for db, palette in zip(dbs, palettes):
    #     f = plot_v2_time_fig3_sub3_line1(db=db, meta=meta['num1'], figspec=fig3, subspec=sub3, linspec=l1, figs=f, palette_name=palette)
    # save_figure(f)
    #
    # f = None
    # for db, palette in zip(dbs, palettes):
    #     f = plot_v2_time_fig3_sub3_line1(db=db, meta=meta['num2'], figspec=fig3, subspec=sub3, linspec=l1, figs=f, palette_name=palette)
    # save_figure(f)
    #
    # f = None
    plot_divg_fig3_sub3_line1.prb = None
    # for db, palette in zip(dbs, palettes):
    #     f = plot_divg_fig3_sub3_line1(db=db, meta=meta['divg-num'], figspec=fig3, subspec=sub3, linspec=l1, figs=f, palette_name=palette)
    # save_figure(f)

    f = None
    for db, palette in zip(dbs, palettes):
        f = plot_v2_lbit_fig3_sub3_line1(db=db, meta=meta['lbit'], figspec=['J', 'r', 'f'], subspec=['L', 'u', 'x:.2f'],
                                         linspec=['w'], figs=f,
                                         palette_name=palette)
    save_figure(f)

    fig3 = ['J', 'x:.2f', 'r']
    sub3 = ['L', 'u', 'f']
    sub2 = ['f', 'u']
    l1 = ['w']

    f = None
    for db, palette in zip(dbs, palettes):
        f = plot_v2_time_fig3_sub3_line1(db=db, meta=meta['ent'], figspec=fig3, subspec=sub3, linspec=l1, figs=f,
                                         palette_name=palette)
    save_figure(f)

    f = None
    for db, palette in zip(dbs, palettes):
        f = plot_v2_time_fig3_sub3_line1(db=db, meta=meta['num1'], figspec=fig3, subspec=sub3, linspec=l1, figs=f,
                                         palette_name=palette)
    save_figure(f)

    f = None
    for db, palette in zip(dbs, palettes):
        f = plot_v2_time_fig3_sub3_line1(db=db, meta=meta['num2'], figspec=fig3, subspec=sub3, linspec=l1, figs=f,
                                         palette_name=palette)
    save_figure(f)

    f = None
    plot_divg_fig3_sub3_line1.prb = None
    for db, palette in zip(dbs, palettes):
        f = plot_divg_fig3_sub3_line1(db=db, meta=meta['divg-num'], figspec=fig3, subspec=sub3, linspec=l1, figs=f,
                                      palette_name=palette)
    save_figure(f)

    plt.show()
    exit(0)

    f = None
    for db, palette in zip(dbs, palettes):
        f = plot_divg_fig3_sub3_line1(db=db, meta=meta['divg-num'], figspec=fig3, subspec=sub3, linspec=l1, figs=f, palette_name=palette)
    save_figure(f)

    f = None
    for db, palette in zip(dbs, palettes):
        f = plot_divg_fig3_sub3_line1(db=db, meta=meta['divg-ent'], figspec=fig3, subspec=sub3, linspec=l1, figs=f, palette_name=palette)
    save_figure(f)

    f = None
    for db, palette in zip(dbs, palettes):
        f = plot_tavg_fig3_sub2_line1(db=db, meta=meta['tavg-ent'], figspec=fig3, subspec=['w'], linspec=['f', 'u'], x1=['L'], figs=f, palette_name=palette)
    save_figure(f)

    f = None
    for db, palette in zip(dbs, palettes):
        f = plot_tavg_fig3_sub2_line1(db=db, meta=meta['tavg-num'], figspec=fig3, subspec=['w'], linspec=['f', 'u'], x1=['L'], figs=f, palette_name=palette)
    save_figure(f)

    f = None
    for db, palette in zip(dbs, palettes):
        f = plot_dist_fig3_sub3_line1(db=db, meta=meta['dist-chi'], figspec=fig3, subspec=sub3, linspec=l1, figs=f, palette_name=palette)
    save_figure(f)

    f = None
    for db, palette in zip(dbs, palettes):
        f = plot_dist_fig3_sub3_line1(db=db, meta=meta['dist-num'], figspec=fig3, subspec=sub3, linspec=l1, figs=f, palette_name=palette)
    save_figure(f)

    f = None
    for db, palette in zip(dbs, palettes):
        f = plot_divg_fig3_sub3_line1(db=db, meta=meta['divg-num'], figspec=fig3, subspec=sub3, linspec=l1, figs=f, palette_name=palette)
    save_figure(f)
    f = None
    for db, palette in zip(dbs, palettes):
        f = plot_v2_time_fig3_sub3_line1(db=db, meta=meta['trn'], figspec=fig3, subspec=sub3, linspec=l1, figs=f, palette_name=palette)
    save_figure(f)

    f = None
    for db, palette in zip(dbs, palettes):
        f = plot_v2_time_fig3_sub3_line1(db=db, meta=meta['chi'], figspec=fig3, subspec=sub3, linspec=l1, figs=f, palette_name=palette)
    save_figure(f)

    f = None
    for db, palette in zip(dbs, palettes):
        f = plot_v2_time_fig3_sub3_line1(db=db, meta=meta['ent'], figspec=fig3, subspec=sub3, linspec=l1, figs=f, palette_name=palette)
    save_figure(f)

    f = None
    for db, palette in zip(dbs, palettes):
        f = plot_v2_time_fig3_sub3_line1(db=db, meta=meta['num1'], figspec=fig3, subspec=sub3, linspec=l1, figs=f, palette_name=palette)
    save_figure(f)

    f = None
    for db, palette in zip(dbs, palettes):
        f = plot_v2_time_fig3_sub3_line1(db=db, meta=meta['num2'], figspec=fig3, subspec=sub3, linspec=l1, figs=f, palette_name=palette)
    save_figure(f)

    f = None
    for db, palette in zip(dbs, palettes):
        f = plot_v2_time_fig3_sub3_line1(db=db, meta=meta['numH1'], figspec=fig3, subspec=sub3, linspec=l1, figs=f, palette_name=palette)
    save_figure(f)

    f = None
    for db, palette in zip(dbs, palettes):
        f = plot_v2_time_fig3_sub3_line1(db=db, meta=meta['numH2'], figspec=fig3, subspec=sub3, linspec=l1, figs=f, palette_name=palette)
    save_figure(f)

    f = None
    for db, palette in zip(dbs, palettes):
        f = plot_v2_time_fig3_sub3_line1(db=db, meta=meta['numa'], figspec=fig3, subspec=sub3, linspec=l1, figs=f, palette_name=palette)
    save_figure(f)

    f = None
    for db, palette in zip(dbs, palettes):
        f = plot_v2_time_fig3_sub3_line1(db=db, meta=meta['numHa'], figspec=fig3, subspec=sub3, linspec=l1, figs=f, palette_name=palette)
    save_figure(f)

    f = None
    for db, palette in zip(dbs, palettes):
        f = plot_v2_time_fig3_sub3_line1(db=db, meta=meta['nument1'], figspec=fig3, subspec=sub3, linspec=l1, figs=f, palette_name=palette)
    save_figure(f)

    f = None
    for db, palette in zip(dbs, palettes):
        f = plot_v2_time_fig3_sub3_line1(db=db, meta=meta['nument2'], figspec=fig3, subspec=sub3, linspec=l1, figs=f, palette_name=palette)
    save_figure(f)

    f = None
    for db, palette in zip(dbs, palettes):
        f = plot_v2_time_fig3_sub3_line1(db=db, meta=meta['tsim'], figspec=fig3, subspec=sub3, linspec=l1, figs=f, palette_name=palette)
    save_figure(f)

    f = None
    for db, palette in zip(dbs, palettes):
        f = plot_v2_time_fig3_sub3_line1(db=db, meta=meta['titr'], figspec=fig3, subspec=sub3, linspec=l1, figs=f, palette_name=palette)
    save_figure(f)

    f = None
    for db, palette in zip(dbs, palettes):
        f = plot_v2_lbit_fig3_sub3_line1(db=db, meta=meta['lbit'], figspec=['J', 'r', 'f'], subspec=['w', 'u', 'x:.2f'], linspec=['L'], figs=f,
                                         palette_name=palette)
    save_figure(f)

    plt.show()
    exit(0)


if __name__ == '__main__':
    args = parse('fLBIT')
    lbit_avg(args)
    lbit_plot(args)
