import glob
import os.path

import matplotlib.pyplot as plt
from dmrg_plot.common.io.h5ops import *
from dmrg_plot.common.io.parse import parse

from database.database import *
from lbit_avg import lbit_avg
from plotting.meta80 import *
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
    metas = None
    for batch in batches:
        batchglob = glob.glob('{}/{}*'.format(projdir, batch))
        batchnum = int(''.join(i for i in batch if i.isdigit()))
        version2 = batchnum <= 59
        version3 = batchnum >= 60
        if batchglob:
            print(batchglob)
            batchdir = batchglob[0]
            avgfile = '{}/analysis/data/averaged.h5'.format(batchdir)
            plotdir = '{}/analysis/plots'.format(batchdir)
            # if meta is None:
            # metas = [get_meta(plotdir), get_meta(plotdir), get_meta(plotdir)]
            # metas[1]['common']['include_vals']['tgw8'] = ['ID']
            # metas[2]['common']['include_vals']['cgw8'] = ['ID']
            metas = [get_meta(plotdir)]

            for meta in metas:
                h5avgs.append(h5py.File(avgfile, 'r'))
                if version2:
                    dbs.append(load_time_database2(h5avgs[-1], meta, algo_filter=algo_filter, model_filter=model_filter,
                                                   state_filter=state_filter, debug=False))
                if version3:
                    dbs.append(load_time_database3(h5avgs[-1], meta, algo_filter=algo_filter, model_filter=model_filter,
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
        "Blues",
        "Oranges",
        "Greens",
        "Purples"
    ]
    palettes = [  # Sequential2 colors
        "summer_r",
        "autumn_r",
        "winter_r",
        "spring_r"
    ]




    figspec_x = ['J', 'w', 'r',  'tstd', 'cstd', 'cgw8', 'tgw8']
    subspec_x = ['u']
    linspec_x = ['f']
    xaxspec_x = ['L']
    figspec = ['w', 'J', 'r', 'cgw8', 'tgw8', 'tstd', 'cstd']
    subspec = ['f', 'u']
    linspec = ['L']

    f = None
    for db, meta, palette in zip(dbs, metas, palettes):
        f = plot_lbit_fig_sub_line(db=db, meta=meta['lbit-avg'], figspec=figspec, subspec=subspec, linspec=linspec, figs=f,
                                   palette_name=palette)
    save_figure(f)
    plt.show()
    exit(0)
    f = None
    for db, meta, palette in zip(dbs, metas, palettes):
        f = plot_rise_fig_sub_line(db=db, meta=meta['rise-num2'], figspec=figspec_x, subspec=subspec_x,
                                   linspec=linspec_x, xaxspec=xaxspec_x, figs=f,
                                   palette_name=palette)
    save_figure(f)



    f = None
    for db, meta, palette in zip(dbs, metas, palettes):
        f = plot_divg_fig_sub_line(db=db, meta=meta['divg-num'], figspec=figspec, subspec=subspec, linspec=linspec,
                                   figs=f,
                                   palette_name=palette)
    save_figure(f)

    f = None
    for db, meta, palette in zip(dbs, metas, palettes):
        f = plot_time_fig_sub_line(db=db, meta=meta['ent'], figspec=figspec, subspec=subspec, linspec=linspec, figs=f,
                                   palette_name=palette)
    save_figure(f)

    f = None
    for db, meta, palette in zip(dbs, metas, palettes):
        f = plot_time_fig_sub_line(db=db, meta=meta['num1'], figspec=figspec, subspec=subspec, linspec=linspec, figs=f,
                                   palette_name=palette)
    save_figure(f)

    f = None
    for db, meta, palette in zip(dbs, metas, palettes):
        f = plot_time_fig_sub_line(db=db, meta=meta['num2'], figspec=figspec, subspec=subspec, linspec=linspec, figs=f,
                                   palette_name=palette)
    save_figure(f)

    f = None
    for db, meta, palette in zip(dbs, metas, palettes):
        f = plot_time_fig_sub_line(db=db, meta=meta['chi'], figspec=figspec, subspec=subspec, linspec=linspec, figs=f,
                                   palette_name=palette)
    save_figure(f)



    f = None
    for db, meta, palette in zip(dbs, metas, palettes):
        f = plot_slope_fig_sub_line(db=db, meta=meta['linearFit-num2'], figspec=figspec_x, subspec=subspec_x,
                                    linspec=linspec_x, xaxspec=xaxspec_x, figs=f,
                                    palette_name=palette)
    save_figure(f)

    f = None
    for db, meta, palette in zip(dbs, metas, palettes):
        f = plot_slope_fig_sub_line(db=db, meta=meta['linearFit-num2'], figspec=['L', 'tstd', 'cstd', 'cgw8', 'tgw8'],
                                    subspec=['J', 'w', 'r'],
                                    linspec=['f'], xaxspec=['u'], figs=f,
                                    palette_name=palette)
    save_figure(f)

    f = None
    for db, meta, palette in zip(dbs, metas, palettes):
        f = plot_svnt_fig_sub_line(db=db, meta=meta['num-svnt'], figspec=figspec, subspec=subspec, linspec=linspec,
                                   figs=f,
                                   palette_name=palette)
    save_figure(f)


    # f = None
    # for db, meta, palette in zip(dbs, metas, palettes):
    #     f = plot_lbit_fig_sub_line(db=db, meta=meta['lbit-typ'], figspec=['J', 'r', 'f'],
    #                                subspec=['w', 'x:.2f', 'ubond', 'L'],
    #                                linspec=['u'], figs=f,
    #                                palette_name=palette)
    # save_figure(f)
    # f = None
    # for db, meta, palette in zip(dbs, metas, palettes):
    #     f = plot_lbit_fig_sub_line(db=db, meta=meta['lbit-avg'], figspec=['J', 'w', 'r', 'ubond', 'f:.0f', 'x:.0f'],
    #                                subspec=['L'],
    #                                linspec=['u'], figs=f,
    #                                palette_name=palette)
    # save_figure(f)
    # f = None
    # for db, meta, palette in zip(dbs, metas, palettes):
    #     f = plot_cls_fig_sub_line(db=db, meta=meta['cls-avg'], figspec=figspec_x, subspec=subspec_x, linspec=linspec_x,
    #                               xaxspec=xaxspec_x, figs=f,
    #                               palette_name=palette)
    # save_figure(f)
    # f = None
    # for db, meta, palette in zip(dbs, metas, palettes):
    #     f = plot_cls_fig_sub_line(db=db, meta=meta['cls-avg'], figspec=figspec_x, subspec=subspec_x, linspec=['L'],
    #                               xaxspec=['u'], figs=f,
    #                               palette_name=palette)
    # save_figure(f)


    plt.show()
    exit(0)

    # PLOT SLOPE VS tstd/cstd RATIO

    f = None
    for db, palette in zip(dbs, palettes):
        f = plot_divg_fig3_sub3_line1(db=db, meta=meta['divg-num'], figspec=fig3, subspec=sub3, linspec=l1, figs=f,
                                      palette_name=palette)
    save_figure(f)

    f = None
    for db, palette in zip(dbs, palettes):
        f = plot_divg_fig3_sub3_line1(db=db, meta=meta['divg-ent'], figspec=fig3, subspec=sub3, linspec=l1, figs=f,
                                      palette_name=palette)
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
    args = parse('fLBIT', ['lbit80'])
    lbit_avg(args)
    lbit_plot(args)
