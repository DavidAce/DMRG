from __future__ import unicode_literals
import errno
import matplotlib.pyplot as plt
import os.path
import numpy as np
import h5py
from src.plotting.style import *
from src.general.natural_sort import *
from src.plotting.entanglement_entropy_chain import *
from src.io.h5ops import *
from src.plotting.multiplot import *
from src.database.database import *
from src.plotting.meta import *

import glob
import pickle
import hashlib
import json

# Adjust pycharm width:
desired_width = 320
np.set_printoptions(linewidth=desired_width)

markerlist = get_markerlist()

projdir = '/mnt/WDB-AN1500/mbl_transition'
batchdir = glob.glob(projdir + '/data170*')[0]
# suffix = '_new'
# suffix = '-data152'
suffix = ''

analysisdir = batchdir + '/analysis'
datadir = analysisdir + '/data'
plotdir = '{}/plots{}-weekly'.format(analysisdir, suffix)
avgfile = '{}/averaged{}.h5'.format(datadir, suffix)
# avgfile = datadir + '/averaged-res-multi6.h5'


if not os.path.exists(analysisdir):
    os.makedirs(analysisdir)
if not os.path.exists(datadir):
    os.makedirs(datadir)
if not os.path.exists(plotdir):
    os.makedirs(plotdir)
vwin = [0, np.inf]

h5_avg = h5open(avgfile, 'r')
meta = get_meta(plotdir)

# Pickle
dbhash = hashlib.md5(json.dumps(meta['common']['include'], sort_keys=True).encode('utf-8')).hexdigest()
pkfile = '{}/h5_avg_{}.pk'.format(datadir, dbhash)
db = None

if os.path.isfile(pkfile):
    try:
        with open(pkfile, 'rb') as pk:
            db = pickle.load(pk)
    except:
        db = None
        pass

if db is None:
    print("No previous .pk file found")
    with open(pkfile, 'wb') as pk:
        db = load_database(h5_avg, meta, debug=False)
        # dump your data into the file
        pickle.dump(db, pk)

# plot_table(h5_avg, db, meta['fes-ent'], sub1=['L:.0f'], l1=['l:.2f'])
plot_table(h5_avg, db, meta['chain-ent'], sub1=['L:.0f'], l1=['l:.2f'])
# plot_table(h5_avg, db, meta['chain-bond'], sub1=['L:.0f'], l1=['l:.2f'])
# plot_histogram(h5_avg, db, meta['hist-var'], sub1=['L:.0f'], l1=['l:.2f'])
# plot_slice(h5_avg, db, meta['mid-ent'], sub1=['d'], l1=['l'], x1=['L'])


# plot_midchain(h5_avg, db, meta['ent'], sub1=['l'], l1=['d'], x1=['L'])
# plot_midchain(h5_avg, db, meta['ent'], sub1=['d'], l1=['l'], x1=['L'])
# plot_midchain(h5_avg, db, meta['ent'], sub1=['d'], l1=['L'], x1=['l'])
# plot_midchain(h5_avg, db, meta['ent'], sub1=['L'], l1=['l'], x1=['d'])
# plot_midchain(h5_avg, db, meta['ent'], sub1=['L'], l1=['d'], x1=['l'])

# plot_midchain(h5_avg, db, meta['fesent'], sub1=['d:.2f'], l1=['l:.2f'], x1=['b:.0f'])
# plot_midchain(h5_avg, db, meta['fesent'], sub1=['L:.0f'], l1=['l:.2f'], x1=['b:.0f'])
# plot_midchain(h5_avg, db, meta['fesent'], sub1=['l:.2f'], l1=['L:.0f'], x1=['b:.0f'])  #  Shows S = log(b) ?


# plot_measurement(h5_avg, db, meta['var'], sub1=['L:.0f'], l1=['d:.2f'], x1=['l:.2f'])
# plot_measurement(h5_avg, db, meta['var'], sub1=['d:.2f'], l1=['L:.0f'], x1=['l:.2f'])
# plot_measurement(h5_avg, db, meta['var'], sub1=['l:.2f'], l1=['L:.0f'], x1=['b:.0f'])
# plot_measurement(h5_avg, db, meta['var'], sub1=['L'], l1=['l'], x1=['b'])
# plot_measurement(h5_avg, db, meta['time'], sub1=['d'], l1=['l'], x1=['L'])
# plot_measurement(h5_avg, db, meta['time'], sub1=['d'], l1=['L'], x1=['l'])
#
# multiplot_var_distribution     (h5_avg, db=db, meta=meta['var'], vwin=vwin)
# multiplot_time_distribution    (h5_avg, db=db, meta=meta['time'])


# multiplot_S_vs_Site_fig1_sub1_l1 (h5_avg, db, meta = meta['ent'], plotdir=plotdir, algo_inc=algo_inc_xed, state_inc=state_inc_xed, fig1=['L'], sub1=['d'], l1=['l'], vwin=vwin)
# multiplot_S_vs_Site_fig1_sub2_l3 (h5_avg, db_SE,meta = meta['ent'], plotdir=plotdir, algo_inc=algo_inc_xed, state_inc=state_inc_xed, fig1=['l'], sub2=['L','d'], vwin=vwin)
# multiplot_S_vs_Site_fig1_sub2_l3 (h5_avg, db_SE,meta = meta['ent'], plotdir=plotdir, algo_inc=algo_inc_xed, state_inc=state_inc_xed, fig1=['l'], sub2=['L','d'], vwin=vwin)
# multiplot_Smid_sub1_l1_x1 (h5_avg, db_SE, meta = meta['ent'], plotdir=plotdir, algo_inc=algo_inc_xed, state_inc=state_inc_xed, sub1=['L'], l1=['l'], x1=['d'], vwin=vwin)
# multiplot_Smid_vs_Delta     (h5_avg, db_SE, plotdir=plotdir, algo_inc=algo_inc_xed, state_inc=state_inc_xed)
# multiplot_Smid_vs_tmax_v2   (h5_avg, db_SE, plotdir=plotdir, algo_inc=algo_inc_xed, state_inc=state_inc_xed)
# multiplot_S_vs_Length       (avgfile, plotdir=plotdir, algo_inc=algo_inc_xed, state_inc=state_inc_xed)
# multiplot_S_distribution       (avgfile, plotdir=plotdir, algo_inc=algo_inc_xed, state_inc=state_inc_xed, bins=40)
# multiplot_energy_distribution  (avgfile, plotdir=plotdir, algo_inc=algo_inc_xed, state_inc=state_inc_xed)
# multiplot_chi_distribution     (avgfile, plotdir=plotdir, algo_inc=algo_inc_x, state_inc=state_inc_x)

h5close(h5_avg)
plt.show()
exit(0)

#
multiplot_S_distribution_diff(avgfile, plotdir=plotdir, algo_inc=algo_inc_x, state_inc=state_inc_x, key_comp='ed-e0.00', normalized=True, bins=25)

# multiplot_Sq_vs_Site           (avgfile, plotdir=plotdir, algo_inc=algo_inc_xed, state_inc=state_inc_xed, q=2)
# multiplot_Sq_vs_Site           (avgfile, plotdir=plotdir, algo_inc=algo_inc_xed, state_inc=state_inc_xed, q=3)
# multiplot_Sq_vs_Site           (avgfile, plotdir=plotdir, algo_inc=algo_inc_xed, state_inc=state_inc_xed, q=4)
# multiplot_Sq_vs_Site           (avgfile, plotdir=plotdir, algo_inc=algo_inc_xed, state_inc=state_inc_xed, q=100)
# multiplot_schmidt              (avgfile, plotdir=plotdir, algo_inc=algo_inc_x, state_inc=state_inc_x)
# multiplot_iter_distribution     (avgfile, plotdir=plotdir, algo_inc=algo_inc_x, state_inc=state_inc_x)
# multiplot_mem_distribution     (avgfile, plotdir=plotdir, algo_inc=algo_inc_x, state_inc=state_inc_x)

# multiplot_Energy_distribution  (h5_filename,plotdir=plotdir,algo_inc=['xDMRG','ed'], state_inc=['state_','states'])
# multiplot_S_vs_ewin_scatter    (h5_filename,plotdir=plotdir,key_list=['orig','proj','ed','chi'])
# multiplot_S_vs_eps_scatter     (h5_filename,plotdir=plotdir,key_list=['orig'])
# multiplot_Sq_vs_Length         (h5_filename,plotdir=plotdir,key_list=['orig','ed-energy-0.00'])
# multiplot_delta_distribution   (h5_filename,plotdir=plotdir,key_list=['orig','proj','ed','chi'])
# multiplot_S_vs_Delta           (h5_filename,plotdir=plotdir,key_list=['orig','proj','ed','chi'])
# multiplot_Sq_vs_Site_log2      (h5_filename,plotdir=plotdir,key_list=['orig','proj','ed-energy-0.00'])

# times = [1,2,3,4,5,6,7,8,10,12,14,16,18,20,24,28,32,36,40,48,58,64,70]
# times = [1,2,3,4,5,6,7,8,10,12,14,16,18,20,24,28,32,36,40,48,58,64,70,80,100,120,140,160,180,200,250,300,350,400,550,600,700,800,900,1200,1600,2000]
# for num, time in enumerate(times):
#     multiplot_S_vs_Site_over_time     (h5_filename,plotdir=plotdir,algo_inc=['xDMRG','ed-e0.00'], state_inc=['state_','states'],time=[0,time],time_num=num)

h5close(h5_avg)
plt.show()
exit(0)

# multiplot_par_distribution(h5_filename,plotdir)

# multiplot_var_vs_Delta(h5_filename,plotdir, type='average' )
# multiplot_chi_vs_Delta(h5_filename,plotdir,type='average' )
# multiplot_time_vs_Delta(h5_filename,plotdir,type='average' )
# multiplot_e_vs_delta_and_lambda_scatter(h5_filename,plotdir)
# multiplot_e_vs_delta_and_lambda_probdist(h5_filename,plotdir)


# plt.show()
# exit(0)
# multiplot_S_vs_Delta(h5_filename,plotdir,type='average')


# plt.show()
# exit(0)


# plot_and_histogram_typical_var_vs_Delta
# plot_and_histogram_typical_entropy_vs_Delta

h5close(h5_avg)
plt.show()
exit(0)
