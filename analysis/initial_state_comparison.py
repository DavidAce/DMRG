from __future__ import unicode_literals
from matplotlib import rc
import matplotlib.pyplot  as plt
import pandas as pd
import os.path
import numpy as np
import seaborn as sns

import itertools

plt.close('all')



####################################################################################
#######           Use these settings by default                            #########
####################################################################################
rc('font', **{'family': 'serif', 'serif': ['Palatino']})
rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})
rc('text', usetex=True)
# sns.set(style="whitegrid", font_scale=1.5 ,rc={'text.usetex' : True})
sns.set(style="darkgrid", font_scale=1.0, font='Helvetica',rc={"lines.linewidth": 1.2})
paper_rc = {'lines.linewidth': 2, 'lines.markersize': 10}
sns.set_style({"axes.facecolor": ".9"},rc={'text.usetex' : True})
sns.set_palette(sns.color_palette("husl", 8))
####################################################################################


####################################################################################
####### Use these settings for presentations or where big fonts are needed #########
####################################################################################
# rc('font', **{'family': 'serif', 'serif': ['Palatino']})
# rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})
# rc('text', usetex=True)
# # sns.set(style="whitegrid", font_scale=1.5 ,rc={'text.usetex' : True})
# sns.set(style="darkgrid", font_scale=1.5, font='Helvetica',rc={"lines.linewidth": 1.2})
# paper_rc = {'lines.linewidth': 2, 'lines.markersize': 16}
# sns.set_style({"axes.facecolor": ".9"},rc={'text.usetex' : True})
# sns.set_palette(sns.color_palette("husl", 5))
####################################################################################

filename = ['../output/state-ghz-1.h5', '../output/state-w-4.h5', '../output/state-rps.h5', '../output/state-rands.h5','../output/state-upup-4.h5', '../output/state-updown.h5']
# filename = ['../output/state-randosm-13.h5', '../output/state-random-product-state.h5']
tag = ['GHZ','W','Random Product State', 'Random State ($\chi$)','$|\\uparrow,\\uparrow\\rangle$ ', '$|\\uparrow,\downarrow\\rangle$ ']
# tag = [ 'Random State ($\chi$)','Random Product State']
markerlist = itertools.cycle(('<', '>', '^','v'))

def plt_graph(table, xkey,ykey, xlabel,ylabel, ax,label='', scale='linear'):
    ax.plot(store[table].get(xkey),store[table].get(ykey), marker=next(markerlist), label=label, alpha=.90, markersize=4)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_yscale(scale)
    ax.legend()

fig, ax = plt.subplots(2, 2, sharex=False)
ax[0, 0].axhline(y=-1.2732395447351625, xmin=0, xmax=1, c="blue", linewidth=0.5, zorder=0,
                 label='$\langle e\\rangle_{exact}, \chi = \infty$')
group = "iDMRG"
plt.suptitle(group)
i = 0;
for name in filename:
    current_tag = tag[i]
    i = i+1
    store = pd.HDFStore(name)
    keys = [s for s in store.keys() if "/"+group in s]
    for key in keys:
        chi_label =  '$\chi=$' + str( store[key].chi_max[0])
        # plt_graph(key,"iteration", "energy1", "L", "E/L",                       ax[0,0], label="$MPO$")
        plt_graph(key,"iteration", "energy2", "L", "E/L",                       ax[0,0], label=current_tag)
        # plt_graph(key,"iteration", "energy3", "L", "E/L",                       ax[0,0], label="$d\log{G_\infty(a)} / d a$")
        # plt_graph(key,"iteration", "energy4", "L", "E/L",                       ax[0,0], label="$d\log{G_\infty(b)} / d b$")
        # plt_graph(key,"iteration", "energy5", "L", "E/L",                       ax[0,0], label="$d\log{G_\infty(c)} / d c$")
        # plt_graph(key,"iteration", "energy6", "L", "E/L",                       ax[0,0], label="$d\log{G_\infty(d)} / d d$")
        plt_graph(key,"iteration", "entropy", "L", "S"  ,                       ax[0,1], label=current_tag)
        # plt_graph(key,"iteration", "truncation_error", "L", "$\sigma^2(E)/L$",  ax[1,0], label="trunc err",                                                                     scale = 'log')
        # plt_graph(key,"iteration", "variance1", "L", "$\sigma^2(E)/L$",         ax[1,0], label="MPO"                         ,                                                  scale = 'log')
        plt_graph(key,"iteration", "variance2", "L", "$\sigma^2(E)/L$",         ax[1,0], label=current_tag,                                                                         scale = 'log')
        # plt_graph(key,"iteration", "variance3", "L", "$\sigma^2(E)/L$",         ax[1,0], label="$log(|G1(a)|^2) / a^2$"      ,                                                  scale = 'log')
        # plt_graph(key,"iteration", "variance4", "L", "$\sigma^2(E)/L$",         ax[1,0], label="$log(|G1(b)|^2) / b^2$"      ,                                                  scale = 'log')
        # plt_graph(key,"iteration", "variance5", "L", "$\sigma^2(E)/L$",         ax[1,0], label="$log(|G1(c)|^2) / c^2$"      ,                                                  scale = 'log')
        # plt_graph(key,"iteration", "variance6", "L", "$\sigma^2(E)/L$",         ax[1,0], label="$log(|G1(d)|^2) / d^2$"      ,                                                  scale = 'log')
        # plt_graph(key,"iteration", "truncation_error", "L", "$\sigma^2(E)/L$",  ax[1,1], label="trunc err",                                                                     scale = 'linear')
        # plt_graph(key,"iteration", "variance1", "L", "$\sigma^2(E)/L$",         ax[1,1], label="MPO"                         ,                                                  scale = 'linear')
        plt_graph(key,"iteration", "variance2", "L", "$\sigma^2(E)/L$",         ax[1,1], label=current_tag,                                                                        scale = 'linear')
        # plt_graph(key,"iteration", "variance3", "L", "$\sigma^2(E)/L$",         ax[1,1], label="$log(|G1(a)|^2) / a^2$"      ,                                                  scale = 'linear')
        # plt_graph(key,"iteration", "variance4", "L", "$\sigma^2(E)/L$",         ax[1,1], label="$log(|G1(b)|^2) / b^2$"      ,                                                  scale = 'linear')
        # plt_graph(key,"iteration", "variance5", "L", "$\sigma^2(E)/L$",         ax[1,1], label="$log(|G1(c)|^2) / c^2$"      ,                                                  scale = 'linear')
        # plt_graph(key,"iteration", "variance6", "L", "$\sigma^2(E)/L$",         ax[1,1], label="$log(|G1(d)|^2) / d^2$"      ,                                                  scale = 'linear')

plt.tight_layout()
plt.subplots_adjust(wspace=.4, hspace=.2)
store.close()
plt.show()
exit()