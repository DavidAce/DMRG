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

filename = '../output/output-445.h5'

if(not os.path.exists(filename)):
    print("File does not exist.")
    exit(1)

store               = pd.HDFStore(filename)
iDMRG_exists        = "iDMRG" in store
fDMRG_exists        = "fDMRG" in store
xDMRG_exists        = "xDMRG" in store
iTEBD_exists        = "iTEBD" in store


markerlist = itertools.cycle(('<', '>', '^','v'))

def plt_graph(table, xkey,ykey, xlabel,ylabel, ax,label='', scale='linear'):
    ax.plot(store[table].get(xkey),store[table].get(ykey), marker=next(markerlist), label=label, alpha=.90, markersize=4)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_yscale(scale)
    ax.legend()


if iDMRG_exists:
    group = "iDMRG"
    keys = [s for s in store.keys() if "/"+group in s]
    fig, ax = plt.subplots(2, 2, sharex=False)
    plt.suptitle(group)
    ax[0,0].axhline(y=-1.2732395447351625, xmin=0, xmax=1, c="blue", linewidth=0.5, zorder=0,label='$\langle e\\rangle_{exact}, \chi = \infty$')
    for key in keys:
        chi_label =  '$\chi=$' + str( store[key].chi_max[0])
        plt_graph(key,"iteration", "energy1", "L", "E/L",                       ax[0,0], label="$MPO$")
        plt_graph(key,"iteration", "energy2", "L", "E/L",                       ax[0,0], label="$\langle h_{even} + h_{odd}  \\rangle$")
        plt_graph(key,"iteration", "energy3", "L", "E/L",                       ax[0,0], label="$d\log{G_\infty(a)} / d a$")
        plt_graph(key,"iteration", "energy4", "L", "E/L",                       ax[0,0], label="$d\log{G_\infty(b)} / d b$")
        plt_graph(key,"iteration", "energy5", "L", "E/L",                       ax[0,0], label="$d\log{G_\infty(c)} / d c$")
        plt_graph(key,"iteration", "energy6", "L", "E/L",                       ax[0,0], label="$d\log{G_\infty(d)} / d d$")
        plt_graph(key,"iteration", "entropy", "L", "S"  ,                       ax[0,1], label=chi_label)
        plt_graph(key,"iteration", "truncation_error", "L", "$\sigma^2(E)/L$",  ax[1,0], label="trunc err",                                                                     scale = 'log')
        plt_graph(key,"iteration", "variance1", "L", "$\sigma^2(E)/L$",         ax[1,0], label="MPO"                         ,                                                  scale = 'log')
        plt_graph(key,"iteration", "variance2", "L", "$\sigma^2(E)/L$",         ax[1,0], label="$\langle (H-E)^2 \\rangle_1$",                                                  scale = 'log')
        plt_graph(key,"iteration", "variance3", "L", "$\sigma^2(E)/L$",         ax[1,0], label="$log(|G1(a)|^2) / a^2$"      ,                                                  scale = 'log')
        plt_graph(key,"iteration", "variance4", "L", "$\sigma^2(E)/L$",         ax[1,0], label="$log(|G1(b)|^2) / b^2$"      ,                                                  scale = 'log')
        plt_graph(key,"iteration", "variance5", "L", "$\sigma^2(E)/L$",         ax[1,0], label="$log(|G1(c)|^2) / c^2$"      ,                                                  scale = 'log')
        plt_graph(key,"iteration", "variance6", "L", "$\sigma^2(E)/L$",         ax[1,0], label="$log(|G1(d)|^2) / d^2$"      ,                                                  scale = 'log')
        plt_graph(key,"iteration", "truncation_error", "L", "$\sigma^2(E)/L$",  ax[1,1], label="trunc err",                                                                     scale = 'linear')
        plt_graph(key,"iteration", "variance1", "L", "$\sigma^2(E)/L$",         ax[1,1], label="MPO"                         ,                                                  scale = 'linear')
        plt_graph(key,"iteration", "variance2", "L", "$\sigma^2(E)/L$",         ax[1,1], label="$\langle (H-E)^2 \\rangle_1$",                                                  scale = 'linear')
        plt_graph(key,"iteration", "variance3", "L", "$\sigma^2(E)/L$",         ax[1,1], label="$log(|G1(a)|^2) / a^2$"      ,                                                  scale = 'linear')
        plt_graph(key,"iteration", "variance4", "L", "$\sigma^2(E)/L$",         ax[1,1], label="$log(|G1(b)|^2) / b^2$"      ,                                                  scale = 'linear')
        plt_graph(key,"iteration", "variance5", "L", "$\sigma^2(E)/L$",         ax[1,1], label="$log(|G1(c)|^2) / c^2$"      ,                                                  scale = 'linear')
        plt_graph(key,"iteration", "variance6", "L", "$\sigma^2(E)/L$",         ax[1,1], label="$log(|G1(d)|^2) / d^2$"      ,                                                  scale = 'linear')

if fDMRG_exists:
    group = "fDMRG"
    keys = [s for s in store.keys() if "/"+group in s]
    fig, ax = plt.subplots(2, 2, sharex=False)
    plt.suptitle(group)
    ax[0, 0].axhline(y=-1.2732395447351625, xmin=0, xmax=1, c="blue", linewidth=0.5, zorder=0,label='$\chi = \infty$')
    for key in keys:
        chi_label =  '$\chi=$' + str( store[key].chi_max[0])
        plt_graph(key, "iteration", "energy1", "L", "E/L", ax[0, 0], label="MPO")
        plt_graph(key, "iteration", "energy2", "L", "E/L", ax[0, 0], label="$\langle H \\rangle$")
        plt_graph(key, "iteration", "energy3", "L", "E/L", ax[0, 0], label="$d\log{G1(a)} / da$")
        plt_graph(key, "iteration", "energy4", "L", "E/L", ax[0, 0], label="$d\log{G1(b)} / db$")
        plt_graph(key, "iteration", "energy5", "L", "E/L", ax[0, 0], label="$d\log{G1(c)} / dc$")
        plt_graph(key, "iteration", "energy6", "L", "E/L", ax[0, 0], label="$d\log{G1(d)} / dd$")
        plt_graph(key, "iteration", "entropy","L", "S", ax[0,1], label=chi_label)
        plt_graph(key, "iteration", "truncation_error", "L", "$\sigma^2(E)/L$", ax[1,0], label="trunc err", scale = 'log')
        plt_graph(key, "iteration", "variance1", "L", "$\sigma^2(E)/L$", ax[1,0], label="MPO"                         , scale = 'log')
        plt_graph(key, "iteration", "variance2", "L", "$\sigma^2(E)/L$", ax[1,0], label="$\langle (H-E)^2 \\rangle_1$", scale = 'log')
        plt_graph(key, "iteration", "variance3", "L", "$\sigma^2(E)/L$", ax[1,0], label="$log(|G1(a)|^2) / a^2$"      , scale = 'log')
        plt_graph(key, "iteration", "variance4", "L", "$\sigma^2(E)/L$", ax[1,0], label="$log(|G1(b)|^2) / b^2$"      , scale = 'log')
        plt_graph(key, "iteration", "variance5", "L", "$\sigma^2(E)/L$", ax[1,0], label="$log(|G1(c)|^2) / c^2$"      , scale = 'log')
        plt_graph(key, "iteration", "variance6", "L", "$\sigma^2(E)/L$", ax[1,0], label="$log(|G1(d)|^2) / d^2$"      , scale = 'log')
        plt_graph(key, "iteration", "truncation_error", "L", "$\sigma^2(E)/L$", ax[1,1], label="trunc err", scale = 'linear')
        plt_graph(key, "iteration", "variance1", "L", "$\sigma^2(E)/L$", ax[1,1], label="MPO"                         , scale = 'linear')
        plt_graph(key, "iteration", "variance2", "L", "$\sigma^2(E)/L$", ax[1,1], label="$\langle (H-E)^2 \\rangle_1$", scale = 'linear')
        plt_graph(key, "iteration", "variance3", "L", "$\sigma^2(E)/L$", ax[1,1], label="$log(|G1(a)|^2) / a^2$"      , scale = 'linear')
        plt_graph(key, "iteration", "variance4", "L", "$\sigma^2(E)/L$", ax[1,1], label="$log(|G1(b)|^2) / b^2$"      , scale = 'linear')
        plt_graph(key, "iteration", "variance5", "L", "$\sigma^2(E)/L$", ax[1,1], label="$log(|G1(c)|^2) / c^2$"      , scale = 'linear')
        plt_graph(key, "iteration", "variance6", "L", "$\sigma^2(E)/L$", ax[1,1], label="$log(|G1(d)|^2) / d^2$"      , scale = 'linear')



if xDMRG_exists:
    group = "xDMRG"
    keys = [s for s in store.keys() if "/"+group in s]
    fig, ax = plt.subplots(2, 2, sharex=False)
    plt.suptitle(group)
    # ax[0, 0].axhline(y=-1.2732395447351625, xmin=0, xmax=1, c="blue", linewidth=0.5, zorder=0,label='$\chi = \infty$')
    for key in keys:
        chi_label =  '$\chi=$' + str( store[key].chi_max[0])
        plt_graph(key,"iteration", "energy1", "L", "E1", ax[0,0], label=chi_label)
        plt_graph(key,"iteration", "energy2", "L", "E2", ax[0,0], label=chi_label)
        plt_graph(key,"iteration", "energy3", "L", "E3", ax[0,0], label=chi_label)
        plt_graph(key,"iteration", "energy4", "L", "E4", ax[0,0], label=chi_label)
        plt_graph(key,"iteration", "energy5", "L", "E5", ax[0,0], label=chi_label)
        plt_graph(key,"iteration", "energy6", "L", "E6", ax[0,0], label=chi_label)
        plt_graph(key,"position", "entropy","L", "S", ax[0,1], label=chi_label)
        plt_graph(key,"chain_length", "variance1", "L", "$\sigma^2(E)/L$", ax[1,0], label="variance", scale = 'log')
        plt_graph(key,"chain_length", "variance2", "L", "$\sigma^2(E)/L$", ax[1,0], label="variance2", scale = 'log')
        plt_graph(key,"chain_length", "variance3", "L", "$\sigma^2(E)/L$", ax[1,0], label="variance3", scale = 'log')


if iTEBD_exists:
    group = "iTEBD"
    keys = [s for s in store.keys() if "/"+group in s]
    fig, ax = plt.subplots(2, 2, sharex=False)
    plt.suptitle(group)
    ax[0, 0].axhline(y=-1.2732395447351625, xmin=0, xmax=1, c="blue", linewidth=0.5, zorder=0,label='$\chi = \infty$')
    for key in keys:
        chi_label =  '$\chi=$' + str( store[key].chi_max[0])
        plt_graph(key, "iteration", "energy1", "L", "E/L", ax[0, 0], label="MPO")
        plt_graph(key, "iteration", "energy2", "L", "E/L", ax[0, 0], label="$\langle H \\rangle$")
        plt_graph(key, "iteration", "energy3", "L", "E/L", ax[0, 0], label="$d \log{G1(a)} / da$")
        plt_graph(key, "iteration", "energy4", "L", "E/L", ax[0, 0], label="$d \log{G1(b)} / db$")
        plt_graph(key, "iteration", "energy5", "L", "E/L", ax[0, 0], label="$d \log{G1(c)} / dc$")
        plt_graph(key, "iteration", "energy6", "L", "E/L", ax[0, 0], label="$d \log{G1(d)} / dd$")
        plt_graph(key,"iteration", "entropy","step", "S", ax[0,1], label=chi_label)
        plt_graph(key,"iteration", "truncation_error", "L", "$\sigma^2(E)/L$", ax[1,0], label="trunc err", scale = 'log')
        plt_graph(key,"iteration", "variance1", "step", "$\sigma^2(E)/L$", ax[1,0], label="MPO"                         , scale = 'log')
        plt_graph(key,"iteration", "variance2", "step", "$\sigma^2(E)/L$", ax[1,0], label="$\langle (H-E)^2 \\rangle_1$", scale = 'log')
        plt_graph(key,"iteration", "variance3", "step", "$\sigma^2(E)/L$", ax[1,0], label="$log(|G1(a)|^2) / a^2$"      , scale = 'log')
        plt_graph(key,"iteration", "variance4", "step", "$\sigma^2(E)/L$", ax[1,0], label="$log(|G1(b)|^2) / b^2$"      , scale = 'log')
        plt_graph(key,"iteration", "variance5", "step", "$\sigma^2(E)/L$", ax[1,0], label="$log(|G1(c)|^2) / c^2$"      , scale = 'log')
        plt_graph(key,"iteration", "variance6", "step", "$\sigma^2(E)/L$", ax[1,0], label="$log(|G1(d)|^2) / d^2$"      , scale = 'log')
        plt_graph(key,"iteration", "truncation_error", "L", "$\sigma^2(E)/L$", ax[1,1], label="trunc err", scale = 'linear')
        plt_graph(key,"iteration", "variance1", "step", "$\sigma^2(E)/L$", ax[1,1], label="MPO"                         , scale = 'linear')
        plt_graph(key,"iteration", "variance2", "step", "$\sigma^2(E)/L$", ax[1,1], label="$\langle (H-E)^2 \\rangle_1$", scale = 'linear')
        plt_graph(key,"iteration", "variance3", "step", "$\sigma^2(E)/L$", ax[1,1], label="$log(|G1(a)|^2) / a^2$"      , scale = 'linear')
        plt_graph(key,"iteration", "variance4", "step", "$\sigma^2(E)/L$", ax[1,1], label="$log(|G1(b)|^2) / b^2$"      , scale = 'linear')
        plt_graph(key,"iteration", "variance5", "step", "$\sigma^2(E)/L$", ax[1,1], label="$log(|G1(c)|^2) / c^2$"      , scale = 'linear')
        plt_graph(key,"iteration", "variance6", "step", "$\sigma^2(E)/L$", ax[1,1], label="$log(|G1(d)|^2) / d^2$"      , scale = 'linear')


plt.tight_layout()
plt.subplots_adjust(wspace=.4, hspace=.2)
store.close()
plt.show()
exit()
