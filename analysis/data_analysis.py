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

filename = '../output/output-576.h5'

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
    group = "iDMRG/iDMRG"
    keys = [s for s in store.keys() if "/"+group in s]
    fig, ax = plt.subplots(2, 2, sharex=False)
    plt.suptitle(group)
    ax[0,0].axhline(y=-1.2732395447351625, xmin=0, xmax=1, c="blue", linewidth=0.5, zorder=0,label='$\langle e\\rangle_{exact}, \chi = \infty$')
    for key in keys:
        chi_label =  '$\chi=$' + str( store[key].chi_max[0])
        plt_graph(key,"iteration", "energy_mpo", "L", "E/L",                    ax[0,0], label="$MPO$")
        plt_graph(key,"iteration", "energy_ham", "L", "E/L",                    ax[0,0], label="$\langle h_{even} + h_{odd}  \\rangle$")
        plt_graph(key,"iteration", "energy_mom", "L", "E/L",                    ax[0,0], label="$d\log{G_\infty(a)} / d a$")
        plt_graph(key,"iteration", "entanglement_entropy",    "L", "S"  ,                    ax[0,1], label=chi_label)
        plt_graph(key,"iteration", "truncation_error", "L", "$\sigma^2(E)/L$",  ax[1,0], label="trunc err",                        scale = 'log')
        plt_graph(key,"iteration", "variance_mpo", "L", "$\sigma^2(E)/L$",      ax[1,0], label="MPO"                         ,     scale = 'log')
        plt_graph(key,"iteration", "variance_ham", "L", "$\sigma^2(E)/L$",      ax[1,0], label="$\langle (H-E)^2 \\rangle_1$",     scale = 'log')
        plt_graph(key,"iteration", "variance_mom", "L", "$\sigma^2(E)/L$",      ax[1,0], label="$log(|G(a)|^2) / a^2$"      ,     scale = 'log')
        plt_graph(key,"iteration", "truncation_error", "L", "$\sigma^2(E)/L$",  ax[1,1], label="trunc err",                        scale = 'linear')
        plt_graph(key,"iteration", "variance_mpo", "L", "$\sigma^2(E)/L$",      ax[1,1], label="MPO"                         ,     scale = 'linear')
        plt_graph(key,"iteration", "variance_ham", "L", "$\sigma^2(E)/L$",      ax[1,1], label="$\langle (H-E)^2 \\rangle_1$",     scale = 'linear')
        plt_graph(key,"iteration", "variance_mom", "L", "$\sigma^2(E)/L$",      ax[1,1], label="$log(|G(a)|^2) / a^2$"      ,     scale = 'linear')

if fDMRG_exists:
    group = "fDMRG/fDMRG"
    keys = [s for s in store.keys() if "/"+group in s]
    exit(0)

    fig, ax = plt.subplots(2, 2, sharex=False)
    plt.suptitle(group)
    ax[0, 0].axhline(y=-1.2732395447351625, xmin=0, xmax=1, c="blue", linewidth=0.5, zorder=0,label='$\chi = \infty$')
    for key in keys:
        chi_label =  '$\chi=$' + str( store[key].chi_max[0])
        plt_graph(key, "iteration", "energy_mpo",          "L", "E/L", ax[0, 0], label="MPO")
        plt_graph(key, "iteration", "energy_ham",          "L", "E/L", ax[0, 0], label="$\langle H \\rangle$")
        plt_graph(key, "iteration", "energy_mom",          "L", "E/L", ax[0, 0], label="$d\log{G(a)} / da$")
        plt_graph(key, "iteration", "entanglement_entropy","L", "S",   ax[0,1], label=chi_label)
        plt_graph(key, "iteration", "truncation_error", "L", "$\sigma^2(E)/L$", ax[1,0], label="trunc err", scale = 'log')
        plt_graph(key, "iteration", "variance_mpo", "L", "$\sigma^2(E)/L$", ax[1,0], label="MPO"                         , scale = 'log')
        plt_graph(key, "iteration", "variance_ham", "L", "$\sigma^2(E)/L$", ax[1,0], label="$\langle (H-E)^2 \\rangle_1$", scale = 'log')
        plt_graph(key, "iteration", "variance_mom", "L", "$\sigma^2(E)/L$", ax[1,0], label="$log(|G(a)|^2) / a^2$"      , scale = 'log')
        plt_graph(key, "iteration", "truncation_error", "L", "$\sigma^2(E)/L$", ax[1,1], label="trunc err", scale = 'linear')
        plt_graph(key, "iteration", "variance_mpo", "L", "$\sigma^2(E)/L$", ax[1,1], label="MPO"                         , scale = 'linear')
        plt_graph(key, "iteration", "variance_ham", "L", "$\sigma^2(E)/L$", ax[1,1], label="$\langle (H-E)^2 \\rangle_1$", scale = 'linear')
        plt_graph(key, "iteration", "variance_mom", "L", "$\sigma^2(E)/L$", ax[1,1], label="$log(|G(a)|^2) / a^2$"      , scale = 'linear')



if xDMRG_exists:
    group = "xDMRG/xDMRG"
    keys = [s for s in store.keys() if "/"+group in s]
    fig, ax = plt.subplots(2, 2, sharex=False)
    plt.suptitle(group)
    # ax[0, 0].axhline(y=-1.2732395447351625, xmin=0, xmax=1, c="blue", linewidth=0.5, zorder=0,label='$\chi = \infty$')
    for key in keys:
        chi_label =  '$\chi=$' + str( store[key].chi_max[0])
        plt_graph(key,"iteration", "energy_mpo", "L", "E1", ax[0,0], label=chi_label)
        plt_graph(key,"iteration", "energy_ham", "L", "E2", ax[0,0], label=chi_label)
        plt_graph(key,"iteration", "energy_mom", "L", "E3", ax[0,0], label=chi_label)
        plt_graph(key,"iteration", "entanglement_entropy","L", "S", ax[0,1], label=chi_label)
        plt_graph(key,"chain_length", "variance_mpo", "L", "$\sigma^2(E)/L$", ax[1,0], label="variance", scale = 'log')
        plt_graph(key,"chain_length", "variance_ham", "L", "$\sigma^2(E)/L$", ax[1,0], label="variance_ham", scale = 'log')
        plt_graph(key,"chain_length", "variance_mom", "L", "$\sigma^2(E)/L$", ax[1,0], label="variance_mom", scale = 'log')


if iTEBD_exists:
    group = "iTEBD/iTEBD"
    keys = [s for s in store.keys() if "/"+group in s]
    fig, ax = plt.subplots(2, 2, sharex=False)
    plt.suptitle(group)
    ax[0, 0].axhline(y=-1.2732395447351625, xmin=0, xmax=1, c="blue", linewidth=0.5, zorder=0,label='$\chi = \infty$')
    for key in keys:
        chi_label =  '$\chi=$' + str( store[key].chi_max[0])
        plt_graph(key,"iteration", "energy_ham", "L", "E/L", ax[0, 0], label="$\langle H \\rangle$")
        plt_graph(key,"iteration", "energy_mom", "L", "E/L", ax[0, 0], label="$d \log{G(a)} / da$")
        plt_graph(key,"iteration", "entanglement_entropy","step", "S", ax[0,1], label=chi_label)
        plt_graph(key,"iteration", "truncation_error", "L", "$\sigma^2(E)/L$", ax[1,0], label="trunc err", scale = 'log')
        plt_graph(key,"iteration", "variance_ham", "step", "$\sigma^2(E)/L$", ax[1,0], label="$\langle (H-E)^2 \\rangle_1$", scale = 'log')
        plt_graph(key,"iteration", "variance_mom", "step", "$\sigma^2(E)/L$", ax[1,0], label="$log(|G(a)|^2) / a^2$"      , scale = 'log')
        plt_graph(key,"iteration", "truncation_error", "L", "$\sigma^2(E)/L$", ax[1,1], label="trunc err", scale = 'linear')
        plt_graph(key,"iteration", "variance_ham", "step", "$\sigma^2(E)/L$", ax[1,1], label="$\langle (H-E)^2 \\rangle_1$", scale = 'linear')
        plt_graph(key,"iteration", "variance_mom", "step", "$\sigma^2(E)/L$", ax[1,1], label="$log(|G(a)|^2) / a^2$"      , scale = 'linear')


plt.tight_layout()
plt.subplots_adjust(wspace=.4, hspace=.2)
store.close()
plt.show()
exit()
