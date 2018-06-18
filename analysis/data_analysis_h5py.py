from __future__ import unicode_literals
from matplotlib import rc
import matplotlib.pyplot  as plt
import pandas as pd
import os.path
import numpy as np
import seaborn as sns
import h5py
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

filename = '../output/output-835.h5'

if(not os.path.exists(filename)):
    print("File does not exist.")
    exit(1)

def print_dataset_names(name):
    print (name)

def print_dataset_names2(name, obj):
    print (name)

import re

def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    '''
    return [ atoi(c) for c in re.split('(\d+)', text) ]



h5file = h5py.File(filename, 'r')

iDMRG_exists        = "iDMRG" in h5file
fDMRG_exists        = "fDMRG" in h5file
xDMRG_exists        = "xDMRG" in h5file
iTEBD_exists        = "iTEBD" in h5file

markerlist = itertools.cycle(('<', '>', '^','v'))
def plt_graph(table, xkey,ykey, xlabel,ylabel, ax,label='', scale='linear'):
    ax.plot(table[xkey],table[ykey], marker=next(markerlist), label=label, alpha=.90, markersize=4)
    if xlabel != '':  ax.set_xlabel(xlabel)
    if ylabel != '':  ax.set_ylabel(ylabel)
    ax.set_yscale(scale)
    if label != '': ax.legend()

def plt_graph_range(table, xkey,ykey,range, xlabel,ylabel, ax,label='', scale='linear'):
    ax.plot(table[xkey][range],table[ykey][range], marker=next(markerlist), label=label, alpha=.90, markersize=4)
    if xlabel != '':  ax.set_xlabel(xlabel)
    if ylabel != '':  ax.set_ylabel(ylabel)
    ax.set_yscale(scale)
    if label != '': ax.legend()




if(iDMRG_exists):
    print("Printing iDMRG")
    table = h5file.get('iDMRG/iDMRG')
    fig, ax = plt.subplots(2, 2, sharex=False)
    plt.suptitle("iDMRG ($\chi = " + str(table["chi_max"][0]) + "$)" )
    ax[0, 0].axhline(y=-1.2732395447351625, xmin=0, xmax=1, c="blue", linewidth=0.5, zorder=0,label='$\chi,L = \infty$')
    plt_graph(table, "iteration", "energy_mpo",              "Step", "$e$ (per site)", ax[0, 0], label="mpo")
    plt_graph(table, "iteration", "energy_ham",              "Step", "$e$ (per site)", ax[0, 0], label="ham")
    plt_graph(table, "iteration", "energy_mom",              "Step", "$e$ (per site)", ax[0, 0], label="mom")
    plt_graph(table, "iteration", "entanglement_entropy",    "Step", "$S_E$"         , ax[0, 1], label="Entanglement entropy")

    plt_graph(table, "iteration", "variance_mpo",            "Step", "$\log \sigma^2(e)$", ax[1, 0], label="mpo", scale='log')
    plt_graph(table, "iteration", "variance_ham",            "Step", "$\log \sigma^2(e)$", ax[1, 0], label="ham", scale='log')
    plt_graph(table, "iteration", "variance_mom",            "Step", "$\log \sigma^2(e)$", ax[1, 0], label="mom", scale='log')
    plt_graph(table, "iteration", "truncation_error",        "Step", "$\log \sigma^2(e)$", ax[1, 0], label="truncation", scale='log')




if(fDMRG_exists):
    print("Printing fDMRG")
    table = h5file.get('fDMRG/fDMRG')
    fig, ax = plt.subplots(2, 3, sharex=False)
    plt.suptitle("fDMRG ($\chi = " + str(table["chi_max"][0]) + "$)" )
    ax[0, 0].axhline(y=-1.2732395447351625, xmin=0, xmax=1, c="blue", linewidth=0.5, zorder=0,label='$\chi = \infty$')
    plt_graph(table, "iteration", "energy_mpo",              "Step", "$e$ (per site)", ax[0, 0], label="mpo")
    plt_graph(table, "iteration", "entanglement_entropy",    "Step", "$S_E$"         , ax[0, 1], label="Entanglement entropy")
    plt_graph(table, "iteration", "variance_mpo",            "Step", "$\log \sigma^2(e)$", ax[1, 0], label="mpo", scale='log')
    plt_graph(table, "iteration", "truncation_error",        "Step", "$\log \sigma^2(e)$", ax[1, 0], label="truncation", scale='log')
    table = h5file.get('fDMRG/fDMRG_chain')
    imax = table["iteration"][-1]
    for i in range(0,imax):
        plt_graph_range(table, "position", "energy"              , table["iteration"] == i,  "Site", "$e$ (per site)"      , ax[1, 1])
        plt_graph_range(table, "position", "entanglement_entropy", table["iteration"] == i,  "Site", "$S_E$"               , ax[0, 2])
        plt_graph_range(table, "position", "truncation_error"    , table["iteration"] == i,  "Site", "$\epsilon$"          , ax[1, 2], scale='log')






# h5file.
if(xDMRG_exists):
    print("Printing xDMRG")
    table = h5file.get('xDMRG/xDMRG')
    fig, ax = plt.subplots(2, 4, sharex=False)
    plt.suptitle("xDMRG ($\chi = " + str(table["chi_max"][0]) + "$)" )
    plt_graph(table, "iteration", "energy_mpo",              "Step", "$e$ (per site)", ax[0, 0], label="mpo")
    plt_graph(table, "iteration", "entanglement_entropy",    "Step", "$S_E$"         , ax[0, 1], label="Entanglement entropy")
    plt_graph(table, "iteration", "variance_mpo",            "Step", "$\log \sigma^2(e)$", ax[1, 0], label="mpo", scale='log')
    plt_graph(table, "iteration", "truncation_error",        "Step", "$\log \sigma^2(e)$", ax[1, 0], label="truncation", scale='log')
    table = h5file.get('xDMRG/xDMRG_chain')
    imax = table["iteration"][-1]
    for i in range(imax-3,imax):
        plt_graph_range(table, "position", "energy"              , table["iteration"] == i,  "Site", "$e$ (per site)"      , ax[1, 1])
        plt_graph_range(table, "position", "entanglement_entropy", table["iteration"] == i,  "Site", "$S_E$"               , ax[0, 2])
        plt_graph_range(table, "position", "truncation_error"    , table["iteration"] == i,  "Site", "$\epsilon$"          , ax[1, 2], scale='log')
    mpos = [mpo for mpo in h5file["xDMRG/chain/MPO"].keys()]
    mpos.sort(key=natural_keys)
    randfield=[h5file["xDMRG/chain/MPO/"+mpo].attrs["random_field"][0] for mpo in mpos]
    ax[0,3].plot(range(len(randfield)),randfield, marker=next(markerlist), alpha=.90, markersize=4)
    ax.set_xlabel('Site')
    ax.set_ylabel('MPO Random field')


plt.tight_layout()
plt.subplots_adjust(wspace=.4, hspace=.2)
h5file.close()
plt.show()
exit()