import os
import sys
here = os.path.dirname(__file__)
sys.path.insert(0, os.path.join(here, '..'))
sys.path.insert(0, os.path.join(here, '../..'))
sys.path.insert(0, os.path.join(here, ''))
import numpy as np
import h5py
import matplotlib.pyplot as plt
from pathlib import Path
from matplotlib import transforms
import matplotlib.patheffects as pe
from plotting.tools import *
import seaborn as sns



def legend_without_duplicate_labels(ax, loc=None):
    handles, labels = ax.get_legend_handles_labels()
    unique = [(h, l) for i, (h, l) in enumerate(zip(handles, labels)) if l not in labels[:i]]
    print(unique)
    ax.legend(*zip(*unique), loc=loc)
def xover(ts, cut):
    x = ts > cut
    return x.argmax() if x.any() else -1

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "sans-serif",
    "font.sans-serif": "STIX",
})
fig1 = plt.figure(figsize=(12,8))
fig1.suptitle('Half-chain number probability $p$ when quenching the Néel state')
pe_fields = [pe.Stroke(linewidth=2, foreground='white'), pe.Normal()]

seed =29
# seed =16
d = 0
fieldinset=True
analysisfile = "/mnt/WDB-AN1500/mbl_transition/lbit93-precision/analysis/data/averaged.h5"
simpath = "L24/J[+0±1_+0±1_+0±1]/x1/rL/u[d16_f0.4_tw1.0ID_cw1.0EX_bond0]/fLBIT"

with h5py.File(analysisfile) as h5file:
    simnode = h5file[f'{simpath}']
    L = simnode['model/model_size'][()]
    hamiltonian = simnode['model/hamiltonian']
    p = h5file[f'{simpath}/state_real/number_probabilities'][()]
    pdvg = np.mean(p, axis=2)
    times = simnode['state_real/cronos/measurements/avg']['physical_time'].astype(float)
    with np.errstate(invalid='ignore'):
        tlnln = np.log(np.log(times))

    db = {'vals': {}}
    db['vals']['L'] = L
    db['vals']['r'] = hamiltonian['J2_span'][()]
    db['vals']['x'] = hamiltonian['xi_Jcls'][()]
    db['vals']['w'] = (hamiltonian['J1_wdth'][()], hamiltonian['J2_wdth'][()], hamiltonian['J3_wdth'][()])
    t = get_timepoints(times, db)

    pidxs = [int(L // 4 + 0), int(L // 4 - 1), int(L // 4 + 1), int(L // 4 - 2), int(L // 4 + 2)]
    colors = ['Blues', 'Oranges', 'Greens', 'Purples', 'Reds']
    labels = ['$p(L/4)$','$p(L/4-1)$','$p(L/4+1)$','$p(L/4-2)$','$p(L/4+2)$']
    # lbit2 = h5file['fLBIT/model/lbits/corrmat'][0,site2 ,:]
    # sites = [ L // 2 - d, L // 2 + d]
    sites = range(0,L)#[ L // 2 - d, L // 2 + d]


    ax = fig1.add_subplot(2, 2, 1)
    ax.set_yscale('log')
    ax.set_xlabel('$p$')
    ax.set_title("Distribution of $p$\'s during $S_N$ growth\n Light $\\rightarrow$ dark color increases time")
    time_indices = range(t.idx_num_lnlnt_begin,t.idx_num_saturated,6)
    for d,palette,label in zip([0, -1, 1, 2,-2], ['Blues', 'Oranges', 'Greens', 'Purples', 'Reds'],['$p(L/4)$','$p(L/4-1)$','$p(L/4+1)$','$p(L/4-2)$','$p(L/4+2)$']):
        ph, pc = [],[]
        ax.set_prop_cycle('color', sns.color_palette(palette=palette, n_colors=len(time_indices)))
        for tidx in time_indices:
            hist, edge = np.histogram(p[int(L//4) + d, tidx, :], bins=100, range=[0,1], density=True)
            binc = [(edge[j] + edge[j + 1]) / 2. for j in range(len(edge) - 1)]
            ph = np.asarray(ph)
            pc = np.asarray(pc)
            ax.step(x=binc, y=hist, where='mid', label=label if tidx == time_indices[-1] else None, linewidth=1.0)
    ax.legend()

    ax = fig1.add_subplot(2, 2, 2)
    # ax1.set_xscale('log')
    # ax1.set_yscale('log')
    # ax.set_ylim(ymin=1e-8, ymax=1e1)
    # ax.set_ylabel('$|O(x)|$')
    ax.set_ylabel('$\overline p$')
    ax.set_xlabel('$\ln \ln t$')
    ax.set_title('Disorder-averaged time dependence of the largest $p$')
    # ax.set_title(f'$\ell$-bits')
    for pidx,label,color in zip(pidxs[:1], labels, colors):
        ax.set_prop_cycle('color', sns.color_palette(palette=color, n_colors=1))
        ax.plot(tlnln, pdvg[pidx], label=label)
        ax.axvline(tlnln[t.idx_num_lnlnt_begin],color='gray')
        ax.axvline(tlnln[t.idx_num_lnlnt_cease],color='gray')
        # ax.plot(tlnln, p[pidx,:, ::100], label=None, color='black', alpha=0.05)

    ax.legend()

    ax = fig1.add_subplot(2, 2, 3)
    ax.set_xlabel('$\ln \ln t$')
    ax.set_ylabel('$|\overline{p^\infty} - \overline p|$')
    ax.set_yscale('log')
    ax.set_title(f'Approach to saturation (disorder average $p$)')
    for pidx,label,color in zip(pidxs, labels, colors):
        ax.set_prop_cycle('color', sns.color_palette(palette=color, n_colors=1))
        pdvg_inf = np.mean(pdvg[pidx,t.idx_num_saturated:])
        ax.plot(tlnln, np.abs(pdvg_inf-pdvg[pidx]), label=label)
        ax.axvline(tlnln[t.idx_num_lnlnt_begin],color='gray')
        ax.axvline(tlnln[t.idx_num_lnlnt_cease],color='gray')
    ax.legend()

    ax = fig1.add_subplot(2, 2, 4)
    ax.set_xlabel('$\ln \ln t$')
    ax.set_ylabel('$\overline{|p^\infty  - p|}$')
    ax.set_yscale('log')
    ax.set_title(f'Approach to saturation (disorder average $|p^\infty - p|$)')
    for pidx,label,color in zip(pidxs, labels, colors):
        ax.set_prop_cycle('color', sns.color_palette(palette=color, n_colors=1))
        p_inf = np.mean(p[pidx,t.idx_num_saturated:, :], axis=0)
        p_approach = np.mean(np.abs(p_inf - p[pidx,:,:]), axis=1)
        ax.plot(tlnln, p_approach, label=label)
        ax.axvline(tlnln[t.idx_num_lnlnt_begin],color='gray')
        ax.axvline(tlnln[t.idx_num_lnlnt_cease],color='gray')
    ax.legend()




plt.tight_layout()
plt.show()