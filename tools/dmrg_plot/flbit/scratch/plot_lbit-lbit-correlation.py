import numpy as np
import h5py
import matplotlib.pyplot as plt
from pathlib import Path
from matplotlib import transforms
import matplotlib.patheffects as pe
from plotting.tools import *


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
pe_fields = [pe.Stroke(linewidth=2, foreground='white'), pe.Normal()]

seed =29
# seed =16
d = 0
fieldinset=True
outdir = "/home/david/GitProjects/DMRG++/output"
lbit_maxes = []
try:
    ax = fig1.add_subplot(2, 2, 1)
    with h5py.File(f'{outdir}/mbl_{seed}_{d}.h5') as h5file:
        L = int(h5file['fLBIT/model/hamiltonian'].attrs.get('model_size'))
        lbits = h5file['fLBIT/model/lbits/corravg'][()]
        # lbit2 = h5file['fLBIT/model/lbits/corrmat'][0,site2 ,:]
        # sites = [ L // 2 - d, L // 2 + d]
        sites = range(0,L)#[ L // 2 - d, L // 2 + d]
        colors = ['green','orange']
        ax.plot(range(L), lbits[:].T)

        # ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_ylim(ymin=1e-8, ymax=1e1)
        ax.set_ylabel('$|O(x)|$')
        ax.set_xlabel('$l$')
        ax.set_title(f'$\ell$-bits')

except Exception as ex:
    print(ex)
    pass


try:
    ax = fig1.add_subplot(2, 2, 2)
    with h5py.File(f'{outdir}/mbl_{seed}_{d}.h5') as h5file:
        L = int(h5file['fLBIT/model/hamiltonian'].attrs.get('model_size'))
        datanode= h5file['fLBIT/model/lbits/corrmat']
        lbits = np.mean(np.abs(h5file['fLBIT/model/lbits/corrmat']),axis=0)
        fields = h5file["fLBIT/model/hamiltonian"]['J1_rand'].astype('float')
        # lbit2 = h5file['fLBIT/model/lbits/corrmat'][0,site2 ,:]
        # sites = [ L // 2 - d, L // 2 + d]
        sites = range(0,L)#[ L // 2 - d, L // 2 + d]
        colors = ['green','orange']
        for site in sites:
            ax.plot(range(L), np.abs(lbits[site ,:]),
                    # color='blue',
                    colors[np.mod(site,2)],
                    label='$\\tau_{' f'{site}' '}^z$'
                    )

        # ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_ylim(ymin=1e-8, ymax=1e1)
        ax.set_ylabel('$|O(x)|$')
        ax.set_xlabel('$l$')
        ax.set_title(f'$\ell$-bits')
        # ax.plot(range(L),np.abs(lbit1), color='blue', label='$\\tau_{' f'{site1}' '}^z$')
        # ax.plot(range(L),np.abs(lbit2), color='red',  label='$\\tau_{' f'{site2}' '}^z$')
        # ax.legend()
        if fieldinset:
            trans = transforms.blended_transform_factory(ax.transData, ax.transAxes)
            fieldcolor = 'steelblue'
            for site in range(L):
                ax.arrow(site, 0.12, 0.0, fields[site] / 15.0, width=0.025, length_includes_head=False,
                         head_width=0.00, head_length=0.00,
                         path_effects=pe_fields,
                         zorder=50,
                         transform=trans,
                         color=fieldcolor,
                         # color='salmon'
                         )
            ax.arrow(0, 0.12, L, 0, width=0.001, length_includes_head=False,
                     head_width=0.00, head_length=0.00,
                     zorder=60, transform=trans,
                     color=fieldcolor
                     # color='salmon'
                     )
        lbavg = get_lbit_avg(corrmat=datanode[()], site=list(sites), mean='arithmetic')
        yfits = []
        for yfull, ystdv in zip(lbavg.full.T, lbavg.stdv.T):
            yfits.append(get_lbit_fit_data(x=sites, y=np.atleast_2d(yfull).T, e=np.atleast_2d(ystdv).T,
                                           beta=True,
                                           ymin=1e-16,
                                           ))  # Fit to get characteristic length-scale
        xis, cs = [],[]
        for ifit, fit in enumerate(yfits):
            # yfitmask = np.ma.masked_invalid(np.log10(np.abs(fit.yfit)))
            # else:
            yfitmask = np.ma.masked_invalid((np.abs(fit.yfit)))
            ax.plot(sites, yfitmask, linewidth=0.4, marker=None,
                    linestyle='dashed', color=colors[np.mod(ifit,2)])

            print(fit.xi)
            xis.append(fit.xi)
            cs.append(fit.C)
        xi = np.mean(xis)
        c = np.mean(cs)
        for l in [0,1,2,3,4,5]:
            ax.axhline(c*np.exp(-l/2.0/xi), linestyle='dashed', linewidth=1)
        lbit_max = np.mean(np.max(lbits, axis=1))
        lbit_maxes.append(lbit_max)


except Exception as ex:
    print(ex)
    pass

seed = 29
try:
    ax = fig1.add_subplot(2, 2, 3)
    for d in [1, 2, 3, 4, 5]:
        with h5py.File(f'{outdir}/mbl_{seed}_{d}.h5') as h5file:
            L = int(h5file['fLBIT/model/hamiltonian'].attrs.get('model_size'))
            lbits = np.mean(np.abs(h5file['fLBIT/model/lbits/corrmat']), axis=0)
            fields = h5file["fLBIT/model/hamiltonian"]['J1_rand'].astype('float')
            # lbit2 = h5file['fLBIT/model/lbits/corrmat'][0,site2 ,:]
            # sites = [ L // 2 - d, L // 2 + d]
            sites = range(0,L)#[ L // 2 - d, L // 2 + d]
            colors = ['green','orange']
            for site in sites:
                if site + d >= L: continue
                ax.plot(range(L), np.abs(lbits[site ,:]),
                        # color='blue',
                        colors[np.mod(site,2)],
                        label='$\\tau_{' f'{site}' '}^z$'
                        )

            # ax.set_xscale('log')
            ax.set_yscale('log')
            ax.set_ylim(ymin=1e-8, ymax=1e1)
            ax.set_ylabel('$|O(x)|$')
            ax.set_xlabel('$l$')
            ax.set_title(f'$\ell$-bits')
            # ax.plot(range(L),np.abs(lbit1), color='blue', label='$\\tau_{' f'{site1}' '}^z$')
            # ax.plot(range(L),np.abs(lbit2), color='red',  label='$\\tau_{' f'{site2}' '}^z$')
            # ax.legend()
            if fieldinset:
                trans = transforms.blended_transform_factory(ax.transData, ax.transAxes)
                fieldcolor = 'steelblue'
                for site in range(L):
                    ax.arrow(site, 0.12, 0.0, fields[site] / 15.0, width=0.025, length_includes_head=False,
                             head_width=0.00, head_length=0.00,
                             path_effects=pe_fields,
                             zorder=50,
                             transform=trans,
                             color=fieldcolor,
                             # color='salmon'
                             )
                ax.arrow(0, 0.12, L, 0, width=0.001, length_includes_head=False,
                         head_width=0.00, head_length=0.00,
                         zorder=60, transform=trans,
                         color=fieldcolor
                         # color='salmon'
                         )
            print(xi)
            for l in [0,1,2,3]:
                ax.axhline(np.exp(-l/2.0/xi), linestyle='dashed', linewidth=1)

            lbit_max = np.mean(np.max(lbits[:-d], axis=1))
            lbit_maxes.append(lbit_max)
            ax.axhline(lbit_max, color='red', linewidth=1)

except Exception as ex:
    print(ex)
    pass


try:
    ax = fig1.add_subplot(2, 2, 4)
    ax.plot([0,1,2,3,4,5], lbit_maxes)

    # ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylim(ymin=1e-8, ymax=1e1)
    ax.set_xlabel('$l$')
    ax.set_title(f'$\ell$-bits')

except Exception as ex:
    print(ex)
    pass


plt.tight_layout()
plt.show()