import numpy as np
import h5py
import matplotlib.pyplot as plt
from scipy.stats import gmean
import seaborn as sns

h5files = ['../../output/pert-L16-g1e-6.h5']
sns.set_palette("Dark2")
fig,axes = plt.subplots(nrows=2, ncols=1)
for h5file, marker, color, label in zip(h5files, ['o', 'v', 'x'], ['blue', 'green', 'red'], ['$L=12$', '$L=14$', '$L=16$']):
    with h5py.File(h5file, 'r') as f:
        axes[1].set_ylabel('Energy variance')
        axes[1].set_xlabel('$\delta$')
        # Get a list of all deltas
        deltas = set()
        for idx, (key, val) in enumerate(f.items()):
            deltas.add(val['delta'][0])
        deltas = np.array(sorted(list(deltas)))
        deltas_interp = []
        for didx in range(len(deltas)):
            deltas_interp.append(deltas[didx])
            if didx+1 < len(deltas):
                deltas_interp.append(np.mean([deltas[didx+1], deltas[didx]]))

        deltas_interp = np.array(deltas_interp)
        print(deltas)
        print(deltas_interp)

        # Set the boundaries for each delta as the halfway point between two deltas
        bounds = []
        for didx in range(len(deltas)):
            lbound = -np.inf if didx == 0 else np.mean([deltas[didx-1], deltas[didx]])
            ubound = +np.inf if didx+1 == len(deltas) else np.mean([deltas[didx], deltas[didx+1]])
            # lbound = -np.inf if didx == 0 else deltas[didx]-0.01
            # ubound = +np.inf if didx+1 == len(deltas) else deltas[didx]+0.01
            bounds.append((lbound,ubound))

        # Sort the variance values
        var_prts = [[] for _ in range(len(deltas))]
        for key, val in f.items():
            delta_avg = val['delta_avg'][()]
            delta = val['delta'][0]
            var_prt = val['variance_pert'][()]
            # Get the delta index
            idx = np.where(deltas == delta)[0][0]
            for bidx, bs in enumerate(bounds):
                didx = np.where((bs[0] <= delta_avg) & (delta_avg < bs[1]))[0]
                if len(didx) > 0:
                    print(f'in delta={delta} within {bs=} going to delta={deltas[bidx]}: {len(didx)=}')
                    var_prts[bidx].extend(var_prt[didx])


        for bidx, var_prt in enumerate(var_prts):
            g_pert = val['g_pert'][0]
            var_avg = np.median(np.abs(var_prt))
            delta = deltas[bidx]
            print(f'{delta=} {len(var_prt)=}')
            line1 = axes[1].errorbar(x=delta, y=var_avg, yerr=None, color='red', marker='v')



        for didx, (key, val) in enumerate(f.items()):
            g_pert = val['g_pert'][0]

            varch_num = len(val['variance'][()])
            var_old = np.abs(val['variance'][()])
            var_prt = np.abs(val['variance_pert'][()])
            varch_avg = np.median(np.abs(var_prt))
            varch_std = np.std(np.abs(var_prt))
            varch_ste = varch_std/np.sqrt(varch_num)
            #
            # enech_num = len(val['energy'][()])
            # ene_old = np.abs(val['energy'][()])
            # ene_prt = np.abs(val['energy_pert'][()])
            # enech_avg = np.median(np.abs(ene_old - ene_prt))
            # enech_ste = np.std(np.abs(ene_old - ene_prt))/np.sqrt(enech_num)


            delta = val['delta'][0]
            delta_avg = val['delta_avg'][()]
            hist, edges = np.histogram(delta_avg, bins=150, range=(-6.5,6.5), density=True)
            bincentres = [(edges[j] + edges[j + 1]) / 2. for j in range(len(edges) - 1)]
            # line0, = axes.step(x=bincentres, y=hist, where='mid', label=None,
                            # linewidth=0.85,
                            # )
            line0, = axes[0].plot(bincentres, hist,label=None)


            line1 = axes[1].errorbar(x=delta, y=varch_avg, yerr=None, color=color, marker=marker)


            # line1 = axes[1].errorbar(x=delta, y=enech_avg, yerr=enech_ste, color=color, marker=marker)
            # print(f'{key}: {varch_avg:.3e} +- {varch_ste:.3e}')

        axes[0].axvline(x=np.log(2), color='black')
        axes[0].axvline(x=-np.log(2), color='black')
        axes[0].set_xlabel('$\delta$')
        axes[0].set_ylabel("$P(\\bar\delta)$")
        axes[0].set_title(f'Distribution of "actual" $\delta$-values')
        # axes[0].set_yscale('log')

        # axes[1].axvline(x=np.log(2))
        # axes[1].axvline(x=-np.log(2))
        # axes[1].set_xlabel('$\delta$')
        # axes[1].set_ylabel('$\\mathrm{median}(|\langle H \\rangle - \langle H^\prime \\rangle|)$')
        # axes[1].set_title(f'Change in energy when $g = 0\\rightarrow 10^{{{int(np.log10(g_pert))}}}$')

        line0.set_label(label)
        line1.set_label(label)
        axes[0].legend(loc='best')
        axes[1].legend(loc='best')
plt.tight_layout()
# plt.legend(loc='best')
plt.show()
