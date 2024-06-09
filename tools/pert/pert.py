import numpy as np
import h5py
import matplotlib.pyplot as plt
from scipy.stats import gmean

h5files = ['../../output/pert-L12-g1e-5.h5', '../../output/pert-L14-g1e-5.h5', '../../output/pert-L16-g1e-5.h5']
fig,axes = plt.subplots(nrows=2, ncols=1)
for h5file, marker, color, label in zip(h5files, ['o', 'v', 'x'], ['blue', 'green', 'red'], ['$L=12$', '$L=14$', '$L=16$']):
    with h5py.File(h5file, 'r') as f:
        for key, val in f.items():
            g_pert = val['g_pert'][0]
            varch_num = len(val['variance'][()])
            var_old = np.abs(val['variance'][()])
            var_prt = np.abs(val['variance_pert'][()])
            varch_avg = np.median(np.abs(var_prt))
            varch_ste = np.std(np.abs(var_prt))/np.sqrt(varch_num)

            enech_num = len(val['energy'][()])
            ene_old = np.abs(val['energy'][()])
            ene_prt = np.abs(val['energy_pert'][()])
            enech_avg = np.median(np.abs(ene_old - ene_prt))
            enech_ste = np.std(np.abs(ene_old - ene_prt))/np.sqrt(enech_num)


            delta = val['delta'][0]



            line0 = axes[0].errorbar(x=delta, y=varch_avg, yerr=None, color=color, marker=marker)


            line1 = axes[1].errorbar(x=delta, y=enech_avg, yerr=enech_ste, color=color, marker=marker)
            print(f'{key}: {varch_avg:.3e} +- {varch_ste:.3e}')

        axes[0].axvline(x=np.log(2))
        axes[0].axvline(x=-np.log(2))
        axes[0].set_xlabel('$\delta$')
        axes[0].set_ylabel("$\\mathrm{median} \left( \\frac{\\mathrm{Var}(H)}{\\mathrm{Var}(H^\prime)} \\right)$")
        axes[0].set_title(f'Change in energy variance when $g = 0\\rightarrow 10^{{{int(np.log10(g_pert))}}}$')
        axes[0].set_yscale('log')

        axes[1].axvline(x=np.log(2))
        axes[1].axvline(x=-np.log(2))
        axes[1].set_xlabel('$\delta$')
        axes[1].set_ylabel('$\\mathrm{median}(|\langle H \\rangle - \langle H^\prime \\rangle|)$')
        axes[1].set_title(f'Change in energy when $g = 0\\rightarrow 10^{{{int(np.log10(g_pert))}}}$')

        line0.set_label(label)
        line1.set_label(label)
        axes[0].legend(loc='best')
        axes[1].legend(loc='best')
plt.tight_layout()
# plt.legend(loc='best')
plt.show()
