import h5py
import numpy as np
import matplotlib.pyplot as plt
from src.h5py_extra import *

plt.close('all')
filename = 'data8.h5'
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'], 'size':12})
rc('text', usetex=True)


with h5py.File(filename, 'r') as f:
    group = 'FES_iDMRG'
    chi          = h5py_get_datasets(f, group, 'chi')
    chi_max      = h5py_get_datasets(f, group, 'chi_max')
    length       = h5py_get_datasets(f, group, 'length')
    energy       = h5py_get_datasets(f, group, 'energy')
    entropy      = h5py_get_datasets(f, group, 'entropy')
    error        = h5py_get_datasets(f, group, 'energy_error')
    trunc_error  = h5py_get_datasets(f, group, 'trunc_error')
    gitrevision  = dict(f.attrs.items())
    print(gitrevision)


# Two subplots, the axes array is 1-d
fig, axarr = plt.subplots(1, sharex = False)
fig.set_size_inches(4, 3, forward=True)

for l, e, c in zip(length, energy, chi_max):
    axarr.scatter(l[0], e[0],marker='o',  s=3, label='$\chi_{max} =' + str(c[0][0]) + '$')
axarr.axhline(y=-1.2732395447351625, xmin=0, xmax=1, c="blue", linewidth=0.8, zorder=0, label='Exact')
axarr.set_ylabel('Ground State Energy')
axarr.set_xlabel('Chain length')
axarr.legend(bbox_to_anchor=(1, 0.1), loc='lower right', ncol=1)

fig.tight_layout()


fig, axarr = plt.subplots(1, sharex = False)
fig.set_size_inches(4, 3, forward=True)

for l, e, c in zip(length, entropy, chi_max):
    axarr.scatter(l[0], e[0],marker='o',  s=3, label='$\chi_{max} =' + str(c[0][0]) + '$')
axarr.set_ylabel('Entanglement Entropy')
axarr.set_xlabel('Chain length')
axarr.legend(bbox_to_anchor=(1, 0.1), loc='lower right', ncol=1)
fig.tight_layout()

plt.show()

#
#
# for l, e, c in zip(length, entropy, chi_max):
#     axarr[2,0].scatter(l[0], e[0],marker='o',  s=2, label='chi=' + str(c[0][0]))
#     axarr[2,0].axhline(y=np.log(c[0][0])/(np.sqrt(24)+1), xmin=0, xmax=1, c="blue", linewidth=0.5, zorder=0, label='Exact ($\chi$ = ' + str(c[0][0]) +')')
#
# axarr[2,0].legend(loc='lower right')
# axarr[2,0].set_ylabel('Entropy')
# axarr[2,0].set_xlabel('Chain Length')




# fig.subplots_adjust(hspace=0.2, wspace=0.2)



