import h5py
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from src.h5py_extra import *
import seaborn as sns
plt.close('all')
# sns.set_style("darkgrid")
sns.set(font_scale=1.5)  # crazy big

# sns.set(style="ticks")
# print(sns.style.available)

filename = '../output/data.h5'
# filename = 'data12.h5'

with h5py.File(filename, 'r') as f:
    group = 'FES_iDMRG'
    chi          = h5py_get_datasets(f, group, 'chi')
    chi_max      = h5py_get_datasets(f, group, 'chi_max')
    length       = h5py_get_datasets(f, group, 'length')
    energy       = h5py_get_datasets(f, group, 'energy')
    entropy      = h5py_get_datasets(f, group, 'entropy')
    variance1      = h5py_get_datasets(f, group, 'variance1')
    variance2      = h5py_get_datasets(f, group, 'variance2')
    variance3      = h5py_get_datasets(f, group, 'variance3')
    error        = h5py_get_datasets(f, group, 'energy_error')
    trunc_error  = h5py_get_datasets(f, group, 'trunc_error')
    gitrevision  = dict(f.attrs.items())
    print(gitrevision)


# Two subplots, the axes array is 1-d
fig, axarr = plt.subplots(3, 2, sharex = True)

for l, e, c in zip(length, energy, chi_max):
    axarr[0,0].scatter(l[0], e[0],marker='o',  s=4, label='$\chi=$' + str(c[0][0]))
axarr[0,0].axhline(y=-1.2732395447351625, xmin=0, xmax=1, c="blue", linewidth=0.5, zorder=0, label='$\chi = \infty$')
axarr[0,0].set_ylabel('Energy')
axarr[0,0].legend(loc='upper right')
axarr[0,0].set_title('iDMRG')

# for l, e, c in zip(length, error, chi_max):
#     axarr[1,0].scatter(l[0], e[0],marker='o',  s=3, label='chi=' + str(c[0][0]))
# #
# axarr[1,0].set_ylabel('Error')
# axarr[1,0].legend(loc='upper right')

for l, e, c in zip(length, entropy, chi_max):
    axarr[1,0].scatter(l[0], e[0],marker='o',  s=3, label='$\chi=$' + str(c[0][0]))
    axarr[1,0].axhline(y=np.log(c[0][0])/(np.sqrt(24)+1), xmin=0, xmax=1, c='blue', linewidth=0.5, zorder=0, label='Analytic ($\chi = $ ' + str(c[0][0]) +')')

axarr[1,0].legend(loc='lower right')
axarr[1,0].set_ylabel('Entropy')
# axarr[1,0].set_xlabel('Chain Length')

for l, v1,v2,v3, c in zip(length, variance1, variance2, variance3, chi_max):
    axarr[2,0].scatter(l[0], v1[0],marker='o',  s=4, label='$\langle H^2 \\rangle - \langle E \\rangle^2$, ' + '$\chi=$' + str(c[0][0]))
    axarr[2,0].scatter(l[0], v3[0],marker='o',  s=4, label='$\log{|G(a)|^2}$, ' + '$\chi=$' + str(c[0][0]))
    # axarr[2,0].axhline(y=np.log(c[0][0])/(np.sqrt(24)+1), xmin=0, xmax=1, c='blue', linewidth=0.5, zorder=0, label='Analytic ($\chi = $ ' + str(c[0][0]) +')')

axarr[2,0].legend(loc='lower right')
axarr[2,0].set_ylabel('Variance')
axarr[2,0].set_xlabel('Chain Length')
fig.tight_layout()
fig.subplots_adjust(wspace=.4, hspace=.2)
plt.show()

exit(0)

with h5py.File(filename, 'r') as f:
    group = 'FES_iTEBD'
    chi          = h5py_get_datasets(f, group, 'chi')
    chi_max      = h5py_get_datasets(f, group, 'chi_max')
    length       = h5py_get_datasets(f, group, 'length')
    time_step    = h5py_get_datasets(f, group, 'time_step')
    wall_time    = h5py_get_datasets(f, group, 'wall_time')
    phys_time    = h5py_get_datasets(f, group, 'phys_time')
    energy       = h5py_get_datasets(f, group, 'energy')
    entropy      = h5py_get_datasets(f, group, 'entropy')
    variance     = h5py_get_datasets(f, group, 'variance')
    error        = h5py_get_datasets(f, group, 'energy_error')
    trunc_error  = h5py_get_datasets(f, group, 'trunc_error')
    gitrevision  = dict(f.attrs.items())
    print(gitrevision)


for t, e, c in zip(phys_time, energy, chi_max):
    axarr[0,1].scatter(t[0], e[0],marker='o',  s=4, label='$\chi=$' + str(c[0][0]))
    axarr[0,1].axhline(y=-1.2732395447351625, xmin=0, xmax=1, c="blue", linewidth=0.5, zorder=0, label='$\chi = \infty$')

axarr[0,1].set_ylabel('Energy')
axarr[0,1].legend(loc='upper right')
axarr[0,1].set_title('iTEBD')

# for t, e, c in zip(phys_time, error, chi_max):
#     axarr[1,1].scatter(t[0], e[0],marker='.',  s=3, label='chi=' + str(c[0][0]))
# axarr[1,1].set_ylabel('Error')
# axarr[1,1].legend(loc='upper right')

for t, e, c in zip(phys_time, entropy, chi_max):
    axarr[1,1].scatter(t[0], e[0],marker='o',  s=3, label='$\chi=$' + str(c[0][0]))
    axarr[1,1].axhline(y=np.log(c[0][0])/(np.sqrt(24)+1), xmin=0, xmax=1, c="blue", linewidth=0.5, zorder=0, label='Analytic ($\chi = $ ' + str(c[0][0]) +')')

axarr[1,1].legend(loc='lower right')
axarr[1,1].set_ylabel('Entropy')
axarr[1,1].set_xlabel('Time')


# fig.subplots_adjust(hspace=0.2, wspace=0.2)
fig.tight_layout()
plt.show()




