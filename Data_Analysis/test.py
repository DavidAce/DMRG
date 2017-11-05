import h5py
import numpy as np
import matplotlib.pyplot as plt
from src.h5py_extra import *
plt.close('all')

filename = 'data7.h5'

with h5py.File(filename, 'r') as f:
    group = 'FES_iDMRG'
    chi          = h5py_get_datasets(f, group, 'chi')
    chi_max      = h5py_get_datasets(f, group, 'chi_max')
    length       = h5py_get_datasets(f, group, 'length')
    # delta_t      = h5py_get_datasets(f, group, 'delta_t')
    # wall_time    = h5py_get_datasets(f, group, 'wall_time')
    # phys_time    = h5py_get_datasets(f, group, 'phys_time')
    energy       = h5py_get_datasets(f, group, 'energy')
    entropy      = h5py_get_datasets(f, group, 'entropy')
    trunc_error  = h5py_get_datasets(f, group, 'trunc_error')
    gitrevision  = dict(f.attrs.items())

print(gitrevision)


# Two subplots, the axes array is 1-d
f1, axarr = plt.subplots(2, sharex = True)

for l, e, c in zip(length, energy, chi_max):
    axarr[0].scatter(l[0], e[0],marker='o',  s=3, label='chi=' + str(c[0][0]))
axarr[0].set_ylabel('Energy')
axarr[0].legend(loc='upper right')


for l, e, c in zip(length, entropy, chi_max):
    axarr[1].scatter(l[0], e[0],marker='o',  s=3, label='chi=' + str(c[0][0]))

axarr[1].legend(loc='lower right')
axarr[1].set_ylabel('Entropy')
axarr[1].set_xlabel('Chain Length')

f1.subplots_adjust(hspace=0.2, wspace=0.2)


with h5py.File(filename, 'r') as f:
    group = 'FES_iTEBD'
    chi          = h5py_get_datasets(f, group, 'chi')
    chi_max      = h5py_get_datasets(f, group, 'chi_max')
    length       = h5py_get_datasets(f, group, 'length')
    delta_t      = h5py_get_datasets(f, group, 'delta_t')
    wall_time    = h5py_get_datasets(f, group, 'wall_time')
    phys_time    = h5py_get_datasets(f, group, 'phys_time')
    energy       = h5py_get_datasets(f, group, 'energy')
    entropy      = h5py_get_datasets(f, group, 'entropy')
    trunc_error  = h5py_get_datasets(f, group, 'trunc_error')
    gitrevision  = dict(f.attrs.items())



# Two subplots, the axes array is 1-d
f2, axarr = plt.subplots(2, sharex = True)


for t, e, c in zip(phys_time, energy, chi_max):
    axarr[0].scatter(t[0], e[0],marker='.',  s=1, label='chi=' + str(c[0][0]))
axarr[0].set_ylabel('Energy')
axarr[0].legend(loc='upper right')


for t, e, c in zip(phys_time, entropy, chi_max):
    axarr[1].scatter(t[0], e[0],marker='.',  s=1, label='chi=' + str(c[0][0]))

axarr[1].legend(loc='lower right')
axarr[1].set_ylabel('Entropy')
axarr[1].set_xlabel('Time')

f2.subplots_adjust(hspace=0.2, wspace=0.2)



plt.show()




