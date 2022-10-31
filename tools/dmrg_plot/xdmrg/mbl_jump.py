from __future__ import unicode_literals
from src.plotting.style import *
import glob
import numpy as np
import matplotlib.pyplot as plt
from src.plotting.style import *

projdir = '/mnt/WDB-AN1500/mbl_transition'
batchdir = glob.glob(projdir + '/data165*')[0]
analysisdir = batchdir + '/analysis'
plotdir = analysisdir + '/plots'
datadir = analysisdir + '/data'

seed = 50002596
data = np.genfromtxt('{}/mbl_{}.txt'.format(datadir, seed), delimiter=',')
bond = data[0]
enpy = data[1]
evar = data[2]
iter = data[3]
satr = data[4]

figrows, figcols = (2, 2)  # Add one for the legend
fig, ax = plt.subplots(nrows=figrows, ncols=figcols, figsize=(4 * figcols, 4 * figrows))
fig.tight_layout(pad=5, w_pad=1.0, h_pad=1.0)
fig.subplots_adjust(wspace=0.2, hspace=0.2)
xticks = [0, 66, 85, 113, 149, 151, 153, 155, 157, 159, 161, 163, 165, 167, 169, 171, 172]
xlabel = [0, 66, 85, 113, 149, '', '', '', '', '', '', '', '', '', '', '', '']
xlim = [100, 165]
# 66,85,113,149,151,153,155,157,159,161,163,165,167,169,171,172


ax[0, 0].plot(iter, bond, label='$\chi_\mathrm{limit}$')
ax[0, 0].set_xlabel('iter')
ax[0, 0].set_yticks([8, 16, 24, 32, 40, 48, 56, 64, 72, 80, 88, 96, 104, 112, 120, 128])
ax[0, 0].set_xticks(xticks)
ax[0, 0].set_xticklabels(xlabel)
# ax[0,0].set_xlim(xlim)
ax[0, 0].legend()

ax[0, 1].plot(iter, evar, label='$\sigma^2(H)$')
ax[0, 1].set_yscale('log')
ax[0, 1].set_xlabel('iter')
ax[0, 1].set_xticks(xticks)
ax[0, 1].set_xticklabels(xlabel)
# ax[0,1].set_xlim(xlim)
ax[0, 1].legend()

ax[1, 0].plot(iter, enpy, label='$S(L/2)$')
ax[1, 0].set_xlabel('iter')
ax[1, 0].set_xticks(xticks)
ax[1, 0].set_xticklabels(xlabel)
# ax[1,0].set_ylim([1.44, 1.48])
# ax[1,0].set_xlim(xlim)
ax[1, 0].legend()

ax[1, 1].plot(iter, satr, label='Saturated for $n$ iters')
ax[1, 1].set_xlabel('iter')
ax[1, 1].set_xticks(xticks)
ax[1, 1].set_xticklabels(xlabel)
ax[1, 1].set_yticks([0, 1, 2])
# ax[1,1].set_xlim(xlim)
ax[1, 1].legend()

fig.suptitle('MBL seed {}'.format(seed))
# fig.delaxes(ax[1,1])

plt.show()
# print(np.genfromtxt('{}/mbl_50002596.txt'.format(datadir), delimiter=','))
# with open('{}/mbl_50002596.txt'.format(datadir), 'r') as f:
#     for line in f.readlines():
#         print(np.genfromtxt(line, ','))
