import numpy as np
from scipy.special import comb
from more_itertools import distinct_permutations
import numba as nb
import matplotlib.pyplot as plt
from sympy.utilities.iterables import multiset_permutations
from scipy.optimize import curve_fit


def flinear(x, a, b):
    with np.errstate(invalid='ignore'):
        return a + b * x

@nb.njit(cache=True)
def get_qvalue(state):
    L = len(state)
    q = 1
    b = None
    for i,j in zip(range(L//2-1,-1,-1), range(L//2, L)):
        if b is None:
            b = state[i]
        if state[i] != state[j]:
            return q
        if state[i] != b:
            return q
        q += 1
    return 0

# @njit(cache=True)
def get_qdistribution(L):
    state_init = nb.typed.List(np.resize([0,1], L))
    # qlist = []
    qlist = np.zeros(20)
    for idx,state in enumerate(distinct_permutations(state_init)):
    # for idx, state in enumerate(multiset_permutations(state_init)):
    # for idx, state in enumerate(perm_unique_fast(state_init)):
            qlist[get_qvalue(state)] += 1
        # qlist.append(get_qvalue(state))
        # if idx >= 100000:
        #     break
    return qlist

        # print(state, ":",get_qvalue(state))

# if
if __name__ == '__main__':
    # ax.set_yscale('log')
    # ax.set_xlim(xmin=1,xmax=10)
    plt.style.use('/home/david/GitProjects/DMRG++/tools/dmrg_plot/common/stylesheets/prb.mplstyle')
    Ls = np.arange(2,32,2)
    qdists = np.zeros((len(Ls),20))
    figsize = (3.404, 3.404)
    for l, L in enumerate(Ls):
        qbinc = get_qdistribution(L)
        qdist = qbinc/np.sum(qbinc)
        qdists[l,:] = qdist
        # ax.scatter(x=[L], y=[qdist[1]], color='blue')
        # if len(qdist) > 2:
    fig1,axes = plt.subplots(figsize=figsize,nrows=1, ncols=2, layout='constrained', sharey='all')
    ax1 = axes[0]
    ax1.set_box_aspect(1)
    # ax1.set_yscale('log')
    # ax1.set_xscale('log')
    ax1.scatter(x=Ls, y=qdists[:, 1], color='royalblue', label="$Q_1$", s = 15, linewidth=0.3, edgecolor='black')
    ax1.scatter(x=Ls, y=qdists[:, 2], color='tomato', label="$Q_2$", s = 15, linewidth=0.3, edgecolor='black')
    ax1.scatter(x=Ls, y=qdists[:, 3], color='orange', label="$Q_3$", s = 15, linewidth=0.3, edgecolor='black')
    ax1.set_xlabel('$L$')
    ax1.set_ylabel('$Q_n$')
    ax1.legend(loc='upper right',frameon=True, framealpha=0.7, facecolor='white')

    # fig2 = plt.figure(figsize=figsize, layout='constrained')
    # ax2 = fig1.add_subplot(1, 2, 2)
    ax2 = axes[1]
    ax2.set_box_aspect(1)
    ax2.set_ylim(ymax=1.05, ymin=-0.05)

    popt1, _ = curve_fit(f=flinear, xdata=1/Ls[-4:-1], ydata=qdists[-4:-1,1])
    popt2, _ = curve_fit(f=flinear, xdata=1/Ls[-4:-1], ydata=qdists[-4:-1,2])
    popt3, _ = curve_fit(f=flinear, xdata=1/Ls[-4:-1], ydata=qdists[-4:-1,3])
    xfit=np.linspace(0,0.5,100)
    yfit1 = flinear(xfit, *popt1)
    yfit2 = flinear(xfit, *popt2)
    yfit3 = flinear(xfit, *popt3)

    ax2.plot(xfit, yfit1, marker=None, linewidth=1.0, linestyle='-', label=f'$Q_1^{{\infty}} = {yfit1[0]:.3f}$', color='royalblue',zorder=1)
    ax2.plot(xfit, yfit2, marker=None, linewidth=1.0, linestyle='-', label=f'$Q_2^{{\infty}} = {yfit2[0]:.3f}$', color='tomato',zorder=1)
    ax2.plot(xfit, yfit3, marker=None, linewidth=1.0, linestyle='-', label=f'$Q_3^{{\infty}} = {yfit3[0]:.3f}$', color='orange',zorder=1)
    ax2.scatter(x=1/Ls[::-1], y=qdists[::-1, 1], color='royalblue', label=None, s=15, linewidth=0.3, edgecolor='black',zorder=2)
    ax2.scatter(x=1/Ls[::-1], y=qdists[::-1, 2], color='tomato', label=None, s=15, linewidth=0.3, edgecolor='black',zorder=2)
    ax2.scatter(x=1/Ls[::-1], y=qdists[::-1, 3], color='orange', label=None, s=15, linewidth=0.3, edgecolor='black',zorder=2)
    ax2.set_xlabel('$1/L$')
    # ax2.set_ylabel('$Q_n$')
    ax2.legend(loc='upper left',frameon=True, framealpha=0.7, facecolor='white')
    plt.show()
    exit(0)
    fig3 = plt.figure(figsize=figsize, layout='constrained')
    ax3 = fig3.add_subplot(1, 1, 1)
    ax3.set_box_aspect(1)
    ax3.scatter(x=Ls, y=qdists[:, 2] / qdists[:, 1], color='royalblue', label="$Q_2/Q_1$")
    ax3.scatter(x=Ls, y=qdists[:, 3] / qdists[:, 2], color='tomato', label="$Q_3/Q_2$")
    ax3.scatter(x=Ls, y=qdists[:, 4] / qdists[:, 3], color='orange', label="$Q_4/Q_3$")
    ax3.set_xlabel('$L$')
    ax3.set_ylabel('$Q$-fractions')
    ax3.legend(loc='center right')

    fig4 = plt.figure(figsize=figsize, layout='constrained')
    ax4 = fig4.add_subplot(1, 1, 1)
    ax4.set_box_aspect(1)
    ax4.scatter(x=Ls, y=qdists[:, 2] / qdists[:, 1], color='royalblue', label="$Q_2/Q_1$")
    ax4.scatter(x=Ls, y=qdists[:, 3] / qdists[:, 1], color='tomato', label="$Q_3/Q_1$")
    ax4.scatter(x=Ls, y=qdists[:, 4] / qdists[:, 1], color='orange', label="$Q_4/Q_1$")
    ax4.set_xlabel('$L$')
    ax4.set_ylabel('$Q$-fractions')
    ax4.legend(loc='center right')

    # plt.tight_layout()
    plt.show()