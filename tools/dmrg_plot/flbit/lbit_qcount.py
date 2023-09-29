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
    Ls = np.arange(2,32,2)
    qdists = np.zeros((len(Ls),20))
    for l, L in enumerate(Ls):
        qbinc = get_qdistribution(L)
        qdist = qbinc/np.sum(qbinc)
        qdists[l,:] = qdist
        # ax.scatter(x=[L], y=[qdist[1]], color='blue')
        # if len(qdist) > 2:
    fig = plt.figure(figsize=(6,6))
    ax1 = fig.add_subplot(2, 2, 1)
    # ax1.set_yscale('log')
    # ax1.set_xscale('log')
    ax1.scatter(x=Ls, y=qdists[:, 1], color='blue', label="$q_1$")
    ax1.scatter(x=Ls, y=qdists[:, 2], color='red', label="$q_2$")
    ax1.scatter(x=Ls, y=qdists[:, 3], color='orange', label="$q_3$")
    ax1.set_xlabel('$L$')
    ax1.set_ylabel('$q_n$')
    ax1.legend(loc='upper right')

    ax2 = fig.add_subplot(2, 2, 2)
    ax2.scatter(x=1/Ls, y=qdists[:, 1], color='blue', label=None)
    ax2.scatter(x=1/Ls, y=qdists[:, 2], color='red', label=None)
    ax2.scatter(x=1/Ls, y=qdists[:, 3], color='orange', label=None)
    popt1, _ = curve_fit(f=flinear, xdata=1/Ls[-4:-1], ydata=qdists[-4:-1,1])
    popt2, _ = curve_fit(f=flinear, xdata=1/Ls[-4:-1], ydata=qdists[-4:-1,2])
    popt3, _ = curve_fit(f=flinear, xdata=1/Ls[-4:-1], ydata=qdists[-4:-1,3])
    xfit=np.linspace(0,0.5,100)
    yfit1 = flinear(xfit, *popt1)
    yfit2 = flinear(xfit, *popt2)
    yfit3 = flinear(xfit, *popt3)
    ax2.plot(xfit, yfit1, marker=None, linewidth=0.4, linestyle='--', label=f'$q_1^{{\infty}} = {yfit1[0]:.3f}$', color='blue')
    ax2.plot(xfit, yfit2, marker=None, linewidth=0.4, linestyle='--', label=f'$q_2^{{\infty}} = {yfit2[0]:.3f}$', color='red')
    ax2.plot(xfit, yfit3, marker=None, linewidth=0.4, linestyle='--', label=f'$q_3^{{\infty}} = {yfit3[0]:.3f}$', color='orange')
    ax2.set_xlabel('$1/L$')
    ax2.set_ylabel('$q_n$')
    ax2.legend(loc='upper right')


    ax3 = fig.add_subplot(2, 2, 3)
    ax3.scatter(x=Ls, y=qdists[:, 2] / qdists[:, 1], color='blue', label="$q_2/q_1$")
    ax3.scatter(x=Ls, y=qdists[:, 3] / qdists[:, 2], color='red', label="$q_3/q_2$")
    ax3.scatter(x=Ls, y=qdists[:, 4] / qdists[:, 3], color='orange', label="$q_4/q_3$")
    ax3.set_xlabel('$L$')
    ax3.set_ylabel('$q$-fractions')
    ax3.legend(loc='center right')

    ax4 = fig.add_subplot(2, 2, 4)
    ax4.scatter(x=Ls, y=qdists[:, 2] / qdists[:, 1], color='blue', label="$q_2/q_1$")
    ax4.scatter(x=Ls, y=qdists[:, 3] / qdists[:, 1], color='red', label="$q_3/q_1$")
    ax4.scatter(x=Ls, y=qdists[:, 4] / qdists[:, 1], color='orange', label="$q_4/q_1$")
    ax4.set_xlabel('$L$')
    ax4.set_ylabel('$q$-fractions')
    ax4.legend(loc='center right')

    plt.tight_layout()
    plt.show()