import numpy as np
import matplotlib.pyplot as plt
from src.measurement.compute_values import *
from src.measurement.compute_statistics import *
import seaborn as sns


def plot_avg_entanglement_entropy(hdf5_set, ax, title_str=''):
    SE_data = mps_entanglement_entropy_statistics(hdf5_set['/data'], full_chain=True)
    # SE_data = SE_data.transpose()
    pos = np.arange(0, np.shape(SE_data)[1])
    # ax.plot(range(0, len(SE_avg)), SE_avg,
    #          'k', color='#1B2ACC', label='$\overline{S_E(l)}$')
    ax.set_xlabel('L')
    ax.set_ylabel('Entanglement Entropy $S_E$')

    ax = sns.boxplot(data=SE_data)

    # ax.boxplot(np.array(SE_data).T.tolist())
    # ax.fill_between(range(0, len(SE_avg)), SE_avg + SE_std / 2, SE_avg - SE_std / 2,
    #                  alpha=0.2, edgecolor='#1B2ACC', facecolor='#089FFF',
    #                  linewidth=4, linestyle='dashdot', antialiased=True,
    #                  label='$\overline{S_E(l)} \pm \\frac{\sigma}{2}$'
    #                  )
    ax.set_title(title_str)
    ax.legend()


def plot_avg_S(S_list, title_str):
    n = len(S_list)
    L = len(S_list[0])
    SE = np.zeros([n, L])
    for i in range(0, n):
        SE[i, :] = mps_entanglement_entropy(S_list[i][:])
        # print(mps_entanglement_entropy(S_lists[:][i]))
    SE = mps_entanglement_entropy_statistics(S_list)
    SE_avg = np.mean(SE, axis=0)
    SE_std = np.std(SE, axis=0)
    plt.plot(range(0, len(SE_avg)), SE_avg,
             'k', color='#1B2ACC', label='$\overline{S_E(l)}$')
    plt.xlabel('L')
    plt.ylabel('Entanglement Entropy $S_E$')
    plt.fill_between(range(0, len(SE_avg)), SE_avg + SE_std / 2, SE_avg - SE_std / 2,
                     alpha=0.2, edgecolor='#1B2ACC', facecolor='#089FFF',
                     linewidth=4, linestyle='dashdot', antialiased=True,
                     label='$\overline{S_E(l)} \pm \\frac{\sigma}{2}$'
                     )
    plt.title(title_str)
    plt.legend()
