from src.io.h5ops import *
import numpy as np
import matplotlib.pyplot as plt


# import seaborn as sns

def get_optimal_subplot_num(numplots):
    r = np.sqrt(numplots)
    cols = int(np.ceil(r))
    rows = int(np.floor(r))
    while cols * rows < numplots:
        if (cols <= rows):
            cols = cols + 1
        else:
            rows = rows + 1
    return rows, cols


def boxplot_SE_vs_site_foreach_Delta(src, plotdir=''):
    print('Plotting:     Boxplot S vs site')
    h5_src = h5open(src, 'r')
    for i, (path_l, node_l) in enumerate(h5py_node_iterator(g=h5_src, keypattern='l_')):
        dsetlist = h5py_node_finder(node_l, keypattern='S_vs_L')
        rows, cols = get_optimal_subplot_num(len(dsetlist))
        fig, axes = plt.subplots(nrows=rows, ncols=cols, figsize=(4 * cols, 4 * rows))
        fig.tight_layout(pad=3, w_pad=1.0, h_pad=1.0)
        fig.subplots_adjust(wspace=0.2, hspace=0.2)
        for ax, dset in zip(np.ravel(axes), dsetlist):
            SE_data = np.array(dset[1])
            ax.boxplot(SE_data.transpose(), whis='range', meanline=True)
            ax.set_title(dset[1].attrs['title'])
            ax.set_xlabel(dset[1].attrs['xlabel'])
            ax.set_ylabel(dset[1].attrs['ylabel'])

        if plotdir != '':
            plt.savefig(plotdir + '/S_vs_L_' + str(i) + '.pdf', format='pdf')
    h5close(h5_src)


def plot_S_vs_Delta(src, plotdir=''):
    print('Plotting:     Entanglement entropy vs Delta')
    h5_src = h5open(src, 'r')
    for i, (path_l, node_l) in enumerate(h5py_node_iterator(g=h5_src, keypattern='l_')):
        dset = h5py_node_finder(node_l, keypattern='S_vs_Delta')[0][1]
        stats = np.array(dset).transpose()
        fig, ax = plt.subplots(nrows=1, ncols=1)
        fig.tight_layout(pad=3, w_pad=1.0, h_pad=1.0)
        fig.subplots_adjust(wspace=0.2, hspace=0.2)
        ax.errorbar(x=stats[:, 0], y=stats[:, 1], yerr=stats[:, 2], label=dset.attrs['col1'])
        ax.plot(np.array(stats[:, 0]), np.array(stats[:, 3]), label=dset.attrs['col3'])
        # ax.plot(np.array(stats[:,0]), np.array(stats[:,4]), label=dset.attrs['col4'])
        ax.set_title('$\lambda = $' + dset.attrs['lambda'])
        ax.set_xlabel(dset.attrs['xlabel'])
        ax.set_ylabel(dset.attrs['ylabel'])
        ax.legend()

        if plotdir != '':
            plt.savefig(plotdir + '/S_vs_Delta_' + str(i) + '.pdf', format='pdf')
    h5close(h5_src)


def histogram_S_foreach_Delta(src, plotdir=''):
    print('Plotting:     Histogram of Entanglement entropy for each Delta')
    h5_src = h5open(src, 'r')
    for i, (path_l, node_l) in enumerate(h5py_node_iterator(g=h5_src, keypattern='l_')):
        grouplist = h5py_node_finder(node_l, keypattern='S_histogram')
        rows, cols = get_optimal_subplot_num(len(grouplist))
        fig, axes = plt.subplots(nrows=rows, ncols=cols, figsize=(4 * cols, 4 * rows))
        fig.tight_layout(pad=3, w_pad=1.0, h_pad=1.0)
        fig.subplots_adjust(wspace=0.2, hspace=0.2)
        for ax, group in zip(np.ravel(axes), grouplist):
            hist = group[1]['S_hist']
            edge = group[1]['S_edges']
            SE_hist = np.array(hist)
            SE_edges = np.array(edge)
            ax.hist(x=SE_hist, bins=SE_edges)
            ax.set_ylim(1e-1, 1e4)
            ax.set_yscale('log', nonpositive='clip')
            ax.set_title(hist.attrs['title'])
            ax.set_xlabel(hist.attrs['xlabel'])
            ax.set_ylabel(hist.attrs['ylabel'])
        if plotdir != '':
            plt.savefig(plotdir + '/S_histogram_' + str(i) + '.pdf', format='pdf')
    h5close(h5_src)


def histogram_var_foreach_Delta(src, plotdir=''):
    print('Plotting:     Histogram of Variance for each Delta')
    h5_src = h5open(src, 'r')
    for i, (path_l, node_l) in enumerate(h5py_node_iterator(g=h5_src, keypattern='l_')):
        grouplist = h5py_node_finder(node_l, keypattern='var_histogram')
        rows, cols = get_optimal_subplot_num(len(grouplist))
        fig, axes = plt.subplots(nrows=rows, ncols=cols, figsize=(4 * cols, 4 * rows))
        fig.tight_layout(pad=3, w_pad=1.0, h_pad=1.0)
        fig.subplots_adjust(wspace=0.2, hspace=0.2)
        for ax, group in zip(np.ravel(axes), grouplist):
            hist = group[1]['var_hist']
            edge = group[1]['var_edges']
            var_hist = np.array(hist)
            var_edges = np.array(edge)
            ax.hist(x=var_hist, bins=var_edges)
            ax.set_ylim(1e-1, 1e4)
            ax.set_yscale('log', nonpositive='clip')
            ax.set_title(hist.attrs['title'])
            ax.set_xlabel(hist.attrs['xlabel'])
            ax.set_ylabel(hist.attrs['ylabel'])
        if plotdir != '':
            plt.savefig(plotdir + '/Var_vs_Delta_' + str(i) + '.pdf', format='pdf')
    h5close(h5_src)


def plot_chi_vs_Delta(src, plotdir=''):
    print('Plotting:     Chi vs Delta')
    h5_src = h5open(src, 'r')
    for i, (path_l, node_l) in enumerate(h5py_node_iterator(g=h5_src, keypattern='l_')):
        dset = h5py_node_finder(node_l, keypattern='Chi_vs_Delta')[0][1]
        stats = np.array(dset).transpose()
        fig, ax = plt.subplots(nrows=1, ncols=1)
        fig.tight_layout(pad=3, w_pad=1.0, h_pad=1.0)
        fig.subplots_adjust(wspace=0.2, hspace=0.2)
        ax.errorbar(x=stats[:, 0], y=stats[:, 1], yerr=stats[:, 2], label=dset.attrs['col1'])
        ax.plot(np.array(stats[:, 0]), np.array(stats[:, 3]), label=dset.attrs['col3'])
        # ax.plot(np.array(stats[:,0]), np.array(stats[:,4]), label=dset.attrs['col4'])
        ax.set_title('$\lambda = $' + dset.attrs['lambda'])
        ax.set_xlabel(dset.attrs['xlabel'])
        ax.set_ylabel(dset.attrs['ylabel'])
        ax.legend()

        if plotdir != '':
            plt.savefig(plotdir + '/Chi_vs_Delta_' + str(i) + '.pdf', format='pdf')
    h5close(h5_src)


def plot_time_vs_Delta(src, plotdir=''):
    print('Plotting:     Time vs Delta')
    h5_src = h5open(src, 'r')
    for i, (path_l, node_l) in enumerate(h5py_node_iterator(g=h5_src, keypattern='l_')):
        dset = h5py_node_finder(node_l, keypattern='Time_vs_Delta')[0][1]
        stats = np.array(dset).transpose()
        fig, ax = plt.subplots(nrows=1, ncols=1)
        fig.tight_layout(pad=3, w_pad=1.0, h_pad=1.0)
        fig.subplots_adjust(wspace=0.2, hspace=0.2)
        ax.errorbar(x=stats[:, 0], y=stats[:, 1], yerr=stats[:, 2], label=dset.attrs['col1'])
        ax.plot(np.array(stats[:, 0]), np.array(stats[:, 3]), label=dset.attrs['col3'])
        # ax.plot(np.array(stats[:,0]), np.array(stats[:,4]), label=dset.attrs['col4'])
        ax.set_ylim(1e-1, 5e4)
        ax.set_yscale('log', nonpositive='clip')
        ax.set_title('$\lambda = $' + dset.attrs['lambda'])
        ax.set_xlabel(dset.attrs['xlabel'])
        ax.set_ylabel(dset.attrs['ylabel'])
        ax.legend()

        if plotdir != '':
            plt.savefig(plotdir + '/Time_vs_Delta_' + str(i) + '.pdf', format='pdf')
    h5close(h5_src)


def plot_var_vs_Delta(src, plotdir=''):
    print('Plotting:     Var vs Delta')
    h5_src = h5open(src, 'r')
    for i, (path_l, node_l) in enumerate(h5py_node_iterator(g=h5_src, keypattern='l_')):
        dset = h5py_node_finder(node_l, keypattern='Var_vs_Delta')[0][1]
        stats = np.array(dset).transpose()
        fig, ax = plt.subplots(nrows=1, ncols=1)
        fig.tight_layout(pad=3, w_pad=1.0, h_pad=1.0)
        fig.subplots_adjust(wspace=0.2, hspace=0.2)
        ax.errorbar(x=stats[:, 0], y=stats[:, 1], yerr=stats[:, 2], label=dset.attrs['col1'])
        ax.plot(np.array(stats[:, 0]), np.array(stats[:, 3]), label=dset.attrs['col3'])
        # ax.plot(np.array(stats[:,0]), np.array(stats[:,4]), label=dset.attrs['col4'])
        ax.set_ylim(1e-14, 3e2)
        ax.set_yscale('log', nonpositive='clip')
        ax.set_title('$\lambda = $' + dset.attrs['lambda'])
        ax.set_xlabel(dset.attrs['xlabel'])
        ax.set_ylabel(dset.attrs['ylabel'])
        ax.legend()

        if plotdir != '':
            plt.savefig(plotdir + '/Time_vs_Delta_' + str(i) + '.pdf', format='pdf')
    h5close(h5_src)
