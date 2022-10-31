from src.plotting.tools import *
from src.io.h5ops import *
import matplotlib.pyplot as plt
from src.plotting.filter import *
from scipy.stats import norm


def multiplot_delta_distribution(src, plotdir='', key_list='', type='typical'):
    print('Plotting: Delta distribution')
    h5_src = h5open(src, 'r')
    path_L = h5py_unique_finder(h5_src, filter='L_', dep=1)
    path_l = h5py_unique_finder(h5_src, filter='l_', dep=2)
    path_J = h5py_unique_finder(h5_src, filter='J_', dep=3)
    path_h = h5py_unique_finder(h5_src, filter='h_', dep=4)
    path_d = []
    for J in path_J:
        for h in path_h:
            path_d.append(J + '/' + h)
    timecounter = 0
    simcounter = 0
    # One figure per unique_l, unique_J and unique_h
    try:
        for l in path_l:
            for d in path_d:
                # In each figure we want one subplot per unique_L
                rows, cols = get_optimal_subplot_num(len(path_L))
                fig, axes = plt.subplots(nrows=rows, ncols=cols, figsize=(7 * cols, 7 * rows))
                fig.tight_layout(pad=5, w_pad=1.0, h_pad=1.0)
                fig.subplots_adjust(wspace=0.3, hspace=0.3)
                used_ax = 0
                delt = 0
                lamb = 0
                for ax, L in zip(np.ravel(axes), path_L):
                    key_num = 0
                    h5keys = [item for item in h5_src[L][l][d].keys() if any(s in item for s in key_list)]
                    key_sorted = sorted(h5keys, key=natural_keys)
                    for key in key_sorted:
                        try:
                            nodeJ = h5_src[L][l][d][key]['J_avg']
                            nodeh = h5_src[L][l][d][key]['h_avg']
                        except Exception as er:
                            continue
                        key_num = key_num + 1
                        idx = get_v_filtered_index_list(h5_src[L][l][d][key], variance_window_limits[0])
                        data_logJ = np.log(np.array(nodeJ['data']))
                        data_logh = np.log(np.array(nodeh['data']))
                        dataD = data_logJ - data_logh
                        # hist, edges = np.histogram(dataD, bins=200, density=True)
                        num = len(dataD)
                        avg = np.mean(dataD)
                        length = nodeJ.attrs['chain_length']
                        delt = nodeJ.attrs['delta']
                        lamb = nodeJ.attrs['lambda']
                        nicename = re.sub(r'[\W_]', ' ', str(key))
                        nicename = nicename + ' (' + str(num) + ')'
                        hist, bincentres, dummy = ax.hist(x=dataD, bins=200, align='mid', color='b', linewidth=0.6, alpha=0.8, histtype='step', density=True,
                                                          label=nicename)
                        # Fit a normal distribution to the data:
                        mu, std = norm.fit(dataD)
                        mu_str = "{:.3f}".format(mu)
                        std_str = "{:.3f}".format(std)
                        # Plot the PDF.
                        p = norm.pdf(bincentres, mu, std)
                        ax.plot(bincentres, p, 'b', linewidth=1.4, alpha=0.8, label="N(" + mu_str + "," + std_str + ")")

                        # h5_loic = h5open('merged/LogNormal-mu-1-sigma-1.h5', 'r')
                        # data_logJ_loic = []
                        # data_logh_loic = []
                        # i = 0
                        # J_i = 0
                        # h_i = 0
                        # for rnd in h5_loic['Distrib'][0:4999]:
                        #     J_i = J_i + rnd
                        #     i = i + 1
                        #     if i == length:
                        #         data_logJ_loic.append(np.log(J_i/length))
                        #         J_i = 0
                        #         i = 0
                        # for rnd in h5_loic['Distrib'][5000:-1]:
                        #     h_i = h_i + rnd
                        #     i = i + 1
                        #     if i == length:
                        #         data_logh_loic.append(np.log(h_i / length))
                        #         h_i = 0
                        #         i = 0

                        # dataD_loic = np.asarray(data_logJ_loic) - np.asarray(data_logh_loic)
                        # nicename='Loic data (5000)'
                        # hist,bincentres, dummy = ax.hist(x=dataD_loic, bins=200, align='mid',color='g', linewidth=0.6,alpha=0.8, histtype='step', density=True,label=nicename)
                        # Fit a normal distribution to the data:
                        # mu, std = norm.fit(dataD_loic)
                        # mu_str = "{:.3f}".format(mu)
                        # std_str = "{:.3f}".format(std)
                        # Plot the PDF.
                        # p = norm.pdf(bincentres, mu, std)
                        # ax.plot(bincentres, p, 'g', linewidth=1.4, alpha=0.6, label="N(" + mu_str + "," + std_str + ") (Loic)")

                    ax.set_xlabel('$\log \\bar{J} - \log \\bar{h}$')
                    ax.set_ylabel('Histogram')
                    ax.set_title('$L = ' + str(length) + '$')
                    # ax.set_xlim(left=0)
                    # ax.set_yscale('log')
                    ax.legend()
                    used_ax = used_ax + 1
                fig.suptitle('Distribution of $\log \\bar{J} - \log \\bar{h}$ @ $\Delta = ' + str(delt) + '\quad \lambda = ' + str(lamb) + '$')
                for ax in np.ravel(axes)[used_ax:]:
                    fig.delaxes(ax)
                if plotdir != '':
                    Jh = re.sub('/', '_', str(d))
                    plt.savefig(plotdir + '/delta_distribution_' + l + '_' + str(Jh) + '.pdf', format='pdf')
                    plt.savefig(plotdir + '/delta_distribution_' + l + '_' + str(Jh) + '.png', format='png')
    except Exception as er:
        print("Could not plot Delta distribution. \nReason: " + str(er))
    h5close(h5_src)
