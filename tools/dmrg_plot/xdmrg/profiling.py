import h5py
import numpy
from src.io.h5ops import *
import matplotlib.pyplot as plt
from ast import literal_eval
from scipy.optimize import curve_fit


def line(x, k, m):
    return k * np.array(x) + m


def plot_timers(filename, label=None, filter=None, spin=2, mpod=25):
    x = []
    y = []
    h5file = h5open(filename, 'r')
    timers = h5file['xDMRG/state_0/tables/timers']
    names = np.char.decode(timers['name'][()].astype(np.bytes_), 'UTF-8')
    idxs = np.where(np.char.find(names, filter) != -1)[0]
    for idx in idxs:
        mps_dim = names[idx].split('mps-', 1)[1].split('_', 1)[0]
        mps_dim = np.array(literal_eval(mps_dim), dtype=int)
        mpo_dim = names[idx].split('mpo-', 1)[1].split('_', 1)[0]
        mpo_dim = np.array(literal_eval(mpo_dim), dtype=int)
        # mps_dims = list(map(int, mps_dims.split(',')))
        # print(type(mps_dim), mps_dim)
        # print(type(mpo_dim), mpo_dim)
        if mpo_dim[0] != mpod:
            continue
        if mps_dim[0] != spin:
            continue
        size = mps_dim[0] * mps_dim[1] * mps_dim[2]

        x.append(size)
        y.append(timers['avg'][idx])

    # plt.bar(x,y,width=10, log=True)
    order = np.argsort(x)
    x = np.log10(np.array(x)[order])
    y = np.log10(np.array(y)[order])
    # x = np.array(x)[order]
    # y = np.array(y)[order]

    plt.scatter(x, y, marker='o', s=2, label=label)

    popt, pcov = curve_fit(f=line, xdata=x, ydata=y)
    plt.plot(x, line(x, popt[0], popt[1]), marker=None, linewidth=0.8,
             linestyle='--', label='{} fit'.format(label))

    plt.legend()
    # plt.yscale('log')
    # plt.xscale('log')
    plt.ylabel('Time [s]')
    plt.xlabel('Size')
    plt.title('expval')


# plot_timers('profdata/expval2.h5', label='znver1 mpo=5', filter='expval_',mpod=5, spin=2)
# plot_timers('profdata/expval2.h5', label='znver1 mpo=25', filter='expval_',mpod=25, spin=2)
plot_timers('profdata/expval3.h5', label='native mpo=5', filter='expval_', mpod=5, spin=2)
plot_timers('profdata/expval3.h5', label='native mpo=25', filter='expval_', mpod=25, spin=2)
plot_timers('profdata/expval4.h5', label='native2 mpo=5', filter='expval_', mpod=5, spin=2)
plot_timers('profdata/expval4.h5', label='native2 mpo=25', filter='expval_', mpod=25, spin=2)
# plot_timers('data/lbit_prof3.h5', 'OpenBLAS 0.3.17', filter='o3.svd.lapacke.zgesdd-')
plt.show()
