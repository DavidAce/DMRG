from src.io.h5ops import *
from src.plotting.tools import *
from matplotlib import pyplot as plt
import os
import glob


def find_saturation_idx1(ydata, ax=None):
    # Consider Y vs X: a noisy signal decaying in the shape of a hockey-club, say.
    # We want to identify the point at which the signal stabilizes. We use the fact that the
    # standard deviation is high if it includes parts of the non-stable signal, and low if
    # it includes only the stable part.
    # Here we monitor the standard deviation of the signal between [start_point, end_point],
    # and move "start_point" towards the end. If the standard deviation goes below a certain
    # threshold, i.e. threshold < max_std, then we have found the stabilization point.
    sdata = []
    for i in range(len(ydata)):
        sdata.append(np.std(ydata[i:]))
    sdiff = -np.log10(np.abs(np.diff(sdata)))
    idx = np.argmax(sdiff)
    if ax:
        ax.plot(sdiff)
        ax.plot(idx, sdiff[idx], 'ro')
    return np.argmax(sdiff)


def get_slope_interval(ydata, start=0, end=-1):
    ydata = ydata[start:end]
    if len(ydata) <= 1:
        return np.nan, np.nan

    xdata = range(1, len(ydata) + 1)
    xmean = np.mean(xdata)
    ymean = np.mean(ydata)
    sxx = np.sum((xdata - xmean) ** 2)
    syy = np.sum((ydata - ymean) ** 2)
    sxy = np.sum((xdata - xmean) * (ydata - ymean))
    slope = sxy / sxx
    residual = syy * (1 - sxy * sxy / sxx / syy)
    return slope, residual


def get_slope_at(ydata, at, w=1):
    min_idx = np.max([at - w, 0])
    max_idx = np.min([at + w, len(ydata) - 1])
    return get_slope_interval(ydata, min_idx, max_idx)


def smoother(ydata, w, reps=1):
    ytmp = ydata
    for r in range(reps):
        sdata = []
        for i, y in enumerate(ydata):
            min_idx = np.max([i - w, 0])
            max_idx = np.min([i + w, len(ydata) - 1])
            sdata.append(np.mean(ytmp[min_idx:max_idx]))
        ytmp = sdata
    return ytmp


def find_saturation_idx5(ydata, ax=None, xstart=None, xstop=None, tag=''):
    if len(ydata) <= 2:
        return
    if not xstart:
        xstart = 0
    if not xstop:
        xstop = len(ydata) - 1
    # ysmth = smoother(ydata,w=4)
    # ylog = -np.log10(ydata)
    # ylog = ylog/ylog[-1]
    # ysmt = smoother(ydata,w=2,reps=1)
    yavg = []
    for i, s in enumerate(ydata):
        yavg.append(np.mean(ydata[i:]))

    # std_data = []
    # ste_data = []
    sta_data = []
    slp_data = []
    avg_data = []
    w = 4
    for i, s in enumerate(yavg):
        min_idx = np.min([len(yavg) - w, i])
        min_idx = np.max([min_idx, 0])

        points = len(yavg[min_idx:])
        # std_data.append(np.std(ylog[min_idx:]))
        # ste_data.append(std_data[-1]/np.sqrt(points))
        sta_data.append(np.std(yavg[min_idx:]))
        # avg_data.append(np.mean(ydata[min_idx:]))
        if points <= 2:
            slp_data.append(slp_data[-1])
        else:
            avg = np.mean(yavg[min_idx:])
            itr = len(yavg[min_idx:])
            slp, res = get_slope_interval(yavg, min_idx, )
            slp_data.append(np.abs(slp / avg))
    ste_thresh = 1e-2
    slp_thresh = 1e-2
    sta_thresh = 1e-2
    for i, s in enumerate(sta_data):
        idx = i
        # if s < ste_thresh and slp_data[i] < slp_thresh:
        #     break
        # if s < 0.1*ste_thresh:
        #     break
        # if slp_data[i] < 0.1*slp_thresh:
        #     break
        med = np.median(yavg[i:])
        wdt = 10 * sta_data[i]
        if sta_data[i] / yavg[i] < sta_thresh and ydata[i] <= med + wdt and ydata[i] >= med - wdt:
            break

    if ax:
        # ax.plot(ylog, label='l', marker='.')
        ax.plot(ydata, label='$\mathrm{Var}H$ ' + tag, marker='.')
        # ax.plot(ylog, label='ylog',marker='.')
        # ax.plot(ysmt, label='ysmt',marker='o')
        ax.plot(yavg, label='yavg', marker='v')
        ax.plot(sta_data, label='ysta', marker='<')
        # ax.plot(std_data, label='Std[$-\log{(\mathrm{Var}H)}]$' + tag,marker='.')
        ax.plot(slp_data, label='slp', marker='x')

        ax.plot(idx, ydata[idx], 'ro')
        ax.axvline(idx, color='grey', marker=None, linestyle=':', alpha=0.5)
        ax.axhline(ydata[idx], color='grey', marker=None, linestyle=':', alpha=0.5)
        # ax.axhline(ste_data[idx],color='grey',marker=None,linestyle=':',alpha=0.5)
        # ax.axhline(slp_data[idx],color='green')
        # ax.plot(idx, thresh, 'bo')
        ax.set_yscale('log')
        ax.set_xlim([xstart, xstop])
    return 0


# probably, there are some repetitions of iter in the hdf5 table

# h5file = h5open('convdata/mbl_1004072.h5', 'r')
# h5file = h5open('convdata/mbl_10004_iter32.h5', 'r')


db = {}
algo = 'xDMRG'
state = 'state_0'
# algo = 'fDMRG'
# state = 'state_emin'
h5patterns = ['convdata/mbl_5011228-20*',
              'convdata/mbl_5011228-19*',
              'convdata/mbl_5011228-18*',
              'convdata/mbl_5011228-17*',
              # 'convdata/mbl_5011228-16*',
              'convdata/mbl_5011228-15*',
              # 'convdata/mbl_5011228-14*',
              # 'convdata/mbl_5011228-13*',
              'convdata/mbl_5011228-12*.h5',
              # 'convdata/mbl_5011228-11*.h5',
              'convdata/mbl_5011228*.h5']
h5patterns = ['convdata/mbl_5011228.h5']

for pattern in h5patterns:
    path = glob.glob(pattern)[0]
    h5file = h5open(path, 'r')
    db[path] = {}
    db[path]['file'] = h5file
    db[path]['path'] = path
    db[path]['name'] = os.path.basename(path)
    db[path]['size'] = h5file['{}/{}/tables'.format(algo, state)].get('measurements')['length'][-1]
    db[path]['variance'] = h5file['{}/{}/tables'.format(algo, state)].get('measurements')['energy_variance'][()]
    db[path]['entropy'] = h5file['{}/{}/tables'.format(algo, state)].get('measurements')['entanglement_entropy_midchain'][()]
    db[path]['iter'] = h5file['{}/{}/tables'.format(algo, state)].get('measurements')['iter'][()]

    db[path]['entropies'] = {
        'iter': [],
        'data': []
    }
    for iterkey, iternode in h5file['{}/{}/checkpoint'.format(algo, state)].items():
        if not "iter" in iterkey:
            continue
        if "iter_last" in iterkey:
            continue
        db[path]['entropies']['data'].append(iternode['entanglement_entropies'][()])
        db[path]['entropies']['iter'].append(iternode.get('status')['iter'][-1])

print(db.keys())
start = 0  # 420
stop = 150
numplots = 12
rows, cols = get_optimal_subplot_num(numplots)
fig, axes = plt.subplots(nrows=rows, ncols=cols, figsize=(3 * cols, 3 * rows))
fig.tight_layout(pad=5, w_pad=1.0, h_pad=1.0)
fig.subplots_adjust(wspace=0.3, hspace=0.3)
for file in db.values():
    if not stop:
        stop = file['iter'][-1]
    if numplots == 1:
        iters = [stop]
    else:
        iters = np.unique(np.linspace(start, stop, numplots, dtype=int, endpoint=True))
    for i, ax in zip(iters, np.ravel(axes)):
        xrange = range(start, i)
        idx = find_saturation_idx5(file['variance'][0:i], ax, start, i, file['name'])
        ax.set_xlabel('iteration')
        ax.legend(framealpha=0.7, fontsize='x-small', labelspacing=0.25)
    fig.suptitle('$\mathrm{Var}(H)$ vs iter')

plt.show()
exit()

start = 400  # 420
numplots = 1
rows, cols = get_optimal_subplot_num(numplots)
fig, axes = plt.subplots(nrows=rows, ncols=cols)  # , figsize=(3 * cols, 3 * rows))
fig.tight_layout(pad=5, w_pad=1.0, h_pad=1.0)
fig.subplots_adjust(wspace=0.3, hspace=0.3)
ymin = None
ymax = None
for file in db.values():
    stop = file['iter'][-1]
    if numplots == 1:
        iters = [stop]
    else:
        iters = np.unique(np.linspace(start, stop, numplots, dtype=int, endpoint=True))
    for i, ax in zip(iters, np.ravel(axes)):
        # idx = find_saturation_idx5(entropy[0:i], ax)
        if (start < i):
            ax.plot(file['entropy'][0:i], label='S', marker='.')
            ax.set_xlim([start, i])
            ymin_tmp = np.min(file['entropy'][start:i]) if start != i else file['entropy'][0]
            ymax_tmp = np.max(file['entropy'][start:i]) if start != i else file['entropy'][i]
            ymin = np.min([ymin, ymin_tmp]) if ymin else ymin_tmp
            ymax = np.max([ymax, ymax_tmp]) if ymax else ymax_tmp
            ax.set_ylim([0.999 * ymin, 1.001 * ymax])
            ax.set_ylabel('S')
            ax.set_xlabel('iter')
            # ax.plot(variance)
            # ax.plot(idx, variance[idx], 'ro')
            # ax.set_yscale('log')
            ax.legend(framealpha=0.7, fontsize='x-small', labelspacing=0.25)
    fig.suptitle('Midchain entropy vs iter')

length = list(db.values())[0]['size']
print(length)
rows, cols = get_optimal_subplot_num(length + 1)
fig, axes = plt.subplots(nrows=rows, ncols=cols, figsize=(3.5 * cols, 3.5 * rows))
fig.tight_layout(pad=5, w_pad=1.0, h_pad=1.0)
fig.subplots_adjust(wspace=0.3, hspace=0.3)

for file in db.values():
    iter = file['entropies']['iter']
    if len(iter) > 2:
        data = np.asarray(file['entropies']['data'])
        size = np.shape(data)[1]
        for i, ax in enumerate(np.ravel(axes)):
            print(i, ax)
            print(np.shape(data))
            if i >= size:
                break
            ax.plot(iter, data[:, i], label=file['name'], marker='.')
            ax.set_ylabel('S')
            ax.set_xlabel('iter')
            ax.set_title('Site ' + str(i))
            ax.legend(framealpha=0.7, fontsize='x-small', labelspacing=0.25)
#
# if len(entropies) > 2:
#     entropies = np.asarray(entropies)
#     iteration = np.asarray(iteration)
#     rows,cols = get_optimal_subplot_num(np.shape(entropies)[1])
#     fig,axes = plt.subplots(nrows=rows, ncols=cols, figsize=(3.5* cols, 3.5* rows))
#     fig.tight_layout(pad=5, w_pad=1.0, h_pad=1.0)
#     fig.subplots_adjust(wspace=0.3, hspace=0.3)
#     for i,ax in enumerate(np.ravel(axes)):
#         print(i,ax)
#         print(np.shape(entropies))
#         if i >= np.shape(entropies)[1]:
#             break
#         ax.plot(iteration,entropies[:,i], label=None, marker='.')
#         ax.set_ylabel('S')
#         ax.set_xlabel('iter')
#         ax.set_title('Site ' + str(i))
#
#     fig.suptitle('Entanglement entropy on all sites')


# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# (x, y) = np.meshgrid(np.arange(entropies.shape[1]), np.arange(entropies.shape[0]))
# print(np.shape(x), np.shape(y))
# surf = ax.plot_surface(x, y, entropies)

plt.show()
