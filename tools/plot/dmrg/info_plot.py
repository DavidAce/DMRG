import h5py
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress
from scipy.integrate import cumulative_trapezoid
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import seaborn as sns


def split_consecutive(data, minsize=None):
    consecutives = np.split(data, np.where(np.diff(data) != 1)[0] + 1)
    if len(consecutives) == 0 or minsize is None:
        return consecutives

    consecutives = [segment for segment in consecutives if len(segment) >= minsize]
    if len(consecutives) == 0:
        return [np.array([])]
    else:
        return consecutives

def get_increasing_points(data1d):
    # incr = np.where((np.diff(data1d, append=data1d[-1]) > 0) | (np.diff(data1d, prepend=data1d[0]) > 0))[0]
    # decr = np.where((np.diff(data1d, append=data1d[-1]) < 0) & (np.diff(data1d, prepend=data1d[0]) < 0))[0]

    decr = np.where(np.gradient(data1d) < 0)[0]
    incr =  np.where((np.diff(data1d, append=data1d[-1]) > 0) | (np.diff(data1d, prepend=data1d[0]) > 0) | (np.gradient(data1d) > 0))[0]
    return np.setdiff1d(incr, decr)

def get_decreasing_points(data1d):
    incr = np.where(np.gradient(data1d) > 0)[0]
    decr =  np.where((np.diff(data1d, append=data1d[-1]) < 0) | (np.diff(data1d, prepend=data1d[0]) < 0) | (np.gradient(data1d) < 0))[0]
    return np.setdiff1d(decr, incr)

def get_single_realization_masscenter(ips, debug=False, plot=False): # ips: info per scale  l x realizations
    if np.ndim(ips) == 1:
        ips = ips[:, np.newaxis]
    ldim, rdim = np.shape(ips)
    mc = np.nan * np.zeros(shape=(rdim,))
    xi =  np.nan * np.zeros(shape=(rdim,)) # Characteristic length scale of localization
    tau = np.nan * np.zeros(shape=(rdim,)) # Characteristic length scale of the topological bit
    bit = np.zeros(shape=(rdim,)) # Whether a topological bit was found (1 or 0)
    for r in range(rdim):
        fit_xi, fit_tau = None, None
        topological_bit_detected=False
        ips_cumsum = np.cumsum(ips[:,r])
        ips_nerror = np.abs(ldim - ips_cumsum[-1]) # Numerical error
        idx_exclude = np.where(ips[:, r] <= ips_nerror)[0]  # Consider only points above the numerical error


        if debug:
            print(f'{ips_cumsum[-1]=:.16f}')
            print(f'{ips_nerror=:.5e}')
            print(f'{ips[:, r]=}')

        # Detect the scales over there is an exponential change in information
        # Recall that xi can be either negative (localized) or positive (delocalized)

        idx_increase = get_increasing_points(ips[:, r])  # Points that are increasing
        idx_decrease = get_decreasing_points(ips[:, r])  # Points that are decreasing
        # idx_decrease = idx_decrease[np.where(idx_decrease <= ldim//2)] # Condisder decrease up to half system size
        idx_xi_decr = np.setdiff1d(idx_decrease, idx_exclude)  # Remove the scales where the information is below error
        idx_xi_incr = np.setdiff1d(idx_increase, idx_exclude)  # Remove the scales where the information is below error
        idx_xi_decr = split_consecutive(idx_xi_decr, minsize=3)[0]  # Keep the first consecutive set with more than 3 scales
        idx_xi_incr = split_consecutive(idx_xi_incr, minsize=3)[0]  # Keep the last consecutive set with more than 3 scales

        if debug:
            print('idx_xi_decr', get_decreasing_points(ips[:, r]), '-->', idx_xi_decr)
            print('idx_xi_incr', get_increasing_points(ips[:, r]), '-->', idx_xi_incr)

        # Determine whether the information grows or decays.
        # When localized, the information decays exp., but it may do so after first growing for a scale or two.
        # When thermal the information grows exp.
        # One way to determine is by comparing the center points of idx_xi_decr and idx_xi_incr.
        # Recall that we already know that both are at least 3 points long
        if len(idx_xi_incr) > 0 and len(idx_xi_decr) > 0:
            info_decays = np.mean(idx_xi_decr) < np.mean(idx_xi_incr)
        else:
            info_decays = len(idx_xi_decr) > len(idx_xi_incr)
        idx_xiexpfit = idx_xi_decr if info_decays else idx_xi_incr

        if len(idx_xiexpfit) >= 3:

            val_xiexpfit = ips[idx_xiexpfit, r]
            log_xiexpfit = np.ma.masked_where(val_xiexpfit <= 0,
                                              np.log(val_xiexpfit))  # Consider only the positive values
            slope, intercept, r_value, p_value, std_err = linregress(idx_xiexpfit, log_xiexpfit)  # Get the slope
            # We did a linear fit to y=kx + m = log(C*exp(i/xi))=i/xi + log(C), so the slope is 1/tau
            xi[r] = -1.0 / slope
            fit_xi = np.exp(idx_xiexpfit * slope + intercept)
            if debug:
                print('xi', xi[r])


        # Detect the scales where the last bit is located within numerical error
        idx_lastbit = np.where(ips_cumsum >= (ldim-1.0)+ips_nerror)[0] # The last bit location within numerical error
        if len(idx_lastbit) == 0: # This is an error, since there should be ldim bits in total
            print("WARNING: idx_lastbit is empty!")
            continue

        # The last bit isn't necessarily a topological bit.
        # We calculate the slope along the increasing scales of the last bit
        val_lastbit = ips[idx_lastbit,r] # The info values where the last bit is located
        idx_lastbit_incr = np.intersect1d(idx_lastbit, get_increasing_points(ips[:, r])) # Keep the increasing indices
        idx_lastbit_incr = np.setdiff1d(idx_lastbit_incr, idx_exclude)  # Remove the scales where the information is below error
        idx_lastbit_incr = split_consecutive(idx_lastbit_incr)[0]  # Keep the first consecutive set of indices
        if len(idx_lastbit_incr) > 1 and idx_lastbit_incr[0] > ldim/2.0: # Require that the bit is in the upper half
            val_lastbit_incr = ips[idx_lastbit_incr,r]   # The info values where the increasing part of the last bit is located
            log_lastbit_incr = np.ma.masked_where(val_lastbit_incr <= 0, np.log(val_lastbit_incr)) # Consider only the positive values
            slope, intercept, r_value, p_value, std_err = linregress(idx_lastbit_incr, log_lastbit_incr) # Get the slope
            topological_bit_detected = slope > 0 and np.sum(val_lastbit_incr) > 0.2 # Require at least 20% of the bit in the increasing region
            if topological_bit_detected:
                bit[r] = np.sum(val_lastbit) # Should be very close to 1
                # We did a linear fit to y=kx + m = log(C*exp(i/tau))=i/tau + log(C), so the slope is 1/tau
                tau[r] = 1.0/slope
                fit_tau = np.exp(idx_lastbit_incr*slope + intercept)
                if debug:
                    print(f'{idx_lastbit_incr=}')
                    print('tau', tau[r], 'bit', bit[r])

        num_interp = ldim*500
        fun_interp = interp1d(x=range(0,ldim),y=ips[:,r], kind='previous')
        idx_interp = np.linspace(0, ldim-1, num_interp, endpoint=True)
        sum_interp = cumulative_trapezoid(fun_interp(idx_interp), idx_interp)
        # sum_interp = np.cumsum(val_interp)  #cumulative_trapezoid(y=val_interp, x=idx_interp, initial=0) #np.interp(idx_interp, xp=range(ldim), fp=ips[:,r])


        # print(idx_interp)
        # print(sum_interp)
        # bits_to_count =  0.368*(ldim-bit[r])
        # frac_initial  = ips[0, r] / ldim
        frac_to_count =  0.5
        frac_bitcount =  frac_to_count*(ldim-bit[r])
        idx_masscenter = np.where(sum_interp >= frac_bitcount)[0][0]
        # print(idx_masscenter)
        # print(idx_interp[idx_masscenter])
        val_mass_center = idx_interp[idx_masscenter] * np.log(2)
        # val_mass_center = -idx_interp[idx_masscenter]  / np.log((1-frac_to_count)/frac_initial)

        mc[r] = val_mass_center
        if debug:
            print(f'mass center      : {val_mass_center:.16f}')
        # print('idx_interp', idx_interp)
        # print('val_interp', val_interp)
        if plot:
            if fit_xi is not None or fit_tau is not None:
                fig, ax = plt.subplots()
                ax.plot(range(ldim), ips[:, r])
                ax.axhline(y=ips_nerror, linestyle='--', label='precision')
                ax.axvline(x=val_mass_center, linestyle='--', label=f'mass center=${val_mass_center:.3f}$')
                ax.set_yscale('log')
            if fit_xi is not None:
                marker = 'v' if info_decays else '^'
                color = 'blue' if info_decays else 'red'
                label = 'decay' if info_decays else 'growth'
                ax.scatter(idx_xiexpfit, val_xiexpfit, marker=marker, color='black', s=30, facecolors='none', label=label)
                ax.plot(idx_xiexpfit, fit_xi, color=color, label='exp fit')
                ax.axvline(x=np.abs(xi[r]), linestyle=':', label=f'$\\xi={xi[r]:.3f}$')

                # tau
            if fit_tau is not None:
                ax.scatter(idx_lastbit, val_lastbit, marker='o', s=50, color='black', facecolors='none',label='last bit')
                ax.plot(idx_lastbit_incr, fit_tau, color='green',label=f'$\\tau={tau[r]:.3f}$')

            plt.legend(loc='best')
            plt.show()
    return mc, xi,tau,bit

def get_matching_row(table, row):
    sourcerowid = row[['iter', 'step', 'position', 'event', 'bond_lim']]
    for idx, tablerowid in enumerate(table['iter', 'step', 'position', 'event', 'bond_lim']):
        if sourcerowid == tablerowid:
            return table[idx][()]

def plot_info_per_scale(ax, file, state, eventnames=None):
    if eventnames is None:
        eventnames = []
    if not isinstance(eventnames, list):
        raise TypeError('expected event:list')

    with h5py.File(file, 'r') as f:
        ax.set_xlabel('$\ell$')
        ax.set_ylabel('Information (bits)')
        if not state in f:
            return
        # scales = [f'scale{l}' for l in range(L)]
        infodset=f'{state}/information_per_scale'
        icomdset=f'{state}/information_center_of_mass'
        statusdset=f'{state}/status'
        msmntsdset=f'{state}/measurements'
        modeldset=f"{state.split('/')[0]}/model/hamiltonian"
        scalekeys = [key for key in f[infodset].dtype.fields.keys() if 'scale' in key]
        ldim = len(scalekeys)
        rows = range(len(f[infodset]))
        info_enum_event = h5py.check_enum_dtype(f[f'{infodset}'].dtype['event'])  # key value pairs defining the enum
        eventnumbers = [ val for key,val in info_enum_event.items() if key in eventnames]
        print(state)
        print('This is the model', modeldset)
        print(f)
        print(f[f'{modeldset}'])
        dval =f[f'{modeldset}'][0]['delta']
        gval = f[f'{modeldset}'][0]['g']
        imin = 1.0
        imax = 1.0
        for row in reversed(rows):
            if not f[infodset]['event'][row] in eventnumbers:
                print(f[infodset]['event'][row], eventnumbers)
                continue
            infoperscale = np.abs(np.array(list(f[infodset][row][scalekeys])))
            for i in range(len(infoperscale)):
                if infoperscale[i] <= 0:
                    infoperscale[i] = 0.5 * (infoperscale[i-1] + infoperscale[i+1]) # Take the average of neighbors if zero
            # Get the smallest positive value
            # infoperscale_mask = np.ma.masked_where(infoperscale, infoperscale < 1e-16)
            # infoperscale = np.abs(infoperscale)
            # infoperscale[infoperscale <= 0] = np.nan
            icom = f[icomdset]['scale'][row]
            # mc, xi,tau,bit = get_single_realization_masscenter(infoperscale)
            bond_lim = f[infodset][row]['bond_lim']
            if bond_lim <= 1:
                continue
            imin = np.nanmin([imin, np.nanmin(np.abs(infoperscale))]) if bond_lim > 1 else imin
            imax = np.nanmax([imax, np.nanmax(np.abs(infoperscale))]) if bond_lim > 1 else imax
            # measurements = get_matching_row(f[msmntsdset], f[infodset][row])
            measurements = f[msmntsdset][row]
            energy = measurements['energy']
            variance = measurements['energy_variance']
            # label = f'$\chi={bond_lim}$ $E={energy:.5f}$ $\mathrm{{Var}}(H)=${variance:.2e}'
            label = f'$\Delta={dval}$'
            p = ax.plot(range(ldim), infoperscale, label=label)
            # ax.axvline(x=icom,ymin=1-(0.1 - 0.005*row), ymax=1.0, color=p[-1].get_color())
    ax.legend(loc='upper right',ncol=2)
    ax.set_ylim(ymin=0.5*imin, ymax=2*imax)
    ax.set_yscale('log')


def plot_tabledata(ax, file, state, xaxis=None, yaxis=None, yscale=None, eventnames=None):
    if eventnames is None:
        eventnames = []
    if not isinstance(eventnames, list):
        raise TypeError('expected event:list')
    if not isinstance(xaxis, tuple):
        raise AssertionError("Expected xaxis == tuple(table column, axis label) e.g. tuple(bond_lim, $\chi$)")
    if not isinstance(yaxis, tuple):
        raise AssertionError("Expected yaxis == tuple(table column, axis label) e.g. tuple(energy, $E$)")
    with h5py.File(file, 'r') as f:
        ax.set_xlabel(xaxis[1])
        ax.set_ylabel(yaxis[1])
        if not state in f:
            return
        # scales = [f'scale{l}' for l in range(L)]
        infodset=f'{state}/information_per_scale'
        statusdset=f'{state}/status'
        msmntsdset=f'{state}/measurements'
        modeldset=f"{state.split('/')[0]}/model/hamiltonian"
        scalekeys = [key for key in f[infodset].dtype.fields.keys() if 'scale' in key]
        ldim = len(scalekeys)
        rows = range(len(f[infodset]))
        info_enum_event = h5py.check_enum_dtype(f[f'{infodset}'].dtype['event'])  # key value pairs defining the enum
        eventnumbers = [val for key, val in info_enum_event.items() if key in eventnames]

        for row in reversed(rows):
            if not f[infodset]['event'][row] in eventnumbers:
                # print(f[infodset]['event'][row], eventnumbers)
                continue
            # infoperscale = np.array(list(f[infodset][row][scalekeys]))
            # mc, xi,tau,bit = get_single_realization_masscenter(infoperscale)
            if xaxis[0] in f[msmntsdset].dtype.fields.keys():
                xdata = f[msmntsdset][row][xaxis[0]]
            if xaxis[0] in f[statusdset].dtype.fields.keys():
                xdata = f[statusdset][row][xaxis[0]]

            # measurements = get_matching_row(f[msmntsdset], f[infodset][row])
            measurements=f[msmntsdset][row]
            ydata = measurements[yaxis[0]]
            bond_lim = f[infodset][row]['bond_lim']
            if bond_lim < 2:
                continue
            # variance = measurements['energy_variance']
            # label = f'$\chi={bond_lim}$ $E={energy:.5f}$ $\mathrm{{Var}}(H)=${variance:.2e}'
            # print(np.sum(infoperscale))
            ax.scatter([xdata], [ydata], label=None)
            # plt.axvline(x=mc,ymin=0, ymax=0.1 + 0.01*row, color=p[-1].get_color())
        if yscale is not None:
            ax.set_yscale(yscale)
def plot_variance_per_scale(ax, file, state):
    with h5py.File(file, 'r') as f:
        ax.set_xlabel('$\chi$')
        ax.set_ylabel('$\mathrm{Var}(H)$')
        ax.set_yscale('log')
        # scales = [f'scale{l}' for l in range(L)]
        if not state in f:
            return
        infodset=f'{state}/information_per_scale'
        statusdset=f'{state}/status'
        msmntsdset=f'{state}/measurements'
        modeldset=f'{state}/../model/hamiltonian'
        scalekeys = [key for key in f[infodset].dtype.fields.keys() if 'scale' in key]
        ldim = len(scalekeys)
        rows = len(f[infodset])
        for row in range(rows):
            infoperscale = np.array(list(f[infodset][row][scalekeys]))
            mc, xi,tau,bit = get_single_realization_masscenter(infoperscale)
            bond_lim = f[infodset][row]['bond_lim']
            if bond_lim < 2:
                continue
            measurements = get_matching_row(f[msmntsdset], f[infodset][row])
            energy = measurements['energy']
            variance = measurements['energy_variance']
            print(np.sum(infoperscale))
            ax.scatter([bond_lim], [variance], label=None)
            # plt.axvline(x=mc,ymin=0, ymax=0.1 + 0.01*row, color=p[-1].get_color())


plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": "cm",
    "mathtext.fontset": "cm",
    'figure.constrained_layout.w_pad':  0.0,   # inches
    'figure.constrained_layout.h_pad':  0.0,   # inches
    'figure.subplot.top'   : 0.97,
    'figure.subplot.bottom': 0.06,
    'figure.subplot.left': 0.06,
    'figure.subplot.right' : 0.98,
    'figure.subplot.wspace' : 0.15,
    'figure.subplot.hspace' : 0.5,
    'savefig.bbox'           : 'standard',
    'savefig.pad_inches'     : 0.5
})
plt.rc('text.latex', preamble=r'\usepackage{amssymb}')

px = 1/plt.rcParams['figure.dpi']  # pixel in inches


sns.set_palette("viridis",  1)
fig1, axes1 = plt.subplot_mosaic('''
                                II
                                EV
                                ''',
                              figsize=(1152*px, 864*px),
                              layout="tight",
                              )
# fig.subplots_adjust(wspace=0.20,right=0.8)
# dir='/home/david/GitProjects/DMRG++/output'
dir ='/mnt/WDB-AN1500/mbl_transition/fdmrg8-trf/output'
# file='mbl_11000010000182-unfinished.h5'
# file='mbl_20000030000012-finished.h5'
# file='mbl_11000010000182-infovarsat-trncreduce.h5'
# file='mbl_11000010000182-unfinished-info-bondreduce.h5'
delta='+5.00'
size='64'
seed='25000500000000'
file=f'L{size}/g0.500/d{delta}/mbl_{seed}.h5'
state="fDMRG/state_emin"

plot_info_per_scale(ax=axes1['I'], file=f'{dir}/{file}', state=state, eventnames=['FINISHED'])
plot_tabledata(ax=axes1['E'],   file=f'{dir}/{file}', state=state, xaxis=('algorithm_time', '$t$'), yaxis=('energy', '$E$'), eventnames=['FINISHED'])
plot_tabledata(ax=axes1['V'], file=f'{dir}/{file}', state=state, xaxis=('algorithm_time', '$t$'),yaxis=('energy_variance', '$\mathrm{Var}(H)$'), eventnames=['FINISHED'])
# plot_info_per_scale(file=f'{dir}/mbl_11000010000181.h5', state='xDMRG/state_emid')
# plot_info_per_scale(file=f'{dir}/mbl_11000010000181-chi40.h5', state='xDMRG/state_emid')
# plot_info_per_scale(file=f'{dir}/mbl_11000010000181-chi48.h5', state='xDMRG/state_emid')
# plot_info_per_scale(file=f'{dir}/mbl_11000010000181-chi64.h5', state='xDMRG/state_emid')
fig1.suptitle(f'fDMRG convergence \n $L={size}$ $g=0.5$ $\Delta={delta}$ rnd.seed={seed}')

sns.set_palette("viridis_r",  10)
fig2, axes2 = plt.subplot_mosaic('''
                                II
                                EV
                                ''',
                              figsize=(1152*px, 864*px),
                              layout="tight",
                              )
plot_info_per_scale(ax=axes2['I'], file=f'{dir}/{file}', state=f'{state}/rbds', eventnames=['RBDS_STEP'])
plot_tabledata(ax=axes2['E'],   file=f'{dir}/{file}', state=f'{state}/rbds', xaxis=('bond_lim', '$\chi$'), yaxis=('energy', '$E$'), eventnames=['RBDS_STEP'])
plot_tabledata(ax=axes2['V'], file=f'{dir}/{file}', state=f'{state}/rbds', xaxis=('bond_lim', '$\chi$'),yaxis=('energy_variance', '$\mathrm{Var}(H)$'), eventnames=['RBDS_STEP'],yscale='log')
fig2.suptitle(f'Reverse bond dimension scaling \n $L={size}$ $g=0.5$ $\Delta={delta}$ rnd.seed={seed}')

sns.set_palette("viridis_r",  10)
fig3, axes3 = plt.subplot_mosaic('''
                                II
                                EV
                                ''',
                              figsize=(1152*px, 864*px),
                              layout="tight",
                              )
plot_info_per_scale(ax=axes3['I'], file=f'{dir}/{file}', state=f'{state}/rtes', eventnames=['RTES_STEP'])
plot_tabledata(ax=axes3['E'],   file=f'{dir}/{file}', state=f'{state}/rtes', xaxis=('trnc_lim', '$\epsilon$'), yaxis=('energy', '$E$'), eventnames=['RTES_STEP'])
plot_tabledata(ax=axes3['V'], file=f'{dir}/{file}', state=f'{state}/rtes', xaxis=('trnc_lim', '$\epsilon$'),yaxis=('energy_variance', '$\mathrm{Var}(H)$'), eventnames=['RTES_STEP'],yscale='log')
fig3.suptitle(f'Reverse truncation error scaling\n $L={size}$ $g=0.5$ $\Delta={delta}$ rnd.seed={seed}')


sns.set_palette("viridis",  5)
fig4, axes4 = plt.subplot_mosaic('''
                                II
                                ''',
                              figsize=(1152*px, 864*px),
                              layout="tight",
                              )
fig4.suptitle(f'Information lattice for increasing delta\n $L={size}$ $g=0.5$')

points = [
    {'L': '64', 'g': '0.500', 'd': '+3.00', 'seed': 23000500000004},
    {'L': '64', 'g': '0.500', 'd': '+4.00', 'seed': 24000500000002},
    # {'L': '64', 'g': '0.500', 'd': '+4.00', 'seed': 24000500000008},
    {'L': '64', 'g': '0.500', 'd': '+5.00', 'seed': 25000500000005},
    {'L': '64', 'g': '0.500', 'd': '+6.00', 'seed': 26000500000000},
    {'L': '64', 'g': '0.500', 'd': '+9.00', 'seed': 29000500000000},
]

for p in points:
    file=f"L{p['L']}/g{p['g']}/d{p['d']}/mbl_{p['seed']}.h5"
    plot_info_per_scale(ax=axes4['I'], file=f'{dir}/{file}', state=state, eventnames=['FINISHED'])
    # plot_tabledata(ax=axes3['E'],   file=f'{dir}/{file}', state=f'{state}/rtes', xaxis=('trnc_lim', '$\epsilon$'), yaxis=('energy', '$E$'), eventnames=['RTES_STEP'])
    # plot_tabledata(ax=axes3['V'], file=f'{dir}/{file}', state=f'{state}/rtes', xaxis=('trnc_lim', '$\epsilon$'),yaxis=('energy_variance', '$\mathrm{Var}(H)$'), eventnames=['RTES_STEP'],yscale='log')





# plt.ylabel('Information')
# plt.xlabel('$\ell$')
plt.show()
