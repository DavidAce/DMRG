import h5py
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

h5files_old = [
    # '../../output/pert-L12-g1e-5.h5',
    # ('../../output/pert-L14-g1e-5.h5', '.', 'gray', '$L=14$'),
    ('../../output/pert-L48-g1e-6.h5', 's', 'red'),
    # ('../../output/pert-L32-g1e-6.h5', 'v', 'green'),
    # ('../../output/pert-L16-g1e-6.h5', 'o', 'blue'),

]

h5files_new = [
    # '../../output/pert-L12-g1e-5.h5',
    # ('../../output/pert-L14-g1e-5.h5', '.', 'gray', '$L=14$'),
    # ('../../output/pert-L48-g1e-6.h5', 's', 'red'),

    # ('../../output/pert-L16-g1e-6.h5', 'o', 'blue'),
    ('../../output/pert-L16-new.h5', 's'),
    # ('../../output/pert-L32-new.h5', 'o'),
    ('../../output/pert-L48-new.h5', 'v'),

]

# fig, axes = plt.subplots(nrows=2, ncols=2)
px = 1 / plt.rcParams['figure.dpi']  # pixel in inches
fig1, axes1 = plt.subplot_mosaic('''
                                AB
                                CD
                                ''',
                               figsize=(1152 * px, 864 * px),
                               layout="tight",
                               )
fig2, axes2 = plt.subplot_mosaic('''
                                AB
                                CD
                                ''',
                               figsize=(1152 * px, 864 * px),
                               layout="tight",
                               )

# deltas= [0.0, 0.25, 0.5, 0.75, 1.0,1.25, 1.5, 2.0]


for h5file, marker, color in h5files_old:
    xvals = []
    nvals = []
    Lvals = []
    yA1vals, eA1vals = [], []
    yB1vals, eB1vals = [], []
    xC1vals, yC1vals, eC1vals = [], [], []
    yA2vals, eA2vals = [], []
    yB2vals, eB2vals = [], []
    xC2vals, yC2vals, eC2vals = [], [], []
    with h5py.File(h5file, 'r') as f:
        for key, val in f.items():
            g_pert = val['g_pert'][0]
            L = val['model_size'][0]
            num = len(val)
            delta = val['delta'][0]
            xvals.append(delta)
            nvals.append(num)
            Lvals.append(L)

            ene_max = np.ones(num)
            if 'energy_upper_bound' in val.dtype.names:
                ene_max = val['energy_change'][()]

                # print(ene_max)
            # else:
            #     continue
            # print(ene_max)
            var_num = len(val['variance'][()])
            var_old = np.abs(val['variance'][()])
            var_prt = np.abs(val['variance_pert'][()])
            var_abs_chn = np.median(np.abs(var_prt - var_old))
            var_rel_chn = np.median(np.abs(var_prt - var_old) / var_prt)
            var_frc_chn = np.median(var_prt / var_old)
            var_abs_chn_ste = np.std(np.abs(var_prt - var_old)) / np.sqrt(var_num)
            var_rel_chn_ste = np.std(np.abs(var_prt - var_old) / var_old) / np.sqrt(var_num)
            var_frc_chn_ste = np.std(var_prt / var_old) / np.sqrt(var_num)

            yA1vals.append(var_abs_chn/L)
            eA1vals.append(var_abs_chn_ste/L)
            yA1titl = f'Absolute change/L in energy variance $g = 0\\rightarrow 10^{{{int(np.log10(g_pert))}}}$'
            yA1labl = "${\\mathrm{median}(\sigma_{H}^2 - \sigma_{H^\prime}^2)}/{L}$"

            yB1vals.append(var_rel_chn/L)
            eB1vals.append(var_rel_chn_ste/L)
            yB1titl = f'Relative change/L in energy variance $g = 0\\rightarrow 10^{{{int(np.log10(g_pert))}}}$'
            yB1labl = "${\\mathrm{median}(|\sigma_{H}^2 - \sigma_{H^\prime}^2|/{\sigma_{H}^2}) L}$"


            yC1vals.append(var_frc_chn/L)
            eC1vals.append(var_frc_chn_ste/L)
            yC1titl = f'Ratio change/L in energy variance $g = 0\\rightarrow 10^{{{int(np.log10(g_pert))}}}$'
            yC1labl = "${\\mathrm{median}(\sigma_{H^\prime}^2 /\sigma_{H}^2 )}$"


            ene_num = len(val['energy'][()])
            ene_old = val['energy'][()]
            ene_prt = val['energy_pert'][()]
            enech_avg = np.median(np.abs(ene_old - ene_prt))
            ene_abs_chn = np.median(np.abs(ene_prt - ene_old))
            ene_rel_chn = np.median(np.abs(ene_prt - ene_old) / np.abs(ene_old))
            ene_rel_ste = np.std(np.abs(ene_old - ene_prt)) / np.sqrt(ene_num)
            ene_abs_ste = np.std(np.abs(ene_prt - ene_old)) / np.sqrt(ene_num)

            yA2vals.append(ene_abs_chn/L)
            eA2vals.append(ene_abs_ste/L)
            yA2titl = f'Absolute change/L in energy $g = 0\\rightarrow 10^{{{int(np.log10(g_pert))}}}$'
            yA2labl = "${\\mathrm{median}(E - \sigma_{H^\prime}^2)}/{\sigma_{H}^2 L}$"


            yB2vals.append(ene_rel_chn)
            eB2vals.append(ene_rel_ste)

#  H -> H + gH'
#  H' = XX + ZZ ~ 2L
#  H²-> H² + gHH' + O(g²)
# gH' ~ 4L

            # if 'hsquared' in val.dtype.names:
            #     hsq_num = len(val['hsquared'][()])
            #     hsq_old = val['hsquared'][()]
            #     hsq_prt = val['hsquared_pert'][()]
            #     hsqch_avg = np.median(np.abs(hsq_old - hsq_prt))
            #     hsq_abs_chn = np.median(np.abs(hsq_prt - hsq_old))
            #     hsq_rel_chn = np.median(np.abs(hsq_prt - hsq_old) / np.abs(hsq_old))
            #     hsq_frc_chn = np.median(hsq_prt / hsq_old)
            #     hsq_rel_ste = np.std(np.abs(hsq_prt - hsq_old) / np.abs(hsq_old)) / np.sqrt(hsq_num)
            #     hsq_abs_ste = np.std(np.abs(hsq_prt - hsq_old)) / np.sqrt(hsq_num)
            #     hsq_frc_ste = np.std(np.abs(hsq_prt / hsq_old)) / np.sqrt(hsq_num)
            #
            #     yCvals.append(hsq_rel_chn/L)
            #     eCvals.append(hsq_rel_ste/L)
            #     xCvals.append(val['delta'][0])
            # print(f'{key}: {y1vals[-1]:.3e} +- {varch_ste:.3e}')

    sort = np.argsort(np.asarray(xvals))
    xvals = np.asarray(xvals)[sort]
    yA1vals = np.asarray(yA1vals)[sort]
    eA1vals = np.asarray(eA1vals)[sort]
    yB1vals = np.asarray(yB1vals)[sort]
    eB1vals = np.asarray(eB1vals)[sort]

    yA2vals = np.asarray(yA2vals)[sort]
    eA2vals = np.asarray(eA2vals)[sort]
    yB2vals = np.asarray(yB2vals)[sort]
    eB2vals = np.asarray(eB2vals)[sort]

    # if(len(yCvals) > 0):
    #     yCvals = np.asarray(yCvals)[sort]
    #     eCvals = np.asarray(eCvals)[sort]
    #     xCvals = np.asarray(xCvals)[sort]

    label = f'$L={{{np.min(Lvals)}}}$ ($n_\min={{{np.min(nvals)}}}$)'

    lineA1 = axes1['A'].errorbar(x=xvals, y=yA1vals, yerr=None, color=color, marker=marker)
    lineA1.set_label(label)
    axes1['A'].axvline(x=np.log(2))
    axes1['A'].axvline(x=-np.log(2))
    axes1['A'].set_xlabel('$\delta$')
    axes1['A'].set_ylabel(yA1labl)
    axes1['A'].set_title(yA1titl)
    # axes1['A'].set_yscale('log')
    axes1['A'].legend(loc='best')

    lineB1 = axes1['B'].errorbar(x=xvals, y=yB1vals, yerr=None, color=color, marker=marker)
    lineB1.set_label(label)
    axes1['B'].axvline(x=np.log(2))
    axes1['B'].axvline(x=-np.log(2))
    axes1['B'].set_xlabel('$\delta$')
    axes1['B'].set_ylabel(yB1labl)
    axes1['B'].set_title(yB1titl)
    axes1['B'].legend(loc='best')

    lineC1 = axes1['C'].errorbar(x=xvals, y=yC1vals, yerr=None, color=color, marker=marker)
    lineC1.set_label(label)
    axes1['C'].axvline(x=np.log(2))
    axes1['C'].axvline(x=-np.log(2))
    axes1['C'].set_xlabel('$\delta$')
    axes1['C'].set_ylabel(yC1labl)
    axes1['C'].set_title(yC1titl)
    axes1['C'].legend(loc='best')




    lineA2 = axes2['A'].errorbar(x=xvals, y=yA2vals, yerr=None, color=color, marker=marker)
    lineA2.set_label(label)
    axes2['A'].axvline(x=np.log(2))
    axes2['A'].axvline(x=-np.log(2))
    axes2['A'].set_xlabel('$\delta$')
    axes2['A'].set_ylabel("$\\mathrm{median}(\\mathrm{Var}(H^\prime)) - \min$")
    # axes2['A'].set_title(f'Change in energy variance when $g = 0\\rightarrow 10^{{{int(np.log10(g_pert))}}}$')
    axes2['A'].set_title(yA2labl)
    # axes2['A'].set_yscale('log')
    axes2['A'].legend(loc='best')

    lineB2 = axes2['B'].errorbar(x=xvals, y=yB2vals, yerr=None, color=color, marker=marker)
    lineB2.set_label(label)
    axes2['B'].axvline(x=np.log(2))
    axes2['B'].axvline(x=-np.log(2))
    axes2['B'].set_xlabel('$\delta$')
    axes2['B'].set_ylabel("$\\mathrm{median}(\\mathrm{Var}(H^\prime)) - \min(\\mathrm{Var}(H^\prime))$")
    axes2['B'].set_ylabel('$\\mathrm{median}(|\langle H \\rangle - \langle H^\prime \\rangle|) - \min$')
    axes2['B'].set_title(f'Change in energy when $g = 0\\rightarrow 10^{{{int(np.log10(g_pert))}}}$')
    axes2['B'].legend(loc='best')

    # print(f'{xvals=}')
    # print(f'{yCvals=}')
    # lineC = axes1['C'].errorbar(x=xCvals, y=yCvals, yerr=None, color=color, marker=marker)
    # lineC.set_label(label)
    # axes1['C'].axvline(x=np.log(2))
    # axes1['C'].axvline(x=-np.log(2))
    # axes1['C'].set_xlabel('$\delta$')
    # # axes['C'].set_ylabel("$\\mathrm{median}(\\mathrm{Var}(H^\prime)) - \min(\\mathrm{Var}(H^\prime))$")
    # axes1['C'].set_ylabel('$\\mathrm{median}(\\frac{|E - E^\prime|}{E\cdot L})$')
    # axes1['C'].set_title(f'Relative change in energy when $g = 0\\rightarrow 10^{{{int(np.log10(g_pert))}}}$')
    # axes1['C'].legend(loc='best')


deltas= [0.0,]
sns.set_palette("viridis",  len(h5files_new))
fig3, axes3 = plt.subplot_mosaic('''
                                AB
                                CD
                                ''',
                               figsize=(1152 * px, 864 * px),
                               layout="tight",
                               )
for h5file, marker in h5files_new:
    for d in reversed(deltas):
        xvals = []
        dvals = []
        gvals = []
        nvals = []
        Lvals = []
        yA3vals, eA3vals = [], []
        yB3vals, eB3vals = [], []
        yC3vals, eC3vals = [], []
        yD3vals, eD3vals = [], []
        with h5py.File(h5file, 'r') as f:
            for key, val in f.items():
                g_pert = val['g_pert'][0]
                L = val['model_size'][0]
                delta = val['delta'][0]
                if delta != d:
                    continue

                num = len(val)
                xvals.append(g_pert)
                gvals.append(g_pert)
                dvals.append(delta)
                nvals.append(num)
                Lvals.append(L)
                label = f'$L={{{np.min(Lvals)}}}$ $\delta = {{{delta:.2f}}}$ ($n_\min={{{np.min(nvals)}}}$)'

                ene_max = np.ones(num)
                # if 'energy_upper_bound' in val.dtype.names:
                #     ene_max = val['energy_change'][()]
                #     print(ene_max)
                # else:
                #     continue
                print(ene_max)
                var_num = len(val['variance'][()])
                var_old = np.abs(val['variance'][()])
                var_prt = np.abs(val['variance_pert'][()])
                var_abs_chn = np.median(np.abs(var_prt - var_old))
                var_rel_chn = np.median(np.abs(var_prt - var_old) / var_old)
                var_frc_chn = np.median(var_prt / var_old)
                var_abs_chn_ste = np.std(np.abs(var_prt - var_old)) / np.sqrt(var_num)
                var_rel_chn_ste = np.std(np.abs(var_prt - var_old) / var_old) / np.sqrt(var_num)
                var_frc_chn_ste = np.std(var_prt / var_old) / np.sqrt(var_num)

                yA3vals.append(var_abs_chn/L)
                eA3vals.append(var_abs_chn_ste/L)
                yA3titl = f'Absolute change/L in energy variance'
                yA3labl = "${\\mathrm{median}(\sigma_{H^\prime}^2 - \sigma_{H}^2)}/{L}$"

                yB3vals.append(var_rel_chn/L)
                eB3vals.append(var_rel_chn_ste/L)
                yB3titl = f'Relative change/L in energy variance'
                yB3labl = "${\\mathrm{median}(\sigma_{H^\prime}^2 - \sigma_{H}^2)}/{\sigma_{H}^2 L}$"


                yC3vals.append(var_frc_chn)
                eC3vals.append(var_frc_chn_ste)
                yC3titl = f'Ratio change/L in energy variance '
                yC3labl = "${\\mathrm{median}(\sigma_{H^\prime}^2 /\sigma_{H}^2 )}$"


                ene_num = len(val['energy'][()])
                ene_old = val['energy'][()]
                ene_prt = val['energy_pert'][()]
                enech_avg = np.median(np.abs(ene_old - ene_prt))
                ene_abs_chn = np.median(np.abs(ene_prt - ene_old))
                ene_rel_chn = np.median(np.abs(ene_prt - ene_old) / np.abs(ene_old))
                ene_rel_ste = np.std(np.abs(ene_old - ene_prt)) / np.sqrt(ene_num)
                ene_abs_ste = np.std(np.abs(ene_prt - ene_old)) / np.sqrt(ene_num)

                yD3vals.append(ene_abs_chn)
                eD3vals.append(ene_abs_ste)
                yD3titl = f'Absolute change/L in energy$'
                yD3labl = "${\\mathrm{median}(|E' - E|)}$"

                # yB3vals.append(ene_rel_chn)
                # eB3vals.append(ene_rel_ste)


                # if 'hsquared' in val.dtype.names:
                #     hsq_num = len(val['hsquared'][()])
                #     hsq_old = val['hsquared'][()]
                #     hsq_prt = val['hsquared_pert'][()]
                #     hsqch_avg = np.median(np.abs(hsq_old - hsq_prt))
                #     hsq_abs_chn = np.median(np.abs(hsq_prt - hsq_old))
                #     hsq_rel_chn = np.median(np.abs(hsq_prt - hsq_old) / np.abs(hsq_old))
                #     hsq_frc_chn = np.median(hsq_prt / hsq_old)
                #     hsq_rel_ste = np.std(np.abs(hsq_prt - hsq_old) / np.abs(hsq_old)) / np.sqrt(hsq_num)
                #     hsq_abs_ste = np.std(np.abs(hsq_prt - hsq_old)) / np.sqrt(hsq_num)
                #     hsq_frc_ste = np.std(np.abs(hsq_prt / hsq_old)) / np.sqrt(hsq_num)
                #
                #     yCvals.append(hsq_rel_chn/L)
                #     eCvals.append(hsq_rel_ste/L)
                #     xCvals.append(val['delta'][0])
                # print(f'{key}: {y1vals[-1]:.3e} +- {varch_ste:.3e}')

        sort = np.argsort(np.asarray(xvals))
        xvals = np.asarray(xvals)[sort]
        yA3vals = np.asarray(yA3vals)[sort]
        eA3vals = np.asarray(eA3vals)[sort]
        yB3vals = np.asarray(yB3vals)[sort]
        eB3vals = np.asarray(eB3vals)[sort]
        yC3vals = np.asarray(yC3vals)[sort]
        eC3vals = np.asarray(eC3vals)[sort]
        yD3vals = np.asarray(yD3vals)[sort]
        eD3vals = np.asarray(yD3vals)[sort]
        # yA3vals = np.asarray(y32vals)[sort]
        # eA3vals = np.asarray(e32vals)[sort]
        # yB3vals = np.asarray(y32vals)[sort]
        # eB3vals = np.asarray(e32vals)[sort]

        # if(len(yCvals) > 0):
        #     yCvals = np.asarray(yCvals)[sort]
        #     eCvals = np.asarray(eCvals)[sort]
        #     xCvals = np.asarray(xCvals)[sort]


        lineA1 = axes3['A'].errorbar(x=xvals, y=yA3vals, yerr=None,  marker=marker)
        lineA1.set_label(label)
        axes3['A'].set_xlabel('$g$')
        axes3['A'].set_ylabel(yA3labl)
        axes3['A'].set_title(yA3titl)
        axes3['A'].set_yscale('log')
        axes3['A'].set_xscale('log')
        axes3['A'].legend(loc='best')

        lineB3 = axes3['B'].errorbar(x=xvals, y=yB3vals, yerr=None, marker=marker)
        lineB3.set_label(label)
        axes3['B'].set_xlabel('$g$')
        axes3['B'].set_ylabel(yB3labl)
        axes3['B'].set_title(yB3titl)
        axes3['B'].set_yscale('log')
        axes3['B'].set_xscale('log')
        axes3['B'].legend(loc='best')

        lineC3 = axes3['C'].errorbar(x=xvals, y=yC3vals, yerr=None,  marker=marker)
        lineC3.set_label(label)
        axes3['C'].set_xlabel('$g$')
        axes3['C'].set_ylabel(yC3labl)
        axes3['C'].set_title(yC3titl)
        axes3['C'].set_yscale('log')
        axes3['C'].set_xscale('log')
        axes3['C'].legend(loc='best')



        lineD3 = axes3['D'].errorbar(x=xvals, y=yD3vals, yerr=None, marker=marker)
        lineD3.set_label(label)
        axes3['D'].set_xlabel('$g$')
        axes3['D'].set_ylabel(yD3labl)
        axes3['D'].set_title(yD3titl)
        axes3['D'].set_yscale('log')
        axes3['D'].set_xscale('log')
        axes3['D'].legend(loc='best')

        # lineB2 = axes2['B'].errorbar(x=xvals, y=yB2vals, yerr=None, color=color, marker=marker)
        # lineB2.set_label(label)
        # axes3['D'].axvline(x=np.log(2))
        # axes3['D'].axvline(x=-np.log(2))
        # axes3['D'].set_xlabel('$\delta$')
        # axes3['D'].set_ylabel("$\\mathrm{median}(\\mathrm{Var}(H^\prime)) - \min(\\mathrm{Var}(H^\prime))$")
        # axes3['D'].set_ylabel('$\\mathrm{median}(|\langle H \\rangle - \langle H^\prime \\rangle|) - \min$')
        # axes3['D'].set_title(f'Change in energy when $g = 0\\rightarrow 10^{{{int(np.log10(g_pert))}}}$')
        # axes3['D'].legend(loc='best')



fig1.suptitle('Energy Variance Sensitivity')
fig2.suptitle('Energy Sensitivity')
fig3.suptitle('Sensitivity with increasing $g$')

plt.tight_layout()
# plt.legend(loc='best')
plt.show()
