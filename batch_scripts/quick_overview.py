import h5py
import os
import warnings
import re
import numpy as np
from datetime import datetime
h5directory = 'output'
regex = re.compile(r'\d+')

timestamp = datetime.today().strftime('%Y-%m-%d-%H:%M:%S')
file = open('overview-' + str(timestamp) + '.log', 'w')

for dirName, subdirList, fileList in os.walk(h5directory):
    if not fileList:
        continue
    subdirList.sort()
    fileList.sort()
    realization_num  = []
    variance         = []
    walltime         = []
    resets           = []
    converged        = []
    saturated        = []
    succeeded        = []

    header = "{:15} {:12} {:>12} {:>12} {:>12} {:>12} {:>12}".format("Realization", "Variance","Time", "Resets", "Converged", "Saturated", "Succeeded")
    print(header)
    file.write(header)
    for h5path in fileList:
        # print("Filepath: ", h5path)
        try:
            h5file = h5py.File(dirName + '/' + h5path, 'r',swmr=True)
            if not isinstance(h5file, h5py.File):
                raise IOError("File not readable")
        except Exception as er:
            print("Could not open file [", h5path,"] Reason: ",er)
            continue

        try:
            realization_name = h5path.replace('.h5', '')
            realization_num .append([int(x) for x in regex.findall(realization_name)][-1])
            variance .append(h5file['xDMRG/measurements/energy_variance_per_site'][-1])
            walltime .append(h5file['xDMRG/sim_status/wall_time'][-1])
            resets   .append(h5file['xDMRG/sim_status/num_resets'][-1])
            converged.append(h5file['xDMRG/sim_status/simulation_has_converged'][-1])
            saturated.append(h5file['xDMRG/sim_status/simulation_has_saturated'][-1])
            succeeded.append(h5file['xDMRG/sim_status/simulation_has_succeeded'][-1])
            entry = "{:<15} {:<12.4f} {:>12.4f} {:>12} {:>12} {:>12} {:>12}".format(
                realization_num[-1],
                np.log10(variance[-1]),
                walltime[-1]/60,
                resets[-1],
                converged[-1],
                saturated[-1],
                succeeded[-1])

            print(entry)
            file.write(entry)

        except Exception as er:
            print("Could not read dataset. Reason: ", er)
            continue
    header="{:15} {:12} {:>12} {:>12} {:>12} {:>12} {:>12}".format("Total sims", "Avg Var","Avg Time", "Avg Resets", "Sum Con", "Sum Sat", "Sum Suc")
    entry="{:<15} {:<12.4f} {:>12.4f} {:>12.3f} {:>12} {:>12} {:>12}".format(
        len(realization_num),
        np.log10(np.nanmean(variance)),
        np.mean(walltime) / 60,
        np.mean(resets),
        np.sum(converged),
        np.sum(saturated),
        np.sum(succeeded))

    print(header)
    print(entry)
    file.write(header)
    file.write(entry)
file.close()
