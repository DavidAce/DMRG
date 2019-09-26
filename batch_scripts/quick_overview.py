import h5py
import os
import warnings
import re
import numpy as np
h5directory = 'output'
regex = re.compile(r'\d+')


for dirName, subdirList, fileList in os.walk(h5directory):
    subdirList.sort()
    fileList.sort()
    print("{:15} {:8} {:>8} {:>8} {:>10} {:>10} {:>10}".format("Realization", "Variance","Time", "Resets", "Converged", "Saturated", "Succeeded"))
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
            realization_num = [int(x) for x in regex.findall(realization_name)][-1]
            variance = h5file['xDMRG/measurements/energy_variance_per_site'][-1]
            walltime = h5file['xDMRG/sim_status/wall_time'][-1]
            resets   = h5file['xDMRG/sim_status/num_resets'][-1]
            converged  = h5file['xDMRG/sim_status/simulation_has_converged'][-1]
            saturated  = h5file['xDMRG/sim_status/simulation_has_saturated'][-1]
            succeeded  = h5file['xDMRG/sim_status/simulation_has_succeeded'][-1]
            print("{:<15} {:<8.4f} {:>8.4f} {:>8} {:>10} {:>10} {:>10}".format(
                realization_num,
                np.log10(variance),
                walltime/60,
                resets,
                converged,
                saturated,
                succeeded))

        except Exception as er:
            print("Could not read dataset. Reason: ", er)
            continue

