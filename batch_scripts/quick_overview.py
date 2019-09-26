import h5py
import os
import warnings
import re
import numpy as np
from datetime import datetime
import argparse

parser = argparse.ArgumentParser(description='Quick overview of batch simulation')
parser.add_argument('-S','--summary', action='store_true', help='Summary only')
parser.add_argument('-s','--save', action='store_true', help='Save to file')
parser.add_argument('-f','--filename', type=str, help='Save to file with filename',default='overview')
parser.add_argument('-t','--timestamp', action='store_true', help='Add timestamp to filename')
parser.add_argument('-d','--directory', type=str, help='Search for hdf5 files in directory', default='output')
args = parser.parse_args()

if args.timestamp:
    args.filename = args.filename + '-' + datetime.today().strftime('%Y-%m-%d-%H:%M:%S')

if args.filename != 'overview':
    args.save = True

regex = re.compile(r'\d+')

if args.save:
    file = open(args.filename, 'w')

for dirName, subdirList, fileList in os.walk(args.directory):
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

    if not args.summary:
        header = "{:15} {:12} {:>12} {:>12} {:>12} {:>12} {:>12}".format("Realization", "Variance","Time", "Resets", "Converged", "Saturated", "Succeeded")
        print(header)
        if args.save:
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

            if not args.summary:
                entry = "{:<15} {:<12.4f} {:>12.4f} {:>12} {:>12} {:>12} {:>12}".format(
                    realization_num[-1],
                    np.log10(variance[-1]),
                    walltime[-1]/60,
                    resets[-1],
                    converged[-1],
                    saturated[-1],
                    succeeded[-1])
                print(entry)
                if args.save:
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
    if args.save:
        file.write(header)
        file.write(entry)
if args.save:
    file.close()
