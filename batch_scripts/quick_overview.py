import h5py
import os
import warnings
import re
import numpy as np
from datetime import datetime
import argparse
import subprocess
import colored
from colored import stylize
parser = argparse.ArgumentParser(description='Quick overview of batch simulation')
parser.add_argument('-S', '--summary', action='store_true', help='Summary only')
parser.add_argument('-s', '--save', action='store_true', help='Save to file')
parser.add_argument('-f', '--filename', type=str, help='Save to file with filename', default='experiment')
parser.add_argument('-t', '--timestamp', action='store_true', help='Add timestamp to filename')
parser.add_argument('-d', '--directory', type=str, help='Search for hdf5 files in directory', default='output')
parser.add_argument('-o', '--outdir', type=str, help='Save output to directory', default='experiments')
parser.add_argument('-x', '--suffix', type=str, help='Append suffix', default='.log')
args = parser.parse_args()


if args.timestamp:
    args.filename = args.filename + '-' + datetime.today().strftime('%Y-%m-%d-%H:%M:%S')

if args.filename != 'experiment' or args.outdir != 'experiments'  :
    args.save = True

try:
    git_rev = "GIT INFO:\n" + str(subprocess.check_output(['git', 'log', '-1']).decode('ascii'))
    print(git_rev)
except:
    git_rev=''

regex = re.compile(r'\d+')

if args.save:
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)
    file = open(args.outdir + '/'+ args.filename + args.suffix, 'w')
    file.write(git_rev + '\n')

for dirName, subdirList, fileList in os.walk(args.directory):
    if subdirList:
        subdirList.sort()
    if not fileList:
        continue
    fileList.sort()
    chainlen  = []
    realization_num = []
    variance  = []
    walltime  = []
    resets    = []
    got_stuck = []
    saturated = []
    converged = []
    succeeded = []
    finished  = []



    if not args.summary:
        header = "{:8} {:15} {:>12} {:>12} {:>12} {:>12} {:>12} {:>12} {:>12} {:>12}".format("Length", "Realization", "Variance", "Time",
                                                                                             "Resets", "Got Stuck", "Saturated" ,"Converged",
                                                                                             "Succeeded", "Finished")
        print(header)
        if args.save:
            file.write(header + '\n')
    for h5path in fileList:
        # print("Filepath: ", h5path)
        try:
            h5file = h5py.File(dirName + '/' + h5path, 'r', swmr=True)
            if not isinstance(h5file, h5py.File):
                raise IOError("File not readable")
        except Exception as er:
            print("Could not open file [", h5path, "] Reason: ", er)
            continue

        try:
            realization_name = h5path.replace('.h5', '')
            chainlen.append(h5file['xDMRG/measurements/length'][-1])
            realization_num.append([int(x) for x in regex.findall(realization_name)][-1])
            variance.append(h5file['xDMRG/measurements/energy_variance_per_site'][-1])
            walltime.append(h5file['xDMRG/sim_status/wall_time'][-1])
            resets.append(h5file['xDMRG/sim_status/num_resets'][-1])
            saturated.append(h5file['xDMRG/sim_status/simulation_has_saturated'][-1])
            converged.append(h5file['xDMRG/sim_status/simulation_has_converged'][-1])
            succeeded.append(h5file['xDMRG/sim_status/simulation_has_succeeded'][-1])
            try:
                finished.append(h5file['common/finOK'][-1])
            except:
                finished.append(0)

            try:
                got_stuck.append(h5file['xDMRG/sim_status/simulation_has_got_stuck'][-1])
            except:
                if finished[-1] == 0 and converged[-1] == 0 and saturated[-1] == 1:
                    got_stuck.append(1)
                else:
                    got_stuck.append(0)

            h5file.close()

            style = ''
            if finished[-1] == 1 and succeeded[-1] == 1:
                style = colored.bg("green_4")
            elif finished[-1] == 1 and converged[-1] == 1:
                style = colored.bg("hot_pink_2")
            elif finished[-1] == 1 and got_stuck[-1] == 1:
                if variance[-1] < 1e-10:
                    style = colored.bg("dark_green_sea")
                else:
                    style = colored.bg("red_3b")
            elif finished[-1] == 0 and got_stuck[-1] == 1:
                style = colored.fg("red_3b")
            elif finished[-1] == 0 and got_stuck[-1] == 0 and converged[-1] == 0:
                style = ''
            elif finished[-1] == 0 and got_stuck[-1] == 0 and converged[-1] == 1:
                style = colored.fg("green_4")



            if not args.summary:
                entry = "{:<8} {:<15} {:>12.4f} {:>12.4f} {:>12} {:>12} {:>12} {:>12} {:>12} {:>12}".format(
                    chainlen[-1],
                    realization_num[-1],
                    np.log10(variance[-1]),
                    walltime[-1] / 60,
                    resets[-1],
                    got_stuck[-1],
                    saturated[-1],
                    converged[-1],
                    succeeded[-1],
                    finished[-1])
                print(stylize(entry, style))

                if args.save:
                    file.write(entry + '\n')

        except Exception as er:
            print("Could not read dataset. Reason: ", er)
            continue
    header = "{:8} {:6} {:>12} {:<12} {:>8} {:>12} {:>12} {:>12} {:>12} {:>12} {:>12}".format("Length","Sims", "log10 <Var>", "<log10 Var>", "Avg Time",
                                                                                              "Avg Resets","Sum Stk", "Sum Sat", "Sum Con",
                                                                                              "Sum Suc", "Sum Fin")
    entry = "{:<8} {:<6} {:>12.4f} {:<12.4f} {:>8.4f} {:>12.3f} {:>12} {:>12} {:>12} {:>12} {:>12}".format(
        np.int(np.mean(chainlen)),
        len(realization_num),
        np.log10(np.nanmean(variance)),
        np.nanmean(np.log10(variance)),
        np.mean(walltime) / 60,
        np.mean(resets),
        np.sum(got_stuck),
        np.sum(saturated),
        np.sum(converged),
        np.sum(succeeded),
        np.sum(finished))

    print(header)
    print(entry)
    if args.save:
        file.write(header + '\n')
        file.write(entry + '\n')
if args.save:
    file.close()

print("Legend:")
print("Success implies converging to variance < 1e-12")
print(stylize(colored.bg("green_4")        , "Finished: success"))
print(stylize(colored.bg("red_3b")         , "Finished: no success"))
print(stylize(colored.bg("dark_green_sea") , "Finished: almost success (variance < 1e-10)"))
print(stylize(colored.bg("hot_pink_2")     , "Finished: meeting criteria fur success but not successfully (probably logic error)"))
print(stylize(colored.fg("green_4")        , "Running, reached success"))
print(stylize(colored.fg("red_3b")         , "Running, currently stuck"))
print("Running")

