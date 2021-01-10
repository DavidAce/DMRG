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
parser.add_argument('-a', '--algorithms',action='append', type=str, help='Consider these algorithms', default=[])
parser.add_argument('-S', '--summary', action='store_true', help='Summary only')
parser.add_argument('-p', '--projection', action='store_true', help='Include projection')
parser.add_argument('-s', '--save', action='store_true', help='Save to file')
parser.add_argument('-f', '--filename', type=str, help='Save to file with filename', default='experiment')
parser.add_argument('-F', '--finished', action='store_true', help='Only consider finished simulations')
parser.add_argument('-r', '--state', type=str, help='Consider only a certain state number', default='')
parser.add_argument('-L', '--length', type=str, help='Consider only a certain chain length', default='')
parser.add_argument('-t', '--timestamp', action='store_true', help='Add timestamp to filename')
parser.add_argument('-d', '--directory', type=str, help='Search for hdf5 files in directory', default='output')
parser.add_argument('-o', '--outdir', type=str, help='Save output to directory', default='experiments')
parser.add_argument('-x', '--suffix', type=str, help='Append suffix', default='.log')
args = parser.parse_args()

def get_data(h5obj,name,idx=-1):
    try:
        return h5obj[name][idx]
    except Exception as er:
        return 0

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

if not args.algorithms:
    args.algorithms = ['xDMRG']

if args.algorithms:
    print("Checking algorithms: ", args.algorithms)

for dirName, subdirList, fileList in os.walk(args.directory):
    if subdirList:
        subdirList.sort()
    if not fileList:
        continue

    if not args.length in dirName:
        continue

    fileList.sort()
    algorithm = []
    state     = []
    chainlen  = []
    seed      = []
    iter      = []
    step      = []
    energy    = []
    variance  = []
    variancel = []
    ententrp  = []
    walltime  = []
    resets    = []
    got_stuck = []
    saturated = []
    converged = []
    succeeded = []
    finished  = []



    if not args.summary:
        header = "{:<10} {:<12} {:<8} {:<6} {:<6} {:<6} {:>12} {:>12} {:>12} {:>12} {:>12} {:>8} {:>5} {:>5} {:>5} {:>5} {:>5}".format("Algorithm", "State" ,"Length", "Seed", "Iter","Step", "Energy", "VarNow", "VarLow","Ent.Entr.", "Time",
                                                                                                           "Resets", "Stk", "Sat" ,"Con", "Suc", "Fin")

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

        # Collect the number of states present.
        try:
            if not "common/state_root" in h5file:
                continue
        except Exception as er:
            print("Could not read [common/state_root]. Reason: ", er)
            continue
        try:
            state_keys = [x for x in h5file["common/state_root"].attrs.keys() if any(algo in x for algo in args.algorithms)]
            state_keys = []
            for candidate in h5file["common/finished"].attrs.keys():
                if args.finished and h5file["common/finished"].attrs[candidate] == 0:
                    continue
                if not args.finished and "savepoint" in candidate and h5file["common/state_root"].attrs[candidate] + "/finished" in h5file:
                    continue
                if not args.finished and "checkpoint" in candidate and h5file["common/state_root"].attrs[candidate] + "/finished" in h5file:
                    continue
                if not args.projection and "projection" in candidate:
                    continue
                if not any(algo in candidate for algo in args.algorithms):
                    continue
                if args.state and not str(args.state) in h5file["common/state_root"].attrs[candidate]:
                    continue
                state_keys.append(candidate)
            state_keys.sort()
        except Exception as er:
            print("Could not gather paths in file [",h5path,"]. Reason: ", er)
            continue
        for state_num,state_prefix in enumerate(state_keys):
            entry = []
            ententrp_zero = []
            try:
                algo_name = h5file["common/algo_type"].attrs[state_prefix]
                state_name = h5file["common/state_name"].attrs[state_prefix]
                finished.append(h5file["common/finished"].attrs[state_prefix])
                if (args.finished and finished[-1] == 0):
                    continue
                msrmnt_last_entry = h5file[state_prefix].get('measurements')[-1]
                status_last_entry = h5file[state_prefix].get('status')[-1]
                realization_name = h5path.replace('.h5', '')
                algorithm.append(algo_name)
                state.append(state_name)
                chainlen.append(msrmnt_last_entry['length'])
                seed.append([int(x) for x in regex.findall(realization_name)][-1])
                iter.append(status_last_entry['iter'])
                step.append(status_last_entry['step'])
                energy.append(msrmnt_last_entry['energy_per_site'])
                variance.append(msrmnt_last_entry['energy_variance_per_site'])
                variancel.append(status_last_entry['energy_variance_lowest'])
                ententrp.append(msrmnt_last_entry['entanglement_entropy_midchain'])
                walltime.append(status_last_entry['wall_time'])
                resets.append(status_last_entry['num_resets'])
                got_stuck.append(status_last_entry['algorithm_has_got_stuck'])
                saturated.append(status_last_entry['algorithm_has_saturated'])
                converged.append(status_last_entry['algorithm_has_converged'])
                succeeded.append(status_last_entry['algorithm_has_succeeded'])

                style = ''
                if finished[-1] == 1:
                    if variance[-1] < 1e-12 or succeeded[-1] == 1:
                        style = colored.bg("green_4")
                    elif variance[-1] < 1e-10:
                        style = colored.bg("dark_green_sea")
                    elif variance[-1] < 1e-8:
                        style = colored.bg("dark_orange")
                    else:
                        style = colored.bg("red_3b")
                elif finished[-1] == 0 and got_stuck[-1] == 1:
                    style = colored.fg("red_3b")
                    if variance[-1] < 1e-10:
                        style = colored.fg("dark_green_sea")
                    else:
                        style = colored.fg("dark_orange")
                elif finished[-1] == 0 and got_stuck[-1] == 0 and converged[-1] == 0:
                    style = ''
                elif finished[-1] == 0 and got_stuck[-1] == 0 and converged[-1] == 1:
                    style = colored.fg("green_4")

                if not args.summary:
                    entry.append(
                        "{:<10} {:<12} {:<8} {:<6} {:<6} {:<6} {:>12.4f} {:>12.4f} {:>12.4f} {:>12.4f} {:>12.4f} {:>8} {:>5} {:>5} {:>5} {:>5} {:>5}".format(
                            algorithm[-1],
                            state[-1],
                            chainlen[-1],
                            seed[-1],
                            iter[-1],
                            step[-1],
                            energy[-1],
                            np.log10(variance[-1]),
                            np.log10(variancel[-1]),
                            ententrp[-1],
                            walltime[-1] / 60,
                            resets[-1],
                            got_stuck[-1],
                            saturated[-1],
                            converged[-1],
                            succeeded[-1],
                            finished[-1]))
                    print(stylize(entry[-1], style))

                    if args.save:
                        file.write(entry[-1] + '\n')

            except Exception as er:
                print("Could not read dataset. Reason: ", er)
                continue
        h5file.close()
    if not fileList:
        continue
    if len(chainlen) == 0:
        continue

    header = "{:<23} {:<8} {:<6} {:<6} {:<6} {:>12} {:>12} {:>12} {:>12} {:>12} {:>8} {:>5} {:>5} {:>5} {:>5} {:>5}".format("","Length","Sims", "<iter>","<step>","<Energy>", "<VarNow>","<VarLow>","<Entgl>","<Time>",
                                                                                                              "Resets","Stk", "Sat", "Con",
                                                                                                              "Suc", "Fin")
    entry = "{:<23} {:<8} {:<6} {:<6.1f} {:<6.1f} {:>12.4f} {:>12.4f} {:>12.4f} {:>12.4f} {:>12.3f} {:>8.1f} {:>5} {:>5} {:>5} {:>5} {:>5}".format( "",
        np.nanmax(chainlen),
        len(seed),
        np.nanmean(iter),
        np.nanmean(step),
        np.nanmean(energy),
        np.nanmean(np.log10(variance)),
        np.nanmean(np.log10(variancel)),
        np.nanmean(ententrp),
        np.nanmean(walltime) / 60,
        np.sum(resets),
        np.sum(got_stuck),
        np.sum(saturated),
        np.sum(converged),
        np.sum(succeeded),
        np.sum(finished))
    print("="*len(header))
    print(header)
    print(entry)
    if args.save:
        file.write(header + '\n')
        file.write(entry + '\n')
if args.save:
    file.close()

print("Legend:")
print(stylize("Finished : success        (variance < 1e-12)"                                        , colored.bg("green_4")))
print(stylize("Finished : almost success (variance < 1e-10)"                                        , colored.bg("dark_green_sea")))
print(stylize("Finished : mediocre run   (variance < 1e-8)"                                         , colored.bg("dark_orange")))
print(stylize("Finished : failed         (variance > 1e-8)"                                         , colored.bg("red_3b")))
print(stylize("Finished : meeting criteria for success but not successfully (probably logic error)" , colored.bg("hot_pink_2")))
print(stylize("Running  : reached success"                                                          , colored.fg("green_4")))
print(stylize("Running  : almost success  (variance < 1e-10)"                                       , colored.fg("dark_green_sea")))
print(stylize("Running  : currently stuck"                                                          , colored.fg("dark_orange")))
print(        "Running  : just running")

