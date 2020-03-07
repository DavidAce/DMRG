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

for dirName, subdirList, fileList in os.walk(args.directory):
    if subdirList:
        subdirList.sort()
    if not fileList:
        continue

    if not args.length in dirName:
        continue

    fileList.sort()
    chainlen  = []
    seed      = []
    iter      = []
    step      = []
    energy    = []
    variance  = []
    variancel = []
    ententrp  = []
    entdiff   = []
    walltime  = []
    resets    = []
    got_stuck = []
    saturated = []
    converged = []
    succeeded = []
    finished  = []



    if not args.summary:
        header = "{:<8} {:<6} {:<6} {:<6} {:>12} {:>12} {:>12} {:>12} {:>12} {:>12} {:>8} {:>5} {:>5} {:>5} {:>5} {:>5}".format("Length", "Seed", "Iter","Step", "Energy", "VarNow", "VarLow","Ent.Entr.", "Ent.Diff", "Time",
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
        state_keys = ""
        try:
            if 'xDMRG' in h5file:
                state_keys = [str(x) for x in h5file['xDMRG'].keys() if 'state' in x]
        except Exception as err:
            print("Could not read keys in h5file. Reason:", err)
        entry = []
        ententrp_zero = []
        for state_num,state_key in enumerate(state_keys):
            if len(args.state) > 0 and not args.state in state_key:
                continue

            table_path = ''
            step_results = 0
            step_journal = 0
            prefix = 'xDMRG/' + state_key
            try:
                if prefix + '/journal/sim_status' in h5file:
                    step_array =  h5file[prefix + '/journal'].get('sim_status')['step']
                    if(len(step_array)>0):
                        step_journal = step_array[-1]
                if prefix + '/results/sim_status' in h5file:
                    step_array = h5file[prefix + '/results'].get('sim_status')['step']
                    if(len(step_array)>0):
                        step_results = step_array[-1]
            except Exception as err:
                print("Could not read sim_status. Reason:", err)
            if step_results == 0 and step_journal == 0:
                continue
            if step_journal >= step_results:
                table_path = prefix + '/journal'
            else:
                table_path = prefix + '/results'

            # print(table_path,step_journal, step_results)
            try:
                try:
                    finished.append(h5file['common/finOK'][-1])
                except:
                    finished.append(0)
                if(args.finished and finished[-1] == 0):
                    continue


                realization_name = h5path.replace('.h5', '')
                chainlen       .append(h5file[table_path].get('measurements')['length'][-1])
                seed           .append([int(x) for x in regex.findall(realization_name)][-1])
                iter           .append(h5file[table_path].get('sim_status')['iteration'][-1])
                step           .append(h5file[table_path].get('sim_status')['step'][-1])
                energy         .append(h5file[table_path].get('measurements')['energy_per_site'][-1])
                variance       .append(h5file[table_path].get('measurements')['energy_variance_per_site'][-1])
                #variancel       .append(h5file[table_path].get('measurements')['energy_variance_per_site_lowest'][-1])
                variancel       .append(get_data(h5file[table_path].get('measurements'),'energy_variance_per_site_lowest'))
                ententrp       .append(h5file[table_path].get('measurements')['entanglement_entropy_midchain'][-1])
                walltime       .append(h5file[table_path].get('sim_status')['wall_time'][-1])
                resets         .append(h5file[table_path].get('sim_status')['num_resets'][-1])
                got_stuck      .append(h5file[table_path].get('sim_status')['simulation_has_got_stuck'][-1])
                saturated      .append(h5file[table_path].get('sim_status')['simulation_has_saturated'][-1])
                converged      .append(h5file[table_path].get('sim_status')['simulation_has_converged'][-1])
                succeeded      .append(h5file[table_path].get('sim_status')['simulation_has_succeeded'][-1])

                ententrp_curr = h5file[table_path]['entanglement_entropies'][()]
                if len(entdiff) == 0 or len(ententrp_zero) == 0:
                    ententrp_zero = h5file[table_path]['entanglement_entropies'][()]
                    entdiff.append(np.nan)
                else:
                    entdiff.append(np.sum((ententrp_curr-ententrp_zero)/np.log(2)))



                style = ''
                if finished[-1] == 1 and succeeded[-1] == 1:
                    style = colored.bg("green_4")
                elif finished[-1] == 1:
                    if variance[-1] < 1e-10:
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
                    entry.append("{:<8} {:<6} {:<6} {:<6} {:>12.4f} {:>12.4f} {:>12.4f} {:>12.4f} {:>12.4f} {:>12.4f} {:>8} {:>5} {:>5} {:>5} {:>5} {:>5}".format(
                        chainlen[-1],
                        seed[-1],
                        iter[-1],
                        step[-1],
                        energy[-1],
                        np.log10(variance[-1]),
                        np.log10(variancel[-1]),
                        ententrp[-1],
                        entdiff[-1],
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
    if np.isnan(entdiff).all():
        entdiff = np.zeros(len(entdiff))

    header = "{:<8} {:<6} {:<6} {:<6} {:>12} {:>12} {:>12} {:>12} {:>12} {:>8} {:>5} {:>5} {:>5} {:>5} {:>5}".format("Length","Sims", "<iter>","<step>","<Energy>", "<VarNow>","<VarLow>","<Entgl>","<Entgl-diff>","<Time>",
                                                                                                              "Resets","Stk", "Sat", "Con",
                                                                                                              "Suc", "Fin")
    entry = "{:<8} {:<6} {:<6.1f} {:<6.1f} {:>12.4f} {:>12.4f} {:>12.4f} {:>12.4f} {:>12.4f} {:>12.3f} {:>8.1f} {:>5} {:>5} {:>5} {:>5} {:>5}".format(
        np.nanmax(chainlen),
        len(seed),
        np.nanmean(iter),
        np.nanmean(step),
        np.nanmean(energy),
        np.nanmean(np.log10(variance)),
        np.nanmean(np.log10(variancel)),
        np.nanmean(ententrp),
        np.nanmean(entdiff),
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

