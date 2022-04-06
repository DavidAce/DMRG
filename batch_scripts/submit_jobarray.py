import argparse
import os
import errno
import subprocess
from pathlib import Path
from itertools import zip_longest
from itertools import chain, islice
from datetime import datetime



def chunks(iterable, n):
   "chunks(ABCDE,2) => AB CD E"
   iterable = iter(iterable)
   while True:
       # store one line in memory,
       # chain it to an iterator on the rest of the chunk
       try:
           yield chain([next(iterable)], islice(iterable, n-1))
       except StopIteration:
           return




def parse(project_name):
    parser = argparse.ArgumentParser(description='SLURM batch submission for {}'.format(project_name))
    parser.add_argument('-b', '--build-type', type=str, help='Build type', default='Release', choices=['Debug', 'Release', 'RelWithDebInfo', 'Profile'])
    parser.add_argument('-c', '--cluster', type=str, help='Comma separated list of Slurm clusters', default=None, choices=['kraken', 'draken', 'kthulu', 'tetralith'])
    parser.add_argument('--config', type=str, help='File or path to files containing simulation config files (suffixed .cfg)', default='input')
    parser.add_argument('--cpus-per-task', type=int, help='Number of cores per task (sims)', default=1)
    parser.add_argument('--ntasks-per-core', type=int, help='Number of tasks (sims) on each core', default=1)
    parser.add_argument('--dryrun', action='store_true', help='Dry run')
    parser.add_argument('--debug', action='store_true', help='Debug this script')
    parser.add_argument('--exclusive', action='store_true', help='Reserve whole node')
    parser.add_argument('--execname', type=str, help='Name of executable', default='DMRG++')
    parser.add_argument('--hint', type=str, default=None, help='Slurm parallelization hint. Cannot be used with --ntasks-per-core', choices=['multithread', 'nomultithread', 'compute_bound', 'memory_bound'])
    parser.add_argument('-j', '--job-file', type=str, help='File containing config-seed pairs. Use for resuming failed runs')
    parser.add_argument('-J', '--job-name', type=str, help='Slurm job name', default='DMRG')
    parser.add_argument('-m', '--mem-per-cpu', type=int, help='Memory per core, e.g 2000, 2000M or 2G', default=1000)
    parser.add_argument('-N', '--sims-per-cfg', type=int, help='Number of simulations per config file. Can be split up into chunks with -n', default=10)
    parser.add_argument('-n', '--sims-per-array', type=int, help='Number of simulations in each job-array (splits -N into -N/-n chunks)', default=1000)
    parser.add_argument('--sims-per-task', type=int, help='Number of simulations per job-array task. This is equivalent to the step, or stride in the array', default=10)
    parser.add_argument('--ntasks', type=int, help='Number of tasks per simulation (e.g. mpi threads)', default=10)
    parser.add_argument('-o', '--other', type=str, help='Other options for sbatch (verbatim)', default=None)
    parser.add_argument('--open-mode', type=str, help='Access mode for logs', default='append', choices=['append','truncate'])
    parser.add_argument('--partition', type=str, help='Partition name', default=None)
    parser.add_argument('--qos', type=str, help='Quality of service', default=None)
    parser.add_argument('--requeue', action='store_true', help='Requeue job in case of failure')
    parser.add_argument('--seedfile', type=str, help='File or path to files containing a list of seeds (suffixed .seed). Basenames should match the corresponding input files in -f', default=None)
    parser.add_argument('--startseed', type=int, help='Starting seed for random number generator', default=0)
    parser.add_argument('--shuffle', action='store_true', help='Shuffle all seeeds and config files')
    parser.add_argument('-t', '--time', type=str, help='Time limit for each job', default='0-01:00:00' )
    parser.add_argument('-v', '--verbose', action='store_true', help='Enable verbose sbatch')
    parser.add_argument('--default-kraken', action='store_true', help='Set defaults for kraken cluster')
    parser.add_argument('--default-tetralith', action='store_true', help='Set defaults for tetralith cluster')


    args = parser.parse_args()
    if args.default_kraken:
        parser.set_defaults(partition='dedicated',qos='lowprio')
        args = parser.parse_args()
    if args.default_tetralith:
        parser.set_defaults(partition='tetralith')
        args = parser.parse_args()

    return args



def run(cmd,env,args, shell=False, stdout=subprocess.PIPE, stderr=subprocess.STDOUT):
    if len(cmd) == 0:
        return
    if isinstance(cmd,str):
        cmd = cmd.split()
    if args.debug or args.dryrun:
        print("Run command:", ' '.join(cmd))

    if args.dryrun:
        return

    with subprocess.Popen(cmd,bufsize=1, shell=shell, stdout=stdout, stderr=stderr, encoding='utf-8',env=env) as p:
        if p.stdout:
            for line in p.stdout:
                print(line, end='')  # process line here
        if p.returncode:
            raise subprocess.CalledProcessError(p.returncode, ' '.join(p.args))


def generate_sbatch_commands(project_name, args):
    sbatch_cmd = []
    sbatch_arg = []
    sbatch_env = os.environ.copy()  # Add environment variables here

    # Find executable
    exec = '../build/{}'.format(args.execname)
    if(os.access(exec, os.X_OK)):
        print('Found executable:', exec)
    else:
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), exec)


    # Make sure the executable would find its dynamically loaded libraries
    ldd_out = subprocess.check_output('ldd {}'.format(exec), shell=True,text=True)
    print(ldd_out)
    if 'not found' in ldd_out:
      raise FileNotFoundError(errno.ENOENT, "Some dynamic libraries were not found. Perhaps a module needs to be loaded.", exec)





#     sbatch $jobname $cluster $partition $qos $mempercpu $requeue $exclusive $time $other $hint $openmode $verbosity $ntasks $cpuspertask $ntaskspercore \
#       --array=1-$numseeds:$simspertask \
#       run_jobarray.sh -e $exec -f $jobfile

    if args.cluster:
        sbatch_arg.extend(['--cluster={}'.format(args.cluster)])
    if args.cpus_per_task:
        sbatch_arg.extend(['--cpus-per-task={}'.format(args.cpus_per_task)])
    if args.ntasks_per_core:
        sbatch_arg.extend(['--ntasks-per-core={}'.format(args.ntasks_per_core)])
    if args.job_name:
        sbatch_arg.extend(['--job-name={}'.format(args.job_name)])
    if args.mem_per_cpu:
        sbatch_arg.extend(['--mem-per-cpu={}'.format(args.mem_per_cpu)])
    if args.ntasks:
        sbatch_arg.extend(['--ntasks={}'.format(args.ntasks)])
    if args.open_mode:
        sbatch_arg.extend(['--open-mode={}'.format(args.open_mode)])
    if args.partition:
        sbatch_arg.extend(['--partition={}'.format(args.partition)])
    if args.qos:
        sbatch_arg.extend(['--qos={}'.format(args.qos)])
    if args.requeue:
        sbatch_arg.extend(['--requeue'])
    if args.hint:
        sbatch_arg.extend(['--hint={}'.format(args.hint)])
    if args.time:
        sbatch_arg.extend(['--time={}'.format(args.time)])
    if args.verbose:
        sbatch_arg.extend(['-v'])


    if args.job_file:
        # If we are using a job-file we can use it directly and then skip generating new ones
        numseeds = sum(1 for line in open(args.job_file))
        jobstem  = Path(args.job_file).stem
        jobext   = Path(args.job_file).suffix
        jobdir   = Path(args.job_file).parent
        if numseeds > args.sims_per_array:
            # Let's split up the job-file into chunks
            with open(args.job_file) as jobfile:
                for i, lines in enumerate(chunks(jobfile, args.sims_per_array)):
                    jobchnk = '{}/{}.{}.job'.format(jobdir, jobstem,i)
                    with open(jobchnk, 'w') as f:
                        f.writelines(lines)
            os.rename(args.job_file, '{}/{}.bak'.format(jobdir,jobstem)) # Back up the job file
    else:
        cfgfiles = list(Path(args.config).glob('*.cfg')) #TODO: may have to sort this list
        if not cfgfiles:
            raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT) + ' (no .cfg files found in)', args.config)



        # Generate a new unique directory for seed files
        jobdir = "jobs-{}".format(datetime.now().strftime("%Y-%m-%dT%H.%M.%S"))
        os.makedirs(jobdir, exist_ok=True)

        # Generate a super job list of .cfg and seed pairs.
        # When --shuffle is given, we assign seeds in round-robin instead of actually
        # shuffling, which would take too long for lists of this size.
        superlist = []
        seedcount = args.startseed
        if args.shuffle:
            for sim in range(args.sims_per_cfg):
                for cfg in cfgfiles:
                    superlist.append(['{} {}\n'.format(cfg,seedcount)])
                    seedcount = seedcount + 1
        else:
            for cfg in cfgfiles:
                for sim in range(args.sims_per_cfg):
                    superlist.append(['{} {}\n'.format(cfg,seedcount)])
                    seedcount = seedcount + 1


        # Generate job files
        print("Generating {} seeds per .cfg file ({} files in total) split into job arrays of size {}".format(args.sims_per_cfg, len(cfgfiles), args.sims_per_array))

        # Write the super job list to job files in chunks of size sims_per_array
        for i,joblist in enumerate(chunks(superlist, args.sims_per_array)):
            bgn = i * args.sims_per_array
            end = i * args.sims_per_array + args.sims_per_array - 1
            end = min(end, len(superlist)) - 1
            jobfile = '{}/part.{:0>4}.[{}-{}].job'.format(jobdir, i, bgn,end)
            with open(jobfile,'w+') as f:
                for job in joblist:
                    f.writelines(job)


    # From this point on we are guaranteed to have a set of job files jobs/part.[###-###].job
    # Each .job file contains 2 columns with the path to a config file and a seed, and corresponds to one job array.
    # Now collect all .job files and create separate sbatch commands for each
    jobfiles = sorted(list(Path(jobdir).glob('*.job')))
    if not jobfiles:
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT) + ' (no job files found in)', jobdir)

    for jobfile in jobfiles:
        numseeds = sum(1 for line in open(jobfile))
        sbatch_cmd.append('sbatch {} --array=1-{}:{} run_jobarray.sh -e {} -f {}'
                          .format(' '.join(sbatch_arg),numseeds, args.sims_per_task, exec, jobfile))

    if args.debug:
        print('sbatch args:', ' '.join(sbatch_arg))
        print('environment:', ' '.join(sbatch_env))
    return sbatch_cmd, sbatch_env




def main():
    project_name = 'DMRG'
    args = parse(project_name=project_name)
    sbatch_cmd, sbatch_env = generate_sbatch_commands(project_name=project_name,args=args)

    if not args.dryrun:
        with open('job_report.txt','a') as f:
            f.write('\n------------------------------------\n')
            f.write('Running sbatch - {}\n'.format(datetime.now().strftime("%Y-%m-%dT%H.%M.%S")))
            f.flush()
            run('git log -1', sbatch_env, args, shell=False, stdout=f, stderr=f)

    for cmd in sbatch_cmd:
        run(cmd=cmd, env=sbatch_env, args=args)

if __name__ == "__main__":
    main()
