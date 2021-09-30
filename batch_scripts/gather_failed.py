import argparse
import subprocess
import os
from pathlib import Path
from copy import deepcopy
import time
parser = argparse.ArgumentParser(description='Gathers failed simulations')
parser.add_argument('-l', '--logdir', type=str, help='Search for log files in directory', default='logs')
parser.add_argument('-o', '--outdir', type=str, help='Save results in directory', default='failed')
parser.add_argument('-J', '--jobname', type=str, help='Filter by jobname', default=None)
parser.add_argument('-S', '--start', type=str, help='Consider jobs started after this date', default=None)
parser.add_argument('-E', '--end', type=str, help='Consider jobs that ended before this date', default=time.strftime('%Y-%m-%dT%H:%M:%S'))
parser.add_argument('-L', '--logscan', action='store_true', help='Scan through logs instead of using sacct')
parser.add_argument('-f', '--failfile', type=str, help='Save list of jobids with failed simulations to this file', default='failed_jobs.txt')
parser.add_argument('-r', '--resfile', type=str, help='Save results, i.e. a list of cfg-seed pairs, to this file', default='resume.job')
parser.add_argument('-u', '--user', type=str, help='Call sacct for this user', default=os.getlogin())


args = parser.parse_args()


if args.outdir:
    Path(args.outdir).mkdir(parents=True, exist_ok=True)

def getlines(file, linelist):
  return [x for i, x in enumerate(file) if i in map(int,linelist)]


# Gather job-id's that may contain failed jobs. Remember, 1 job-id can correspond to a job-array
# which may contain many simulations (seeds)
joblist = []
jobitem = {
    'jobid': None,
    'jobidraw': None,
    'jobname' : None,
    'exitcode': None,
    'outfile' : None,
}

with open("{}/{}".format(args.outdir, args.failfile), "w") as output:
    sacct_command = ["sacct","-X", "--parsable2", "--noheader", '--format=jobid,jobidraw,jobname,exitcode,state']
    if args.logscan:
        sacct_command.append("--state=failed,timeout,deadline,node_fail,completed")
    else:
        sacct_command.append("--state=failed,timeout,deadline,node_fail")

    sacct_command.extend(["-u", args.user])
    if args.jobname:
        sacct_command.extend(['--name',args.jobname])
    if args.start:
        sacct_command.extend(["-S", args.start])
    if args.end:
        sacct_command.extend(["-E", args.end])

    with subprocess.Popen(sacct_command,bufsize=1, shell=False, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, encoding='utf-8') as p:
        for out in p.stdout:
            print(out, end='')  # process line here
            if out and not out.isspace():
                output.write(out)
                line = out.split('|')
                item = {}
                item['jobid']    = line[0]
                item['jobidraw'] = line[1]
                item['jobname'] = line[2]
                item['exitcode'] = line[3]
                item['outfile'] = '{}/{}-{}.out'.format(args.logdir,line[2], line[0]) # This needs to match the sbatch file
                joblist.append(item)

        print("Command:", ' '.join(p.args))
        if p.returncode != 0 and p.returncode != None:
            raise subprocess.CalledProcessError(p.returncode, ' '.join(p.args))



with open("{}/{}".format(args.outdir,args.resfile), "w") as resfile:
    print('Creating resume file:', resfile)
count = 0
for jobitem in joblist:
    if not os.path.isfile(jobitem['outfile']):
        raise FileNotFoundError("File does not exist: {}".format(jobitem['outfile']))

    # Found a logfile that may contain a failed simulation (or not, if doing logscan)
    with open(jobitem['outfile'], "r") as log, open("{}/{}".format(args.outdir,args.resfile), "a") as resfile:
        cfgline = None
        sequence = None
        jobfile = None
        success = set([])
        jobid   = None
        for line in log:
            if not line.startswith(('JOB ID', 'JOB FILE', 'EXIT', 'TASK ID SEQUENCE')):
                continue
            keyval = [item.strip() for item in line.split(':',1)]
            if keyval[0] == 'JOB FILE':
                jobfile = keyval[1] # Contains the array jobarray chunk file.
            if keyval[0] == 'TASK ID SEQUENCE':
                sequence = set(map(int,keyval[1].split(' ')))
            if keyval[0] == 'JOB ID':
                jobid = int(keyval[1])
            if keyval[0] == "EXIT CODE" and keyval[1] == "0":
                success.add(jobid)

        # We can now compare the ids in the task id sequence that should have run,
        # to the ones that actually suceeded.
        failed = sequence - success

        # We can now extract the config lines that should have executed
        cfglines = getlines(open(jobfile,"r"), failed)
        for cfgline in cfglines:
            print("Found failed simulation | exit {:>3} | {}".format(keyval[1],cfgline))
            resfile.write("{}\n".format(cfgline))
            count = count + 1
#
#             if keyval[0] == 'CONFIG LINE' or keyval[0] == 'JOB FILE LINE(S)':
#                 cfgline = keyval[1]
#             if keyval[0] == "EXIT CODE" and keyval[1] != "0":
#                 if not cfgline:
#                     raise Exception('cfgline is empty')
#                 print("Found failed simulation | exit {:>3} | {}".format(keyval[1],cfgline))
#                 resfile.write("{}\n".format(cfgline))
#                 count = count + 1

print('Found {} jobids'.format(len(joblist)))
print('Found {} failed simulations'.format(count))