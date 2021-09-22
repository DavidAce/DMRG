import argparse
import subprocess
import os
from pathlib import Path

parser = argparse.ArgumentParser(description='Gathers failed simulations')
parser.add_argument('-l', '--logdir', type=str, help='Search for log files in directory', default='logs')
parser.add_argument('-o', '--outdir', type=str, help='Save results in directory', default='failed')
parser.add_argument('-J', '--jobname', type=str, help='Filter by jobname', default=None)
parser.add_argument('-S', '--start', type=str, help='Consider jobs started after this date', default=None)
parser.add_argument('-E', '--end', type=str, help='Consider jobs that ended before this date', default=None)
parser.add_argument('-f', '--failfile', type=str, help='Save failed job list to this filename', default='failed_jobs.txt')
parser.add_argument('-r', '--resfile', type=str, help='Save resumable job list to this filename', default='resume.job')
parser.add_argument('-u', '--user', type=str, help='Call sacct for this user', default=os.getlogin())


args = parser.parse_args()


if args.outdir:
    Path(args.outdir).mkdir(parents=True, exist_ok=True)



with open("{}/{}".format(args.outdir, args.failfile), "w") as output:
    env=dict(os.environ, SACCT_FORMAT="jobid%16,jobname%16,user,state%14,maxrss%10,exitcode%10,account%20,cluster%20")
    sacct_command = ["sacct","-X"]
    sacct_command.extend(["-u", args.user])
    if args.jobname:
        sacct_command.extend(['--name',args.jobname])
    if args.start:
        sacct_command.extend(["-S", args.start])
    if args.end:
        sacct_command.extend(["-E", args.end])


    col_command = ['column', '-t']
    tr_command = ['tr', '-s', ' ']
    egrep_command = ['egrep', 'OUT|FAI|CAN']

    sacct = subprocess.Popen(sacct_command, stdout=subprocess.PIPE)
    col = subprocess.Popen(col_command, stdin=sacct.stdout, stdout=subprocess.PIPE, shell=False, encoding='utf-8')
    tr = subprocess.Popen(tr_command, stdin=col.stdout, stdout=subprocess.PIPE, shell=False, encoding='utf-8')
    egrep = subprocess.Popen(egrep_command, stdin=tr.stdout, stdout=subprocess.PIPE, shell=False, encoding='utf-8')

    sacct.stdout.close()
    col.stdout.close()
    tr.stdout.close()
    out, _ = egrep.communicate()
    output.write(out)

# Gather job-id's that contain failed jobs. Remember, 1 job-id corresponds to a job-array
# which may contain many simulations (seeds)
jobids = []
with open("{}/{}".format(args.outdir, args.failfile), "r") as output:
    for line in output:
        jobids.append(line.split(' ',1)[0])


with open("{}/{}".format(args.outdir,args.resfile), "w") as resfile:
    print('Creating resume file:', resfile)
count = 0
for file in os.listdir(args.logdir):
    filepath = "{}/{}".format(args.logdir,file)
    if not ".out" in file:
        continue
    if not os.path.isfile(filepath):
        continue
    if not any("{}.out".format(id) in file for id in jobids):
        continue

    # Found a logfile that contains at least one failed simulation
    with open(filepath, "r") as log, open("{}/{}".format(args.outdir,args.resfile), "a") as resfile:
        cfgline = ''
        for line in log:
            if not line.startswith(('CONFIG', 'JOB', 'EXIT')):
                continue
            keyval = [item.strip() for item in line.split(':',1)]
            if keyval[0] == 'CONFIG LINE' or keyval[0] == 'JOB FILE LINE(S)':
                cfgline = keyval[1]
            if keyval[0] == "EXIT CODE" and keyval[1] != "0":
                if len(cfgline) == 0:
                    raise Exception('cfgline is empty')
                print("Found failed simulation | exit {:>3} | {}".format(keyval[1],cfgline))
                resfile.write("{}\n".format(cfgline))
                count = count + 1

print('Found {} jobids with failures'.format(len(jobids)))
print('Found {} failed simulations'.format(count))