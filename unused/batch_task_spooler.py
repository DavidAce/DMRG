import subprocess
import os
cwd = os.path.dirname(os.path.realpath(__file__))


executable = './../cmake-build-release/DMRG++'
prefix = "mbl_"
suffix = ".cfg"
jobnum = range(0,5)


for num in jobnum:
    label  = prefix + str(num)
    command =  'tsp -L ' + label + ' ' +executable + ' ' + prefix + str(num) + suffix
    print (command)
    process = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True, cwd=cwd)
    output, error = process.communicate()
    print(output, error)



