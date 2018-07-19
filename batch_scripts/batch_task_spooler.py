import subprocess
import os
cwd = os.path.dirname(os.path.realpath(__file__))






command = 'ls'


process = subprocess.Popen(command.split(), stdout=subprocess.PIPE, stdin=subprocess.PIPE, shell=True, cwd=cwd)
output, error = process.communicate()

print(output, error)