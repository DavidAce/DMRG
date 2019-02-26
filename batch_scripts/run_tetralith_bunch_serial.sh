#!/bin/bash



#SBATCH --job-name=DMRG
#SBATCH --time=0-15:00:00
#SBATCH --mem-per-cpu=2000
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-core=1
#SBATCH --kill-on-invalid-dep=yes

echo "Running on tetralith"
echo "$@"
exec=$1
bunch_filename="$2"

while read -r infile
do
    $exec $infile
done < $bunch_filename

# Explanation

# This proces will read the lines in bunch_filename and feed them to the executable as arguments, one by one.
# It will wait for the executable to terminate before going to the next iteration.