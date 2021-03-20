#!/bin/bash
#SBATCH -p blade
#SBATCH -o slurmlog/count-%a.out
#SBATCH -e slurmlog/count-%a.err
#SBATCH -a 0-15
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

TOTALPARTS=16

hostname

source ~/.bashrc
module load hhvm

mkdir -p $OutDIR

hhvm findprob.php ${TOTALPARTS} ${SLURM_ARRAY_TASK_ID} 
#hhvm findprob.php ${TOTALPARTS} ${SLURM_ARRAY_TASK_ID} 0 libdef.austro_krypto.txt
#hhvm findprob.php ${TOTALPARTS} ${SLURM_ARRAY_TASK_ID} 0 libdef.NFZ.txt
exit

