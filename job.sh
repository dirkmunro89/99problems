#!/bin/sh

#SBATCH --job-name=NAME_X
#SBATCH --partition=compute
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=500MB
#SBATCH --account=research-3me-pme

module load 2022r2
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

source ../99problems_0/venv/bin/activate

python ../99problems_0/main.py PROBLEM_X
deactivate
