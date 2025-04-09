#!/bin/bash
#SBATCH --job-name=test01
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=100M

#SBATCH --output=slurm_%j.out
#SBATCH --error=slurm_%j.err
#SBATCH --time=0-01:00:00
#SBATCH --partition=main

#SBATCH --mail-user=noah01@physik.fu-berlin.de
#SBATCH --mail-type=END

scontrol show job $SLURM_JOBID

module purge
module load Julia/1.10.59
julia -t ${SLURM_CPUS_PER_TASK} Tpmfrg_xyz.jl