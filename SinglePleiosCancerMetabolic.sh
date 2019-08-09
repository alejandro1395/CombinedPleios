#!/bin/bash

#set the job name
#SBATCH --job-name=Haps_Onset_Immune
#SBATCH --cpus-per-per-task = 1
#SBATCH -o slurm.%j.out
#SBATCH -e slurm.%j.err
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --time=01:00:00

#run the application
#PATHS
OUTPUT=/homes/users/avalenzuela/scratch/PhD_EvoGenomics/1st_year/StratiPleios_Jul2019/Combinations/results/

module load Python

### Python script to count pleiotropies in haplotypes
mkdir -p ${OUTPUT}SinglePleios_CancerMetabolic/
for i in {10..60};
   do mkdir -p ${OUTPUT}SinglePleios_CancerMetabolic/Age_threeshold_${i}/;
	python3 SinglePleiosCancerMetabolic.py ${i} ${OUTPUT}SinglePleios_CancerMetabolic/Age_threeshold_${i}/
    done;
