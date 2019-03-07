#!/bin/bash
#SBATCH --job-name=lfa
#SBATCH --account=rrg-hsn
#SBATCH --time=23:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=240000M
#SBATCH --mem-per-cpu=120000M
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
module load nixpkgs/16.09  gcc/7.3.0 r/3.5.2

echo "running analysis"

Rscript -e 'library(rmarkdown); rmarkdown::render("1KGP_lfa.Rmd", "all")'
