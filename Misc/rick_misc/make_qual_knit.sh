#!/bin/bash
#SBATCH --job-name=make_qual
#SBATCH --account=rrg-hsn
#SBATCH --time=3:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem= 32400M

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
module load nixpkgs/16.09  gcc/7.3.0 r/3.5.2

echo "running analysis"

Rscript -e 'library(rmarkdown); rmarkdown::render("make_quality_table.Rmd", "all")'
