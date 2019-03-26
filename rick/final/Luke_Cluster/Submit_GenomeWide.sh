#!/bin/bash

execPath="/lb/project/gravel/luke_projects/1000Genomes"

for f in `seq 1 22`
do
echo -e " #!/bin/bash\n \
\n \
#PBS -l pmem=7700m\n \
#PBS -l nodes=1:ppn=16\n \
#PBS -l walltime=24:00:00\n \
#PBS -N Chrom_${f}\n \
\n \
export MUGQIC_INSTALL_HOME=/cvmfs/soft.mugqic/CentOS6
\n \
module use ${MUGQIC_INSTALL_HOME}/modulefiles
\n \
module load mugqic/R_Bioconductor/3.5.0_3.7

Rscript ${execPath}/RickRegression/3_fit_model_script.R \
${execPath}/RickRegression/ \
${execPath}/plink/ \
${f}" \
 > ${execPath}/RickRegression/job${f}_LM.pbs
qsub ${execPath}/RickRegression/job${f}_LM.pbs

done

exit
