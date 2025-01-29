#!/bin/bash

#SBATCH --time=5:00:00
#SBATCH --account=st-kdkortha-1
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=64G
#SBATCH --job-name=00_process_shen2018_metadata
#SBATCH -e error_and_output/00_process_shen2018_metadata.err
#SBATCH -o error_and_output/00_process_shen2018_metadata.out 
#SBATCH --mail-user=ning.shen@stat.ubc.ca
#SBATCH --mail-type=ALL
 
################################################################################

module load apptainer

cd /scratch/st-kdkortha-1/nshen7/cfMeDIP_deconv/cfMeDIP-deconv-experiments/
apptainer exec \
-B /scratch/st-kdkortha-1/nshen7 \
-B /home/nshen7 \
/arc/project/st-kdkortha-1/nshen7/apptainers/bioc_cfmedip_deconv.sif \
Rscript code/case_studies/shen2018_new/00_process_shen2018_metadata.R

