#!/bin/bash

#SBATCH --time=12:00:00
#SBATCH --account=st-kdkortha-1
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=96G
#SBATCH --job-name=03
#SBATCH -e error_and_output/03_summarize_results.err
#SBATCH -o error_and_output/03_summarize_results.out 
#SBATCH --mail-user=ning.shen@stat.ubc.ca
#SBATCH --mail-type=ALL

################################################################################

module load apptainer

cd /scratch/st-kdkortha-1/nshen7/cfMeDIP_deconv/cfMeDIP-deconv-experiments/
apptainer exec \
-B /scratch/st-kdkortha-1/nshen7 \
-B /home/nshen7 \
/arc/project/st-kdkortha-1/nshen7/apptainers/bioc_cfmedip_deconv.sif \
Rscript code/case_studies/shen2018/03_summarize_results.R
