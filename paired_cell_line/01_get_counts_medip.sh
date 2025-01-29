#!/bin/bash

#SBATCH --time=3:00:00
#SBATCH --account=st-kdkortha-1
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem=128G
#SBATCH --job-name=01_get_counts_medip
#SBATCH -e sbatch_error_and_output/01_get_counts_medip.err
#SBATCH -o sbatch_error_and_output/01_get_counts_medip.out 
#SBATCH --mail-user=ning.shen@stat.ubc.ca
#SBATCH --mail-type=ALL
 
################################################################################

module load apptainer

cd /scratch/st-kdkortha-1/nshen7/cfMeDIP_deconv/cfMeDIP-deconv-experiments/
apptainer exec \
-B /scratch/st-kdkortha-1/nshen7 \
-B /home/nshen7 \
/arc/project/st-kdkortha-1/nshen7/apptainers/bioc_cfmedip_deconv.sif \
Rscript code/medip_vs_wgbs_ref_regions_v4/01_get_counts_medip.R -b 1
