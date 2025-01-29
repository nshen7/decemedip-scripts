#!/bin/bash

#SBATCH --time=10:00:00
#SBATCH --account=st-kdkortha-1
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=64G
#SBATCH --job-name=04
#SBATCH -e error_and_output/04_get_read_counts_refv2.err
#SBATCH -o error_and_output/04_get_read_counts_refv2.out 
#SBATCH --mail-user=ning.shen@stat.ubc.ca
#SBATCH --mail-type=ALL
 
################################################################################

module load apptainer

cd /scratch/st-kdkortha-1/nshen7/cfMeDIP_deconv/cfMeDIP-deconv-experiments/
apptainer exec \
-B /scratch/st-kdkortha-1/nshen7 \
-B /home/nshen7 \
/arc/project/st-kdkortha-1/nshen7/apptainers/bioc_cfmedip_deconv.sif \
Rscript code/process_tissue_samples/berchuck2022/04_get_read_counts_refv2.R -b 1

