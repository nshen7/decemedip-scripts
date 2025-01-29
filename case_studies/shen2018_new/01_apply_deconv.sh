#!/bin/bash

#SBATCH --time=3:00:00
#SBATCH --account=st-kdkortha-1
#SBATCH --nodes=1
#SBATCH --ntasks=6
#SBATCH --mem=64G
#SBATCH --job-name=01
#SBATCH -e error_and_output/01_apply_deconv.err
#SBATCH -o error_and_output/01_apply_deconv.out 
#SBATCH --mail-user=ning.shen@stat.ubc.ca
#SBATCH --mail-type=FAIL,TIME_LIMIT

################################################################################

module load apptainer

cd /scratch/st-kdkortha-1/nshen7/cfMeDIP_deconv/cfMeDIP-deconv-experiments/
apptainer exec \
-B /scratch/st-kdkortha-1/nshen7 \
-B /home/nshen7 \
/arc/project/st-kdkortha-1/nshen7/apptainers/bioc_cfmedip_deconv.sif \
Rscript code/case_studies/shen2018_new/01_apply_deconv.R -i ${index}
