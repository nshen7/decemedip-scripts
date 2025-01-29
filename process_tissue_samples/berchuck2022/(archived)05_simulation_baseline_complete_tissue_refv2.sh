#!/bin/bash

#SBATCH --time=12:00:00
#SBATCH --account=st-kdkortha-1
#SBATCH --nodes=1
#SBATCH --ntasks=6
#SBATCH --mem=96G
#SBATCH --job-name=05_tissue_refv2
#SBATCH -e error_and_output/05_simulation_baseline_complete_tissue_refv2.err
#SBATCH -o error_and_output/05_simulation_baseline_complete_tissue_refv2.out 
#SBATCH --mail-user=ning.shen@stat.ubc.ca
#SBATCH --mail-type=ALL

################################################################################

module load apptainer

cd /scratch/st-kdkortha-1/nshen7/cfMeDIP_deconv/cfMeDIP-deconv-experiments/
apptainer exec \
-B /scratch/st-kdkortha-1/nshen7 \
-B /home/nshen7 \
/arc/project/st-kdkortha-1/nshen7/apptainers/bioc_cfmedip_deconv.sif \
Rscript code/process_tissue_samples/berchuck2022/05_simulation_baseline_complete_tissue_refv2.R

