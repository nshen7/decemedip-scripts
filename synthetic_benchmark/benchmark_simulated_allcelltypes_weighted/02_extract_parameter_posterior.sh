#!/bin/bash

#SBATCH --time=3:00:00
#SBATCH --account=st-kdkortha-1
#SBATCH --nodes=1
#SBATCH --ntasks=6
#SBATCH --mem=72G
#SBATCH --job-name=01
#SBATCH -e error_and_output/02_extract_parameter_posterior_refv4.err
#SBATCH -o error_and_output/02_extract_parameter_posterior_refv4.out 
#SBATCH --mail-user=ning.shen@stat.ubc.ca
#SBATCH --mail-type=FAIL,TIME_LIMIT

################################################################################

module load apptainer

cd /scratch/st-kdkortha-1/nshen7/cfMeDIP_deconv/cfMeDIP-deconv-experiments/
apptainer exec \
-B /scratch/st-kdkortha-1/nshen7 \
-B /home/nshen7 \
/arc/project/st-kdkortha-1/nshen7/apptainers/bioc_cfmedip_deconv.sif \
Rscript code/synthetic_benchmark_new/benchmark_simulated_allcelltypes_weighted/02_extract_parameter_posterior_refv4.R 