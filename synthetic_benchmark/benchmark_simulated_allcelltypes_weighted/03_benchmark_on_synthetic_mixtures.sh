#!/bin/bash

#SBATCH --time=12:00:00
#SBATCH --account=st-kdkortha-1
#SBATCH --nodes=1
#SBATCH --ntasks=6
#SBATCH --mem=72G
#SBATCH --job-name=03
#SBATCH -e error_and_output/03_benchmark_on_synthetic_mixtures.err
#SBATCH -o error_and_output/03_benchmark_on_synthetic_mixtures.out 
#SBATCH --mail-user=ning.shen@stat.ubc.ca
#SBATCH --mail-type=FAIL,TIME_LIMIT

################################################################################

module load apptainer

cd /scratch/st-kdkortha-1/nshen7/cfMeDIP_deconv/cfMeDIP-deconv-experiments/
apptainer exec \
-B /scratch/st-kdkortha-1/nshen7 \
-B /home/nshen7 \
/arc/project/st-kdkortha-1/nshen7/apptainers/bioc_cfmedip_deconv.sif \
Rscript code/synthetic_benchmark_new/benchmark_simulated_allcelltypes_weighted/03_benchmark_on_synthetic_mixtures.R \
-c ${coverage} \
-p ${proportion} \
-s ${seed} 

