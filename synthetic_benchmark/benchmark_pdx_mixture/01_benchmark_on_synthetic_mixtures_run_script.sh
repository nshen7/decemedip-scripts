#!/bin/bash

cd /scratch/st-kdkortha-1/nshen7/cfMeDIP_deconv/cfMeDIP-deconv-experiments/code/synthetic_benchmark_new/benchmark_pdx_mixture

# declare -a coverages=(1)
declare -a coverages=(1 2)
nc=${#coverages[@]}

# declare -a proportions=(0.5)
declare -a proportions=(0.01 0.05 0.10 0.15 0.20 0.25 0.3 0.35 0.40 0.45 0.5)
np=${#proportions[@]}

# declare -a seeds=(2)
declare -a seeds=(1 2 3 4 5 6 7 8 9 10)
ns=${#seeds[@]}

for (( i=0; i<${nc}; i++ )); do
  for (( j=0; j<${np}; j++ )); do
    for (( k=0; k<${ns}; k++ )); do
      echo ${coverages[$i]} 
      echo ${proportions[$j]} 
      echo ${seeds[$k]} 
      my_command="sbatch --export=ALL,coverage=${coverages[$i]},proportion=${proportions[$j]},seed=${seeds[$k]} 01_benchmark_on_synthetic_mixtures.sh"
      eval $my_command
    done
  done
done