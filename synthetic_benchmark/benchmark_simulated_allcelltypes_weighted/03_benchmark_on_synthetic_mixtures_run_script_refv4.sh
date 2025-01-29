#!/bin/bash

cd /scratch/st-kdkortha-1/nshen7/cfMeDIP_deconv/cfMeDIP-deconv-experiments/code/synthetic_benchmark_new/benchmark_simulated_allcelltypes_weighted

declare -a coverages=(1)
# declare -a coverages=(1 2 3)
nc=${#coverages[@]}

declare -a proportions=(0.25)
# declare -a proportions=(0.01 0.05 0.10 0.15 0.20 0.25 0.3 0.35 0.40 0.45 0.5)
np=${#proportions[@]}

declare -a seeds=(10)
# declare -a seeds=(1 2 3 4 5 6 7 8 9 10)
ns=${#seeds[@]}

for (( i=0; i<${nc}; i++ )); do
  for (( j=0; j<${np}; j++ )); do
    for (( k=0; k<${ns}; k++ )); do
      echo "coverage=${coverages[$i]},proportion=${proportions[$j]},seed=${seeds[$k]}"
      my_command="sbatch --export=ALL,coverage=${coverages[$i]},proportion=${proportions[$j]},seed=${seeds[$k]} 03_benchmark_on_synthetic_mixtures_refv4.sh"
      eval $my_command
    done
  done
done
# done