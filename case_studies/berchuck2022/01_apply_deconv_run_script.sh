#!/bin/bash

cd /scratch/st-kdkortha-1/nshen7/cfMeDIP_deconv/cfMeDIP-deconv-experiments/code/case_studies/berchuck2022/

declare -a indexs=(`seq 1 32`)
nidx=${#indexs[@]}

for (( i=0; i<${nidx}; i++ )); do
  echo ${indexs[$i]}
  my_command="sbatch --export=ALL,index=${indexs[$i]} 01_apply_deconv.sh"
  eval $my_command
done