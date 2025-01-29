#!/bin/bash

cd /scratch/st-kdkortha-1/nshen7/cfMeDIP_deconv/cfMeDIP-deconv-experiments/code/case_studies/shen2018_new/

declare -a indexs=(`seq 1 188`)
nidx=${#indexs[@]}

for (( i=0; i<${nidx}; i++ )); do
  echo ${indexs[$i]}
  my_command="sbatch --export=ALL,index=${indexs[$i]} 01_apply_deconv_refv4.sh"
  eval $my_command
done