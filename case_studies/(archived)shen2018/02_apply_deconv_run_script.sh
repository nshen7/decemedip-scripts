#!/bin/bash

cd /scratch/st-kdkortha-1/nshen7/cfMeDIP_deconv/cfMeDIP-deconv-experiments/code/case_studies/shen2018/

declare -a groups=('AML' 'BL' 'BRCA' 'Control' 'CRC' 'LUC' 'PDAC' 'RCC')
ngr=${#groups[@]}

for (( i=0; i<${ngr}; i++ )); do
  echo ${groups[$i]}
  my_command="sbatch --export=ALL,group=${groups[$i]} 02_apply_deconv.sh"
  eval $my_command
done