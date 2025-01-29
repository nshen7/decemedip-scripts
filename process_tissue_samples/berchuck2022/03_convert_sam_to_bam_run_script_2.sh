#!/bin/bash

# Navigate to the job's working directory
cd /scratch/st-kdkortha-1/nshen7/cfMeDIP_deconv/cfMeDIP-deconv-experiments/code/process_tissue_samples/berchuck2022/
 
while read line; do
	my_command="sbatch --export=ALL,line=${line} 03_convert_sam_to_bam.sh"
  	eval $my_command
done < ../../../data/raw/Berchuck2022_LuCap_PDX_MeDIP/sample_name_list_2.txt

