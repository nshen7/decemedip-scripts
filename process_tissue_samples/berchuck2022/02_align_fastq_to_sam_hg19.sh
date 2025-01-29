#!/bin/bash
#SBATCH --time=10:00:00        
#SBATCH --account=st-kdkortha-1   
#SBATCH --nodes=1
#SBATCH --ntasks=6          
#SBATCH --mem=64G               
#SBATCH --job-name=02   # Specify the job name
#SBATCH --output=error_and_output/02_align_fastq_to_sam_hg19_output.txt    
#SBATCH --error=error_and_output/02_align_fastq_to_sam_hg19_error.txt      
#SBATCH --mail-user=ning.shen@stat.ubc.ca  
#SBATCH --mail-type=ALL        
 
##################################################################################################
## Navigate to the job's working directory
cd /scratch/st-kdkortha-1/nshen7/cfMeDIP_deconv/cfMeDIP-deconv-experiments/data/raw/Berchuck2022_LuCap_PDX_MeDIP/fastq
mkdir -p ../sam/
 
## Load conda environment
source ~/.bashrc
 
## Activate conda environment
conda activate /arc/project/st-kdkortha-1/nshen7/conda-envs/conda-env-cfMeDIP-deconv

## (hg19 reference genome downloaded beforehand. Script: code/download-ref-genome-hg19.sh)

## Align to SAM
bowtie2 -x ../../hg19/hg19 -p 20 -1 ${line}_1_val_1.fq.gz -2 ${line}_2_val_2.fq.gz -S ../sam/${line}.sam

# Deactivate environment
conda deactivate
