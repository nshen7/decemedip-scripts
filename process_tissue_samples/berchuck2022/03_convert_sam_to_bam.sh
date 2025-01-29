#!/bin/bash
#SBATCH --time=10:00:00        
#SBATCH --account=st-kdkortha-1   
#SBATCH --nodes=2
#SBATCH --ntasks=6           
#SBATCH --mem=64G               
#SBATCH --job-name=03_convert_sam_to_bam   # Specify the job name
#SBATCH --output=error_and_output/03_convert_sam_to_bam_output.txt    
#SBATCH --error=error_and_output/03_convert_sam_to_bam_error.txt      
#SBATCH --mail-user
#SBATCH --mail-type=ALL        
 
##################################################################################################
## Navigate to the job's working directory
cd /scratch/st-kdkortha-1/nshen7/cfMeDIP_deconv/cfMeDIP-deconv-experiments/data/raw/Berchuck2022_LuCap_PDX_MeDIP/
mkdir -p bam
cd bam

## Load conda environment
source ~/.bashrc
 
## Activate conda environment
conda activate /arc/project/st-kdkortha-1/nshen7/conda-envs/conda-env-cfMeDIP-deconv

## (hg19 reference genome downloaded beforehand. Script: code/download-ref-genome-hg19.sh)
convertToBAM() {
	local acc=$1
	samtools view -@ 20 -bS "../sam/${acc}.sam" > "${acc}.bam"
    samtools sort -@ 20 "${acc}.bam" -o "${acc}_sorted.bam"
    samtools index -@ 20 "${acc}_sorted.bam"
}

convertToBAM ${line}

# Deactivate environment
conda deactivate
