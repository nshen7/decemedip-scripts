#!/bin/bash
#SBATCH --time=6:00:00        
#SBATCH --account=st-kdkortha-1   
#SBATCH --nodes=2
#SBATCH --ntasks=8            
#SBATCH --mem=64G               
#SBATCH --job-name=01_trim_fastq_files
#SBATCH --output=error_and_output/01_trim_fastq_files_output.txt    
#SBATCH --error=error_and_output/01_trim_fastq_files_error.txt      
#SBATCH --mail-user=ning.shen@stat.ubc.ca  
#SBATCH --mail-type=ALL        
##################################################################################################

# Navigate to the job's working directory
cd /scratch/st-kdkortha-1/nshen7/cfMeDIP_deconv/cfMeDIP-deconv-experiments/data/raw/Berchuck2022_LuCap_PDX_MeDIP/fastq
 
# Load conda environment
source ~/.bashrc
 
# Activate conda environment
conda activate /arc/project/st-kdkortha-1/nshen7/conda-envs/conda-env-cfMeDIP-deconv
 
# Add your commands here
trim_galore --paired --gzip ${line}_1.fq.gz ${line}_2.fq.gz

# Deactivate environment
conda deactivate
