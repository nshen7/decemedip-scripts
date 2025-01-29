#!/bin/bash
#SBATCH --time=6:00:00        
#SBATCH --account=st-kdkortha-1   
#SBATCH --nodes=1
#SBATCH --ntasks=4           
#SBATCH --mem=64G               
#SBATCH --job-name=00_process_berchuck2022_metadata
#SBATCH --output=error_and_output/00_process_berchuck2022_metadata_output.txt    
#SBATCH --error=error_and_output/00_process_berchuck2022_metadata_error.txt      
#SBATCH --mail-user=ning.shen@stat.ubc.ca  
#SBATCH --mail-type=ALL        
##################################################################################################

module load apptainer

cd /scratch/st-kdkortha-1/nshen7/cfMeDIP_deconv/cfMeDIP-deconv-experiments/
apptainer exec \
-B /scratch/st-kdkortha-1/nshen7 \
-B /home/nshen7 \
/arc/project/st-kdkortha-1/nshen7/apptainers/bioc_cfmedip_deconv.sif \
Rscript code/process_tissue_samples/berchuck2022/00_process_berchuck2022_metadata.R
