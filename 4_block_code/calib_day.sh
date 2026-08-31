# Calibration model run for 1D
#!/bin/bash
#$ -N calval
#$ -wd /exports/eddie/scratch/s2441782/calibration_power/
#$ -o /exports/eddie/scratch/s2441782/calibration_power/jobfiles/
#$ -e /exports/eddie/scratch/s2441782/calibration_power/jobfiles/
##$ -l h_rt=4:00:0,h_vmem=16G
#$ -l amd=true,intel=true
#$ -pe sharedmem 8
#$ -M s2441782@ed.ac.uk
#$ -m bea
#$ -t 1-150

# Initialise modules
source /etc/profile.d/modules.sh

# Load R
module load R/4.5

# Run resolution code
Rscript 4_block_code/0_main_calib_valid.R $SGE_TASK_ID 50 TRUE TRUE 3 batchY25d150 FALSE 0.3
# file name/ day id/ mesh edge length / recreate files / 
# rerun st model / days in traning samp / batch folder / save models