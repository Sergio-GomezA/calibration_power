# Calibration model run for 1D
#!/bin/bash
#$ -N c12l0A0.5
#$ -wd /exports/eddie/scratch/s2441782/calibration_power/
#$ -o /exports/eddie/scratch/s2441782/calibration_power/jobfiles/
#$ -e /exports/eddie/scratch/s2441782/calibration_power/jobfiles/
##$ -l h_rt=4:00:0,h_vmem=16G
##$ -l amd=true,intel=true
#$ -pe sharedmem 8
#$ -M s2441782@ed.ac.uk
#$ -m bea
#$ -t 25,75,125

# Initialise modules
source /etc/profile.d/modules.sh

# Load R
module load R/4.5

# Run resolution code
Rscript 4_block_code/0_main_calib_valid.R $SGE_TASK_ID 50 TRUE TRUE 120 batchY25d150_v12 FALSE 0.5
# file name/ day id/ mesh edge length / recreate files / 
# rerun st model / days in traning samp / batch folder / save models / anomaly perc exclusion


