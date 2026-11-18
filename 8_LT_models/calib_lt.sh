# Calibration model run for whole data 1 model
#!/bin/bash
#$ -N LT1y3m10s
#$ -wd /exports/eddie/scratch/s2441782/calibration_power/
#$ -o /exports/eddie/scratch/s2441782/calibration_power/jobfiles/
#$ -e /exports/eddie/scratch/s2441782/calibration_power/jobfiles/
#$ -l h_rt=24:00:0,h_vmem=16G
#$ -pe sharedmem 48
#$ -M s2441782@ed.ac.uk
#$ -m bea


# Initialise modules
source /etc/profile.d/modules.sh

# Load R
module load R/4.5

# Run resolution code
Rscript 8_LT_models/0_lt_main.R 1 TRUE 365 90 10 batchLT FALSE 1
# file name/ day id/ recreate files / 
# days in traning samp / batch folder / save models / anomaly perc exclusion