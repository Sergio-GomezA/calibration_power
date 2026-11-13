# Calibration model run for whole data 1 model
#!/bin/bash
#$ -N clt
#$ -wd /exports/eddie/scratch/s2441782/calibration_power/
#$ -o /exports/eddie/scratch/s2441782/calibration_power/jobfiles/
#$ -e /exports/eddie/scratch/s2441782/calibration_power/jobfiles/
##$ -l h_rt=4:00:0,h_vmem=16G
##$ -l amd=true,intel=true
#$ -pe sharedmem 12
#$ -M s2441782@ed.ac.uk
#$ -m bea
##$ -t 25,75,125

# Initialise modules
source /etc/profile.d/modules.sh

# Load R
module load R/4.5

# Run resolution code
Rscript 8_LT_models/0_lt_main.R $SGE_TASK_ID TRUE 2190 batchLT_test FALSE 1
# file name/ day id/ recreate files / 
# days in traning samp / batch folder / save models / anomaly perc exclusion