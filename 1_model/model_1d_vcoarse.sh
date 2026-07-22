# Calibration model run for 1D
#!/bin/bash
#$ -N calvcoarse
#$ -wd /exports/eddie/scratch/s2441782/calibration_power/
#$ -o /exports/eddie/scratch/s2441782/calibration_power/jobfiles/
#$ -e /exports/eddie/scratch/s2441782/calibration_power/jobfiles/
##$ -l h_rt=4:00:0,h_vmem=16G
#$ -pe sharedmem 12
#$ -M s2441782@ed.ac.uk
#$ -m bea
#$ -t 1-15

# Initialise modules
source /etc/profile.d/modules.sh

# Load R
module load R/4.5

# Run resolution code
Rscript 1_model/inlabru_1d.R $SGE_TASK_ID 50 TRUE TRUE 3 batch2025 TRUE
