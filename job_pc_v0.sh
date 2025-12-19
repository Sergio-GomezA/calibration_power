#!/bin/bash

#$ -N PCmodels
#$ -wd /exports/eddie/scratch/s2441782/calibration_power/
#$ -o /exports/eddie/scratch/s2441782/calibration_power/jobfiles/
#$ -e /exports/eddie/scratch/s2441782/calibration_power/jobfiles/
#$ -M s2441782@ed.ac.uk
#$ -m bea

# Initialise modules
source /etc/profile.d/modules.sh

# Load R
module load R/4.5

# Run resolution code
Rscript power_curve_models.R
