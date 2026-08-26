# Coverage check for the simulation study
#!/bin/bash
#$ -N cover_sim
#$ -wd /exports/eddie/scratch/s2441782/calibration_power/
#$ -o /exports/eddie/scratch/s2441782/calibration_power/jobfiles/
#$ -e /exports/eddie/scratch/s2441782/calibration_power/jobfiles/
#$ -l h_rt=4:00:0,h_vmem=16G,h_rss=16G
#$ -pe sharedmem 8
#$ -M s2441782@ed.ac.uk
#$ -m bea
#$ -t 1-150

# Initialise modules
source /etc/profile.d/modules.sh

# Load R
module load R/4.5

# Run resolution code
Rscript 5_extra/0_spat_aggr_elexon.R $SGE_TASK_ID 0.5
# nloc / n simulations / oos percentage
