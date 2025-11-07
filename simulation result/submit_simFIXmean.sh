#!/bin/bash
#SBATCH --job-name=bin2normSimFixmean         # Job name
#SBATCH --account=def-ubcxzh           # Your Compute Canada account
#SBATCH --cpus-per-task=10              # Number of CPU cores per task
#SBATCH --mem=100G                     
#SBATCH -t 7-00:00:00
#SBATCH --output=logs/job_%A_%a.out    # Standard output log
#SBATCH --error=logs/job_%A_%a.err     # Standard error log

# Load R module
module load r/4.4.0



Rscript simulationFIX2001A.R
Rscript simulationFIX2001B.R
Rscript simulationFIX2001C.R
Rscript simulationFIX5001A.R
Rscript simulationFIX5001B.R
Rscript simulationFIX5001C.R
Rscript simulationFIX2010A.R
Rscript simulationFIX2010B.R
Rscript simulationFIX2010C.R
Rscript simulationFIX5010A.R
Rscript simulationFIX5010B.R
Rscript simulationFIX5010C.R