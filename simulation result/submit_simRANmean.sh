#!/bin/bash
#SBATCH --job-name=bin2normSimRanmean         # Job name
#SBATCH --account=def-ubcxzh           # Your Compute Canada account
#SBATCH --cpus-per-task=10              # Number of CPU cores per task
#SBATCH --mem=100G                     
#SBATCH -t 7-00:00:00
#SBATCH --output=logs/job_%A_%a.out    # Standard output log
#SBATCH --error=logs/job_%A_%a.err     # Standard error log

# Load R module
module load r/4.4.0



Rscript simulationRAN200101A.R
Rscript simulationRAN200101B.R
Rscript simulationRAN200101C.R
Rscript simulationRAN200105A.R
Rscript simulationRAN200105B.R
Rscript simulationRAN200105C.R
Rscript simulationRAN200109A.R
Rscript simulationRAN200109B.R
Rscript simulationRAN200109C.R
Rscript simulationRAN201001A.R
Rscript simulationRAN201001B.R
Rscript simulationRAN201001C.R
Rscript simulationRAN201005A.R
Rscript simulationRAN201005B.R
Rscript simulationRAN201005C.R
Rscript simulationRAN201009A.R
Rscript simulationRAN201009B.R
Rscript simulationRAN201009C.R

Rscript simulationRAN500101A.R
Rscript simulationRAN500101B.R
Rscript simulationRAN500101C.R
Rscript simulationRAN500105A.R
Rscript simulationRAN500105B.R
Rscript simulationRAN500105C.R
Rscript simulationRAN500109A.R
Rscript simulationRAN500109B.R
Rscript simulationRAN500109C.R
Rscript simulationRAN501001A.R
Rscript simulationRAN501001B.R
Rscript simulationRAN501001C.R
Rscript simulationRAN501005A.R
Rscript simulationRAN501005B.R
Rscript simulationRAN501005C.R
Rscript simulationRAN501009A.R
Rscript simulationRAN501009B.R
Rscript simulationRAN501009C.R













