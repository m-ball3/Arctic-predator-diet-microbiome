#!/bin/bash
#SBATCH --account=coenv
#SBATCH --partition=cpu-g2 # Queue selection
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --mem=10G
#SBATCH --time=10:00:00 # Time limit hrs:min:sec\
#SBATCH --job-name=WADE003-16SP2-dada2_QAQC # Job name
#SBATCH --mail-type=ALL # Mail events (BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=mball3@uw.edu # Where to send mail\
#SBATCH --output=WADE003-16P2_dada2_QAQC.out # Standard output/error\
#SBATCH --chdir=/gscratch/coenv/mball3
#SBATCH --export=all
#SBATCH --error=WADE003-16P2_dada2_QAQC.err

module load apptainer

apptainer run --bind /gscratch /gscratch/coenv/containerstidyverse_4.0.1.sif Rscript /gscratch/coenv/mball3/WADE003-arctic-pred/scripts/dada2QAQC-arcticpred.R