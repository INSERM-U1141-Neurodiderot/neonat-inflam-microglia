#!/bin/bash
#SBATCH -J scpa_sub
#SBATCH -c 8
#SBATCH --mem=160G
#SBATCH --output=/home/adufour/work/logs/SCPA.log

cd /home/adufour/work
source /home/adufour/.bashrc
conda activate singlecell
Rscript /home/adufour/work/SCPA_sub.R