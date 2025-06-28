#!/bin/bash

#SBATCH --mem=5GB
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --account=BIOL-SPECGEN-2018
#SBATCH --time=0-00:30:00

################################
# Load libraries and set paths #
################################
module load SLiM/5.0-GCC-13.3.0
module load msprime/1.0.1-foss-2020a-Python-3.8.2

slim simulations.slim
python3 ./Tree_to_fasta.py > simu.fasta

