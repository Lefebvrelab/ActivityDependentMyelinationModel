#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=20:00:00
#SBATCH --job-name=Fig5__Epsilon0p2_job
#SBATCH --output=Fig5__Epsilon0p2_output_%j.txt
#SBATCH --mail-type=FAIL

cd $SLURM_SUBMIT_DIR

./Fig5__Epsilon0p2
