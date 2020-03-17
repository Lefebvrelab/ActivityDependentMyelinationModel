#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=20:00:00
#SBATCH --job-name=Fig7__Epsilon0p2_sr1p0_job
#SBATCH --output=Fig7__Epsilon0p2_sr1p0_output_%j.txt
#SBATCH --mail-type=FAIL

cd $SLURM_SUBMIT_DIR

./Fig7__Epsilon0p2_sr1p0

