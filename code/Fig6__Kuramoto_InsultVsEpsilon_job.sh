#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=20:00:00
#SBATCH --job-name=Fig6__Kuramoto_InsultVsEpsilon_job
#SBATCH --output=Fig6__Kuramoto_InsulteVsEpsilon_output_%j.txt
#SBATCH --mail-type=FAIL

cd $SLURM_SUBMIT_DIR

./Fig6__Kuramoto_InsultVsEpsilon
