#!/bin/sh
#SBATCH --nodes=3
#SBATCH --ntasks-per-node=45
#SBATCH --time=00:05:00
#SBATCH --job-name=spectrum_processing
#SBATCH --error=error_log.log
#SBATCH --output=spectrum_processing.log
#SBATCH --partition=debug
#SBATCH --mail-type=ALL
#SBATCH --mail-user kaustav_b@ph.iitr.ac.in

cd /home/kaustav_b_ph.iitr/Lunar-Mapping/Code
module load python/python3.10
conda activate conda_env
python test.py