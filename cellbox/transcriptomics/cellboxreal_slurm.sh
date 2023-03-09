#!/bin/bash
#SBATCH -A emsls60202
#SBATCH -t 1-0
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -o cellboxreal_out
#SBATCH -e cellboxreal_err
#SBATCH --mail-user=liam.mackey@pnnl.gov
#SBATCH --mail-type END

module purge
module load pnnl_env


export PATH=/home/mack378/miniconda37/envs/py36/bin:$PATH

cd /home/mack378/CellBox

python scripts/main.py -config=configs/cellboxreal_config.json -i=121522
