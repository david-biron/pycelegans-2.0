#!/bin/sh
#SBATCH --exclusive
#SBATCH --constraint=ib
#SBATCH --nodes=1
source /etc/profile
module load pycelegans/0.61
echo "start time: `date`"
preprocessajp.py --ajpdir ./data
background.py --windowsize 10000 --Nimg 100
echo "end   time: `date`"