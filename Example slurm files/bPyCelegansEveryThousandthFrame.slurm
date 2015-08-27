#!/bin/sh

#SBATCH --exclusive
#SBATCH --constraint=ib
#SBATCH --nodes=4

source /etc/profile
module load pycelegans/0.61
module load mpi4py

echo "start time: `date`"

export PATH=/scratch/midway/goessling/trunk/src:$PATH

mpirun -np 32 processcelegans.py --bgimage background_max_window_0000000_0010000.jpg --framerange 0 10000 --doviz yes  --imgdir ./input --thresh 140 --threshdiff 70 --filter 999.jpg --angle-step-size 4 --detmethod 1
mpirun -np 32 processcelegans.py --bgimage background_max_window_0010000_0020000.jpg --framerange 10000 20000 --doviz yes  --imgdir ./input --thresh 140 --threshdiff 70 --filter 999.jpg --angle-step-size 4 --detmethod 1
mpirun -np 32 processcelegans.py --bgimage background_max_window_0020000_0030000.jpg --framerange 20000 30000 --doviz yes  --imgdir ./input --thresh 140 --threshdiff 70 --filter 999.jpg --angle-step-size 4 --detmethod 1
mpirun -np 32 processcelegans.py --bgimage background_max_window_0030000_0040000.jpg --framerange 30000 40000 --doviz yes  --imgdir ./input --thresh 140 --threshdiff 70 --filter 999.jpg --angle-step-size 4 --detmethod 1
mpirun -np 32 processcelegans.py --bgimage background_max_window_0040000_0050000.jpg --framerange 40000 50000 --doviz yes  --imgdir ./input --thresh 140 --threshdiff 70 --filter 999.jpg --angle-step-size 4 --detmethod 1
mpirun -np 32 processcelegans.py --bgimage background_max_window_0050000_0060000.jpg --framerange 50000 60000 --doviz yes  --imgdir ./input --thresh 140 --threshdiff 70 --filter 999.jpg --angle-step-size 4 --detmethod 1
mpirun -np 32 processcelegans.py --bgimage background_max_window_0060000_0070000.jpg --framerange 60000 70000 --doviz yes  --imgdir ./input --thresh 140 --threshdiff 70 --filter 999.jpg --angle-step-size 4 --detmethod 1
mpirun -np 32 processcelegans.py --bgimage background_max_window_0080000_0090000.jpg --framerange 80000 90000 --doviz yes  --imgdir ./input --thresh 140 --threshdiff 70 --filter 999.jpg --angle-step-size 4 --detmethod 1
mpirun -np 32 processcelegans.py --bgimage background_max_window_0090000_0100000.jpg --framerange 90000 100000 --doviz yes  --imgdir ./input --thresh 140 --threshdiff 70 --filter 999.jpg --angle-step-size 4 --detmethod 1
mpirun -np 32 processcelegans.py --bgimage background_max_window_0100000_0110000.jpg --framerange 100000 110000 --doviz yes  --imgdir ./input --thresh 140 --threshdiff 70 --filter 999.jpg --angle-step-size 4 --detmethod 1
mpirun -np 32 processcelegans.py --bgimage background_max_window_0110000_0120000.jpg --framerange 110000 120000 --doviz yes  --imgdir ./input --thresh 140 --threshdiff 70 --filter 999.jpg --angle-step-size 4 --detmethod 1
mpirun -np 32 processcelegans.py --bgimage background_max_window_0120000_0130000.jpg --framerange 120000 130000 --doviz yes  --imgdir ./input --thresh 140 --threshdiff 70 --filter 999.jpg --angle-step-size 4 --detmethod 1
mpirun -np 32 processcelegans.py --bgimage background_max_window_0130000_0140000.jpg --framerange 130000 140000 --doviz yes  --imgdir ./input --thresh 140 --threshdiff 70 --filter 999.jpg --angle-step-size 4 --detmethod 1
mpirun -np 32 processcelegans.py --bgimage background_max_window_0140000_0150000.jpg --framerange 140000 150000 --doviz yes  --imgdir ./input --thresh 140 --threshdiff 70 --filter 999.jpg --angle-step-size 4 --detmethod 1
mpirun -np 32 processcelegans.py --bgimage background_max_window_0150000_0160000.jpg --framerange 150000 160000 --doviz yes  --imgdir ./input --thresh 140 --threshdiff 70 --filter 999.jpg --angle-step-size 4 --detmethod 1
mpirun -np 32 processcelegans.py --bgimage background_max_window_0160000_0170000.jpg --framerange 160000 170000 --doviz yes  --imgdir ./input --thresh 140 --threshdiff 70 --filter 999.jpg --angle-step-size 4 --detmethod 1
mpirun -np 32 processcelegans.py --bgimage background_max_window_0160000_0170000.jpg --framerange 170000 180000 --doviz yes  --imgdir ./input --thresh 140 --threshdiff 70 --filter 999.jpg --angle-step-size 4 --detmethod 1

collectoutput.py scratch

echo "end   time: `date`"
