#!/bin/sh

#SBATCH --exclusive
#SBATCH --constraint=ib
#SBATCH --nodes=4

source /etc/profile
module load pycelegans/0.61
module load mpi4py

echo "start time: `date`"

export PATH=/scratch/midway/goessling/trunk/src:$PATH

mpirun -np 32 processcelegans.py --bgimage background_max_window_0000000_0010000.jpg --framerange 0 10000 --doviz yes  --imgdir ./input --thresh 170 --threshdiff 30 --filter 999.jpg --angle-step-size 6 --detmethod 0
mpirun -np 32 processcelegans.py --bgimage background_max_window_0010000_0020000.jpg --framerange 10000 20000 --doviz yes  --imgdir ./input --thresh 170 --threshdiff 30 --filter 999.jpg --angle-step-size 6 --detmethod 0
mpirun -np 32 processcelegans.py --bgimage background_max_window_0020000_0030000.jpg --framerange 20000 30000 --doviz yes  --imgdir ./input --thresh 170 --threshdiff 30 --filter 999.jpg --angle-step-size 6 --detmethod 0
mpirun -np 32 processcelegans.py --bgimage background_max_window_0030000_0040000.jpg --framerange 30000 40000 --doviz yes  --imgdir ./input --thresh 170 --threshdiff 30 --filter 999.jpg --angle-step-size 6 --detmethod 0
mpirun -np 32 processcelegans.py --bgimage background_max_window_0040000_0050000.jpg --framerange 40000 50000 --doviz yes  --imgdir ./input --thresh 170 --threshdiff 30 --filter 999.jpg --angle-step-size 6 --detmethod 0
mpirun -np 32 processcelegans.py --bgimage background_max_window_0050000_0060000.jpg --framerange 50000 60000 --doviz yes  --imgdir ./input --thresh 170 --threshdiff 30 --filter 999.jpg --angle-step-size 6 --detmethod 0
mpirun -np 32 processcelegans.py --bgimage background_max_window_0060000_0070000.jpg --framerange 60000 70000 --doviz yes  --imgdir ./input --thresh 170 --threshdiff 30 --filter 999.jpg --angle-step-size 6 --detmethod 0
mpirun -np 32 processcelegans.py --bgimage background_max_window_0070000_0080000.jpg --framerange 70000 80000 --doviz yes  --imgdir ./input --thresh 170 --threshdiff 30 --filter 999.jpg --angle-step-size 6 --detmethod 0
mpirun -np 32 processcelegans.py --bgimage background_max_window_0080000_0090000.jpg --framerange 80000 90000 --doviz yes  --imgdir ./input --thresh 170 --threshdiff 30 --filter 999.jpg --angle-step-size 6 --detmethod 0
mpirun -np 32 processcelegans.py --bgimage background_max_window_0090000_0100000.jpg --framerange 90000 100000 --doviz yes  --imgdir ./input --thresh 170 --threshdiff 30 --filter 999.jpg --angle-step-size 6 --detmethod 0
mpirun -np 32 processcelegans.py --bgimage background_max_window_0100000_0110000.jpg --framerange 100000 110000 --doviz yes  --imgdir ./input --thresh 170 --threshdiff 30 --filter 999.jpg --angle-step-size 6 --detmethod 0
mpirun -np 32 processcelegans.py --bgimage background_max_window_0110000_0120000.jpg --framerange 110000 120000 --doviz yes  --imgdir ./input --thresh 170 --threshdiff 30 --filter 999.jpg --angle-step-size 6 --detmethod 0
mpirun -np 32 processcelegans.py --bgimage background_max_window_0120000_0130000.jpg --framerange 120000 130000 --doviz yes  --imgdir ./input --thresh 170 --threshdiff 30 --filter 999.jpg --angle-step-size 6 --detmethod 0
mpirun -np 32 processcelegans.py --bgimage background_max_window_0130000_0140000.jpg --framerange 130000 140000 --doviz yes  --imgdir ./input --thresh 170 --threshdiff 30 --filter 999.jpg --angle-step-size 6 --detmethod 0
mpirun -np 32 processcelegans.py --bgimage background_max_window_0140000_0150000.jpg --framerange 140000 150000 --doviz yes  --imgdir ./input --thresh 170 --threshdiff 30 --filter 999.jpg --angle-step-size 6 --detmethod 0
mpirun -np 32 processcelegans.py --bgimage background_max_window_0150000_0160000.jpg --framerange 150000 160000 --doviz yes  --imgdir ./input --thresh 170 --threshdiff 30 --filter 999.jpg --angle-step-size 6 --detmethod 0
mpirun -np 32 processcelegans.py --bgimage background_max_window_0160000_0170000.jpg --framerange 160000 170000 --doviz yes  --imgdir ./input --thresh 170 --threshdiff 30 --filter 999.jpg --angle-step-size 6 --detmethod 0
mpirun -np 32 processcelegans.py --bgimage background_max_window_0170000_0180000.jpg --framerange 170000 180000 --doviz yes  --imgdir ./input --thresh 170 --threshdiff 30 --filter 999.jpg --angle-step-size 6 --detmethod 0
mpirun -np 32 processcelegans.py --bgimage background_max_window_0180000_0190000.jpg --framerange 180000 190000 --doviz yes  --imgdir ./input --thresh 170 --threshdiff 30 --filter 999.jpg --angle-step-size 6 --detmethod 0
mpirun -np 32 processcelegans.py --bgimage background_max_window_0190000_0200000.jpg --framerange 190000 200000 --doviz yes  --imgdir ./input --thresh 170 --threshdiff 30 --filter 999.jpg --angle-step-size 6 --detmethod 0
mpirun -np 32 processcelegans.py --bgimage background_max_window_0200000_0210000.jpg --framerange 200000 210000 --doviz yes  --imgdir ./input --thresh 170 --threshdiff 30 --filter 999.jpg --angle-step-size 6 --detmethod 0
mpirun -np 32 processcelegans.py --bgimage background_max_window_0210000_0220000.jpg --framerange 210000 220000 --doviz yes  --imgdir ./input --thresh 170 --threshdiff 30 --filter 999.jpg --angle-step-size 6 --detmethod 0
mpirun -np 32 processcelegans.py --bgimage background_max_window_0220000_0230000.jpg --framerange 220000 230000 --doviz yes  --imgdir ./input --thresh 170 --threshdiff 30 --filter 999.jpg --angle-step-size 6 --detmethod 0
mpirun -np 32 processcelegans.py --bgimage background_max_window_0230000_0240000.jpg --framerange 230000 240000 --doviz yes  --imgdir ./input --thresh 170 --threshdiff 30 --filter 999.jpg --angle-step-size 6 --detmethod 0
mpirun -np 32 processcelegans.py --bgimage background_max_window_0240000_0250000.jpg --framerange 240000 250000 --doviz yes  --imgdir ./input --thresh 170 --threshdiff 30 --filter 999.jpg --angle-step-size 6 --detmethod 0
mpirun -np 32 processcelegans.py --bgimage background_max_window_0250000_0260000.jpg --framerange 250000 260000 --doviz yes  --imgdir ./input --thresh 170 --threshdiff 30 --filter 999.jpg --angle-step-size 6 --detmethod 0
mpirun -np 32 processcelegans.py --bgimage background_max_window_0260000_0270000.jpg --framerange 260000 270000 --doviz yes  --imgdir ./input --thresh 170 --threshdiff 30 --filter 999.jpg --angle-step-size 6 --detmethod 0
mpirun -np 32 processcelegans.py --bgimage background_max_window_0270000_0280000.jpg --framerange 270000 280000 --doviz yes  --imgdir ./input --thresh 170 --threshdiff 30 --filter 999.jpg --angle-step-size 6 --detmethod 0
mpirun -np 32 processcelegans.py --bgimage background_max_window_0280000_0290000.jpg --framerange 280000 290000 --doviz yes  --imgdir ./input --thresh 170 --threshdiff 30 --filter 999.jpg --angle-step-size 6 --detmethod 0
mpirun -np 32 processcelegans.py --bgimage background_max_window_0290000_0300000.jpg --framerange 290000 300000 --doviz yes  --imgdir ./input --thresh 170 --threshdiff 30 --filter 999.jpg --angle-step-size 6 --detmethod 0
mpirun -np 32 processcelegans.py --bgimage background_max_window_0300000_0310000.jpg --framerange 300000 310000 --doviz yes  --imgdir ./input --thresh 170 --threshdiff 30 --filter 999.jpg --angle-step-size 6 --detmethod 0
mpirun -np 32 processcelegans.py --bgimage background_max_window_0310000_0320000.jpg --framerange 310000 320000 --doviz yes  --imgdir ./input --thresh 170 --threshdiff 30 --filter 999.jpg --angle-step-size 6 --detmethod 0
mpirun -np 32 processcelegans.py --bgimage background_max_window_0320000_0330000.jpg --framerange 320000 330000 --doviz yes  --imgdir ./input --thresh 170 --threshdiff 30 --filter 999.jpg --angle-step-size 6 --detmethod 0
mpirun -np 32 processcelegans.py --bgimage background_max_window_0330000_0340000.jpg --framerange 330000 340000 --doviz yes  --imgdir ./input --thresh 170 --threshdiff 30 --filter 999.jpg --angle-step-size 6 --detmethod 0
mpirun -np 32 processcelegans.py --bgimage background_max_window_0340000_0350000.jpg --framerange 340000 360000 --doviz yes  --imgdir ./input --thresh 170 --threshdiff 30 --filter 999.jpg --angle-step-size 6 --detmethod 0

collectoutput.py scratch

echo "end   time: `date`"
