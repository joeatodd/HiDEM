#!/bin/bash -l                                                                     

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH -J Calv
#SBATCH -p medium
#SBATCH -t 1:00:00
#SBATCH --account=project_2000881
#SBATCH -o /scratch/project_2000881/acrawford/JI_HiDEM/%J.out
#SBATCH -e /scratch/project_2000881/acrawford/JI_HiDEM/%J.err

echo "Running HiDEM"
export OMP_NUM_THREADS=1

mkdir -p ./Results ./work ./logs ./errlogs

srun -n 128 ./HiDEM
mv *err errlogs
mv *out logs
