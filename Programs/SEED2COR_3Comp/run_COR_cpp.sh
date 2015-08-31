#!/bin/bash

#SBATCH -J S2Cpp
#SBATCH -o S2Cpp_%j.out
#SBATCH -e S2Cpp_%j.err
#SBATCH -N 1
#SBATCH --ntasks-per-node=12
#SBATCH --time=1:00:00
#SBATCH --mem=MaxMemPerNode

#module load gcc/gcc-4.8.0 
#module load fftw/fftw-3.3.3_openmpi-1.7.4_gcc-4.8.2_double_ib
. ~/.my.bashforSEED2COR


cppexe=/projects/life9360/code/Programs/back_up/SEED2CORpp_1.08/Seed2Cor
cd /lustre/janus_scratch/life9360/TEST_correction001
#cd /lustre/janus_scratch/life9360/TEST_001
export OMP_NUM_THREADS=12; $cppexe parameters.txt <<- END
Y
END
