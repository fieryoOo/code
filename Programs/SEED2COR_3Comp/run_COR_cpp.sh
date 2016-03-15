#!/bin/bash

#SBATCH -J S2Cpp
#SBATCH -o S2Cpp_%j.out
#SBATCH -e S2Cpp_%j.err
#SBATCH -N 1
#SBATCH --ntasks-per-node=12
#SBATCH --time=24:00:00
#SBATCH --mem=MaxMemPerNode

module load gcc/gcc-4.9.1
module load fftw
#. ~/.my.bashforSEED2COR


cppexe=/projects/life9360/code/Programs/SEED2CORpp_1.10/Seed2Cor
cd /lustre/janus_scratch/life9360/COR_TEST_forYe_norm3Comp
#cd /lustre/janus_scratch/life9360/TEST_001
export OMP_NUM_THREADS=12; $cppexe parameters_1.10.txt <<- END
Y
END
