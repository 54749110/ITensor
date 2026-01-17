#!/usr/bin/env bash

#SBATCH --job-name=lieb-hubbard-11-16
#SBATCH --partition=64c512g
#SBATCH -N 4
#SBATCH --ntasks-per-node=64
#SBATCH --output=%j.out
#SBATCH --error=%j.err

module load intel-oneapi-compilers/2021.4.0
module load intel-mkl/2020.4.304

./liebhubbard9 liebinput2