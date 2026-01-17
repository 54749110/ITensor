#!/usr/bin/env bash

#SBATCH --job-name=lieb-hubbard-U_0.3_MU_2_N_303_E_0.001
#SBATCH --partition=64c512g
#SBATCH -N 1
#SBATCH --ntasks-per-node=8
#SBATCH --output=run_outputs2/outputlog/lieb-hubbard-U_0.3_MU_2_N_303_E_0.001.out
#SBATCH --error=run_outputs2/outputlog/lieb-hubbard-U_0.3_MU_2_N_303_E_0.001.err

module load intel-oneapi-compilers/2021.4.0
module load intel-mkl/2020.4.304

./liebhubbard9 slurm_jobs2/liebinput2_U_0.3_MU_2_N_303_E_0.001
