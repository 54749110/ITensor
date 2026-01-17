#!/usr/bin/env bash

#SBATCH --job-name=lieb-hubbard-U_1_MU_0.1_N_303_E_0.1
#SBATCH --partition=64c512g
#SBATCH -N 1
#SBATCH --ntasks-per-node=8
#SBATCH --output=run_outputs/outputlog/lieb-hubbard-U_1_MU_0.1_N_303_E_0.1.out
#SBATCH --error=run_outputs/outputlog/lieb-hubbard-U_1_MU_0.1_N_303_E_0.1.err

module load intel-oneapi-compilers/2021.4.0
module load intel-mkl/2020.4.304

./liebhubbard9 slurm_jobs/liebinput2_U_1_MU_0.1_N_303_E_0.1
