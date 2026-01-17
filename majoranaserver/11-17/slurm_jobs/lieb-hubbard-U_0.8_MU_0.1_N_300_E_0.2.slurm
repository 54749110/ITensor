#!/usr/bin/env bash

#SBATCH --job-name=lieb-hubbard-U_0.8_MU_0.1_N_300_E_0.2
#SBATCH --partition=64c512g
#SBATCH -N 1
#SBATCH --ntasks-per-node=8
#SBATCH --output=run_outputs/outputlog/lieb-hubbard-U_0.8_MU_0.1_N_300_E_0.2.out
#SBATCH --error=run_outputs/outputlog/lieb-hubbard-U_0.8_MU_0.1_N_300_E_0.2.err

module load intel-oneapi-compilers/2021.4.0
module load intel-mkl/2020.4.304

./liebhubbard9 slurm_jobs/liebinput2_U_0.8_MU_0.1_N_300_E_0.2
