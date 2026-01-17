#!/usr/bin/env bash

#SBATCH --job-name=lieb-hubbard-U_1.20_MU_1.0_N_300_E_0
#SBATCH --partition=64c512g
#SBATCH -N 1
#SBATCH --ntasks-per-node=16
#SBATCH --output=run_outputs/outputlog/lieb-hubbard-U_1.20_MU_1.0_N_300_E_0.out
#SBATCH --error=run_outputs/outputlog/lieb-hubbard-U_1.20_MU_1.0_N_300_E_0.err

module load intel-oneapi-compilers/2021.4.0
module load intel-mkl/2020.4.304

./liebhubbard9 slurm_jobs/liebinput2_U_1.20_MU_1.0_N_300_E_0
