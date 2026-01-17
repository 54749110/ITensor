#!/usr/bin/env bash

#SBATCH --job-name=lieb-hubbard-U_0.75_MU_0_N_90_E_0
#SBATCH --partition=64c512g
#SBATCH --nodes 1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --output=run_outputs/outputlog/lieb-hubbard-U_0.75_MU_0_N_90_E_0.out
#SBATCH --error=run_outputs/outputlog/lieb-hubbard-U_0.75_MU_0_N_90_E_0.err

module load intel-oneapi-compilers/2021.4.0
module load intel-mkl/2020.4.304

./liebhubbard9 slurm_jobs/liebinput2_U_0.75_MU_0_N_90_E_0
