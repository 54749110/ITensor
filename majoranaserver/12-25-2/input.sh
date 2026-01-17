#!/usr/bin/env bash
set -euo pipefail

# 参数值
#U_VALUES=(-2 -1.5 -1 -0.5 -0.3 -0.2 -0.1 -0.01 -0.001 0 0.001 0.01 0.1 0.2 0.3 0.5 0.8 1 1.2 1.5 2 3 5)
#MU_VALUES=(0 0.001 0.01 0.1 0.3 0.5 1.0 2.0)
#E_VALUES=(0.005 0.01 0.03 0.05 0.1 0.2)
N_VALUES=(150 90)
U_VALUES=(0.7 0.71 0.72 0.73 0.74 0.75 0.76 0.77 0.78 0.79 0.8 0.81 0.82 0.83 0.84 0.85 0.86 0.87 0.88 0.89 0.9 0.91 0.92 0.93 0.94 0.95 0.96 0.97 0.98 0.99 1 1.01 1.02 1.03 1.04 1.05 1.06 1.07 1.08 1.09 1.1 1.12 1.14 1.16 1.18 1.2 1.22 1.24 1.26 1.28 1.3 1.35 1.4 1.45 1.5 1.6 1.7 1.8 2 3 4 5 10 )
MU_VALUES=(0)
E_VALUES=(0)

# 检查输入文件
INPUT="liebinput2"
if [ ! -f "$INPUT" ]; then
  echo "Input file $INPUT not found." >&2
  exit 1
fi


OUTDIR2=run_outputs

# 创建输出目录
mkdir -p "$OUTDIR2"
mkdir -p "$OUTDIR2/outputpsi"
mkdir -p "$OUTDIR2/outputcheck"
mkdir -p "$OUTDIR2/outputgamma3gammaj"
mkdir -p "$OUTDIR2/outputgammajgammaN"
mkdir -p "$OUTDIR2/outputgammaigammaj"
mkdir -p "$OUTDIR2/outputcicj"
mkdir -p "$OUTDIR2/outputcicdagj"
mkdir -p "$OUTDIR2/outputni"
mkdir -p "$OUTDIR2/outputninj"
mkdir -p "$OUTDIR2/outputEE"
mkdir -p "$OUTDIR2/outputdEdU"
mkdir -p "$OUTDIR2/outputdensitycorr"
mkdir -p "$OUTDIR2/outputCDWorder"
mkdir -p "$OUTDIR2/outputlog"
mkdir -p "$OUTDIR2/outputGSenergy"

# 输出目录
OUTDIR="slurm_jobs"
mkdir -p "$OUTDIR"

# 遍历所有参数组合
for U in "${U_VALUES[@]}"; do
  for MU in "${MU_VALUES[@]}"; do
    for E in "${E_VALUES[@]}"; do
      for N in "${N_VALUES[@]}"; do
        # 修改 liebinput2 文件中的参数值
        TASK_INPUT="$OUTDIR/liebinput2_U_${U}_MU_${MU}_N_${N}_E_${E}"

        awk -v UVAL="$U" -v MUVAL="$MU" -v NVAL="$N" -v EVAL="$E" '
        BEGIN{replaced_U=0; replaced_MU=0; replaced_N=0; replaced_E=0}
        { 
          if(!replaced_U && $0 ~ /^\s*U\s*[:= ]/) { print "U = " UVAL; replaced_U=1 } 
          else if(!replaced_MU && $0 ~ /^\s*mu\s*[:= ]/) { print "mu = " MUVAL; replaced_MU=1 } 
          else if(!replaced_N && $0 ~ /^\s*N\s*[:= ]/) { print "N = " NVAL; replaced_N=1 } 
          else if(!replaced_E && $0 ~ /^\s*E\s*[:= ]/) { print "E = " EVAL; replaced_E=1 } 
          else print $0 
        }
        END{ 
          if(!replaced_U) print "U = " UVAL 
          if(!replaced_MU) print "mu = " MUVAL 
          if(!replaced_N) print "N = " NVAL 
          if(!replaced_E) print "E = " EVAL 
        }
        ' "$INPUT" > "$TASK_INPUT"

        # 生成对应的 B.slurm 文件
        JOB_NAME="lieb-hubbard-U_${U}_MU_${MU}_N_${N}_E_${E}"
        SLURM_FILE="$OUTDIR/${JOB_NAME}.slurm"
        cat > "$SLURM_FILE" <<EOF
#!/usr/bin/env bash

#SBATCH --job-name=$JOB_NAME
#SBATCH --partition=64c512g
#SBATCH --nodes 1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --output=$OUTDIR2/outputlog/${JOB_NAME}.out
#SBATCH --error=$OUTDIR2/outputlog/${JOB_NAME}.err

module load intel-oneapi-compilers/2021.4.0
module load intel-mkl/2020.4.304

./liebhubbard9 $TASK_INPUT
EOF

        # 提交任务
        sbatch "$SLURM_FILE"
      done
    done
  done
done
