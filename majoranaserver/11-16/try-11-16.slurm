#!/usr/bin/env bash

#SBATCH --job-name=lieb-hubbard-11-16
#SBATCH --partition=64c512g
#SBATCH -N 4
#SBATCH --ntasks-per-node=64
#SBATCH --output=%j.out
#SBATCH --error=%j.err

module load intel-oneapi-compilers/2021.4.0
module load intel-mkl/2020.4.304

set -euo pipefail

# run_three_U.sh
# Usage: ./run_three_U.sh [inputfile]
# If inputfile omitted, default is 'liebinput2' in current dir.

INPUT=${1:-liebinput2}
BINARY=./liebhubbard9
OUTDIR=run_outputs

# 参数值 22*8*6*2 =2112 个组合
U_VALUES=(-2 -1.5 -1 -0.5 -0.3 -0.2 -0.1 -0.01 -0.001 0 0.001 0.01 0.1 0.2 0.3 0.5 0.8 1 1.2 1.5 2 3 5)
MU_VALUES=(0 0.001 0.01 0.1 0.3 0.5 1.0 2.0)
E_VALUES=(0.005 0.01 0.03 0.05 0.1 0.2)
N_VALUES=(60 63)

# 检查输入文件和可执行文件
if [ ! -f "$INPUT" ]; then
  echo "Input file $INPUT not found." >&2
  exit 2
fi
if [ ! -x "$BINARY" ]; then
  if [ -f "$BINARY" ]; then
    echo "$BINARY exists but is not executable. Please build and/or chmod +x $BINARY" >&2
    exit 3
  else
    echo "$BINARY not found. Please compile liebhubbard9 first." >&2
    exit 4
  fi
fi

# 创建输出目录
mkdir -p "$OUTDIR"
mkdir -p "outputpsi"
mkdir -p "outputcheck"
mkdir -p "outputgamma3gammaj"
mkdir -p "outputgammajgammaN"
mkdir -p "outputgammaigammaj"
mkdir -p "outputcicj"
mkdir -p "outputcicdagj"
mkdir -p "outputni"
mkdir -p "outputninj"
mkdir -p "outputEE"
mkdir -p "outputdEdU"
mkdir -p "outputdensitycorr"
mkdir -p "outputCDWorder"

# 定义运行函数
run_simulation() {
  local U=$1
  local MU=$2
  local E=$3
  local N=$4

  LOGOUT="$OUTDIR/run_U_${U}_MU_${MU}_N_${N}_E_${E}.log"
  LOGERR="$OUTDIR/run_U_${U}_MU_${MU}_N_${N}_E_${E}.err"

  # 修改输入文件中的 U、MU、N 和 E 值
  awk -v UVAL="$U" -v MUVAL="$MU" -v NVAL="$N" -v EVAL="$E" '
  BEGIN{replaced_U=0; replaced_MU=0; replaced_N=0; replaced_E=0}
  { 
    if(!replaced_U && $0 ~ /^\s*U\s*[:= ]/) { print "U = " UVAL; replaced_U=1 } 
    else if(!replaced_MU && $0 ~ /^\s*MU\s*[:= ]/) { print "MU = " MUVAL; replaced_MU=1 } 
    else if(!replaced_N && $0 ~ /^\s*N\s*[:= ]/) { print "N = " NVAL; replaced_N=1 } 
    else if(!replaced_E && $0 ~ /^\s*E\s*[:= ]/) { print "E = " EVAL; replaced_E=1 } 
    else print $0 
  }
  END{ 
    if(!replaced_U) print "U = " UVAL 
    if(!replaced_MU) print "MU = " MUVAL 
    if(!replaced_N) print "N = " NVAL 
    if(!replaced_E) print "E = " EVAL 
  }
  ' "$INPUT" > "${INPUT}.tmp" && mv "${INPUT}.tmp" "$INPUT"

  # 运行程序
  echo "Running with U=$U, MU=$MU, N=$N, E=$E -> input $INPUT (logs: $LOGOUT, $LOGERR)"
  (time "$BINARY" "$INPUT") > "$LOGOUT" 2> "$LOGERR" || echo "Run with U=$U, MU=$MU, N=$N, E=$E exited non-zero; see $LOGERR"
}

export -f run_simulation

# 使用 GNU Parallel 并行运行
parallel -j "$(nproc)" run_simulation ::: "${U_VALUES[@]}" ::: "${MU_VALUES[@]}" ::: "${E_VALUES[@]}" ::: "${N_VALUES[@]}"

echo "All runs completed. Logs in $OUTDIR/"