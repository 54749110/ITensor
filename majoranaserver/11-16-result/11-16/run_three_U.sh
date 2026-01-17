#!/usr/bin/env bash
set -euo pipefail

# run_three_U.sh
# Usage: ./run_three_U.sh [inputfile]
# If inputfile omitted, default is 'liebinput' in current dir.

INPUT=${1:-liebinput2}
BINARY=./liebhubbard9
OUTDIR=run_outputs

# 22 U , 10 mu , 10 E values ,2 N values
#U_VALUES=( -2 -1.5 -1 -0.5 -0.3 -0.2 -0.1 -0.01 -0.001 0 0.001 0.01 0.1 0.2 0.3 0.5 0.8 1 1.2 1.5 2 3 5 )
#MU_VALUES=( 0 0.001 0.01 0.1 0.2 0.3 0.5 0.8 1.0 2.0)
E_VALUES=(0.01  0.05  0.2 )
N_VALUES=(300 303)

# justify whether input file exists 
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

# create output directory
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

# for E in "${E_VALUES[@]}"; do
#   INFILE="${INPUT}.E_${E}"
#   # create modified input: replace a line 'E = <val>' or 'E <val>' or 'E:' style; we do a safe sed replace
#   # Try to replace a line starting with 'E' (possibly with whitespace and separators)
#   awk -v EVAL="$E" '
#     BEGIN{replaced=0}
#     { if(!replaced && $0 ~ /^\s*E\s*[:= ]/) { print "E = " EVAL; replaced=1 } else print $0 }
#     END{ if(!replaced) print "E = " EVAL }
#   ' "$INPUT" > "$INFILE"

#   LOGOUT="$OUTDIR/run_E_${E}.log"
#   LOGERR="$OUTDIR/run_E_${E}.err"
#   echo "Running with E=$E -> input $INFILE (logs: $LOGOUT, $LOGERR)"
#   # Run and capture stdout/stderr
#   (time "$BINARY" "$INFILE") > "$LOGOUT" 2> "$LOGERR" || echo "Run with E=$E exited non-zero; see $LOGERR"
# done

# echo "All runs completed. Logs in $OUTDIR/"



for N in "${N_VALUES[@]}"; do
  for E in "${E_VALUES[@]}"; do
    (
      INFILE="${INPUT}.N_${N}.E_${E}"
      LOGOUT="$OUTDIR/run_N_${N}_E_${E}.log"
      LOGERR="$OUTDIR/run_N_${N}_E_${E}.err"

      # 修改输入文件中的 N 和 E 值
      awk -v NVAL="$N" -v EVAL="$E" '
      BEGIN{replaced_N=0; replaced_E=0}
      { 
        if(!replaced_N && $0 ~ /^\s*N\s*[:= ]/) { print "N = " NVAL; replaced_N=1 } 
        else if(!replaced_E && $0 ~ /^\s*E\s*[:= ]/) { print "E = " EVAL; replaced_E=1 } 
        else print $0 
      }
      END{ 
        if(!replaced_N) print "N = " NVAL 
        if(!replaced_E) print "E = " EVAL 
      }
      ' "$INPUT" > "$INFILE"

      # 运行程序
      echo "Running with N=$N, E=$E -> input $INFILE (logs: $LOGOUT, $LOGERR)"
      (time "$BINARY" "$INFILE") > "$LOGOUT" 2> "$LOGERR" || echo "Run with N=$N, E=$E exited non-zero; see $LOGERR"
    ) &
  done
done

# 等待所有后台任务完成
wait

echo "All runs completed. Logs in $OUTDIR/"