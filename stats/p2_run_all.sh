#!/bin/bash
# Paper 2 (public-data version): main Part A analysis plus the pre-specified sensitivity analyses.
# Usage: bash stats/p2_run_all.sh      (run from the project root after p2_metrics + visual QC)
set -e
P2_PROJECT_ROOT="${P2_PROJECT_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
cd "$P2_PROJECT_ROOT"
mkdir -p logs
Rscript stats/p2_partA.R results/p2 > logs/p2_partA_main.log 2>&1 &
run() {  # name, then env assignments
  name=$1; shift
  env P2_BOOT=0 "$@" Rscript stats/p2_partA.R results/p2/sens_$name > logs/p2_partA_sens_$name.log 2>&1 &
}
run kernel_consistent P2_FILTER=kernel_consistent
run three_rounds      P2_FILTER=three_rounds
run gap90             P2_FILTER=gap90
run autochthon        P2_SEG=aut
run resected_only     P2_SETS=p2
run unresected_only   P2_SETS=p2b
wait
echo "Part A main + 6 sensitivity runs done $(date)"
