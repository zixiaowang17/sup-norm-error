#!/bin/bash
set -e
cd "$(dirname "$0")"

# Run all gamma values in parallel (batched)
BATCH=8
count=0
for gamma in $(Rscript -e "source('config_dd.R'); cat(sprintf('%.2f\n', GAMMA_GRID))"); do
  Rscript 01_run_dd.R "$gamma" &
  count=$((count + 1))
  if [ "$count" -ge "$BATCH" ]; then
    wait
    count=0
  fi
done
wait

Rscript 02_plot_dd.R
