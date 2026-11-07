#!/bin/bash
set -e
cd "$(dirname "$0")"

Rscript 01_run_sim.R 0.5 &
Rscript 01_run_sim.R 2.0 &
wait

Rscript 02_plot.R
