#!/usr/bin/env bash
# Runs 03_preprocess.R over all subjects in two parallel workers.
cd "$(dirname "$0")/.."
mkdir -p data/logs
S=$(Rscript -e 'cat(read.csv("results/audit_subjects.csv")$subject)')
A=$(echo $S | tr ' ' '\n' | awk 'NR%2==1' | tr '\n' ' ')
B=$(echo $S | tr ' ' '\n' | awk 'NR%2==0' | tr '\n' ' ')
Rscript code/03_preprocess.R $A > data/logs/w1.log 2>&1 &
Rscript code/03_preprocess.R $B > data/logs/w2.log 2>&1 &
wait
