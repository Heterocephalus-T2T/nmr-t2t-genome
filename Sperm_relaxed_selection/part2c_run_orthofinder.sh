#!/bin/bash
# Part 2c: OrthoFinder on the cleaned proteomes.
# Flags are frozen; downstream steps pin the Results_* directory this creates.

nohup orthofinder -f results/02b_cleaned_proteins/ -o results/02c_orthofinder \
    -t 44 -a 44 -M msa -S diamond_ultra_sens > orthofinder.log 2>&1 &
echo "PID $! - tail -f orthofinder.log to watch"
