#!/bin/bash

# load output of gmma synthetic with scan var, extract results
scan_out=$1
out="scan_dat_pruned.txt" 

# header
cat $scan_out | grep "LAB" -m1 | cut -c 19- >> $out
# data and remove preceeding text on each line
cat $scan_out | grep "VAL" | cut -c 19- | rev | cut -c 2- | rev >> $out  




