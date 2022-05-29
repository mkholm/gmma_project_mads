#!/bin/bash

Rscript gmma00_fetch.r
Rscript gmma01_structure.r podgornaia15.txt > gmma01_structure.out
# Rscript gmma02_fit_individual.r > gmma02_fit_individual.out
Rscript gmma02_init_zero.r > gmma02_init.out
Rscript gmma03_graph.r gmma_structured.rda > gmma03_graph.out
# Rscript gmma04_global_estimation_baselines.r -2 > gmma04_global_estimation.out
Rscript gmma04_global_estimation.r > gmma04_global_estimation.out
Rscript gmma05_analysis.r gmma_fit_global.rda > gmma05_analysis.out
