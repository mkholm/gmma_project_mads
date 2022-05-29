#!/bin/bash


out_file="scan_out.txt"
[[ -f $out_file ]] && {
    while true; do
	read -p "Do you wish to overwrite output of a previous run in this directory?\n" yn
	case $yn in
            [Yy]* ) echo Typed yes, overwrite and continue; rm $out_file; break;;
            [Nn]* ) echo Typed no, exit run; exit;;
            * ) echo "Please answer yes or no.";;
	esac
    done
}


# remove output plots from possible previous run
# directories with all plots
rm -r eval_synth
rm -r global_fit
rm -r lib_compo
rm -r signal_distributions
rm -r stability_distributions


for scan_var in  0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57; do
    f=$(printf "Scanning variable %d" $scan_var)

    # Make synthetic data
    Rscript make_synthetic_binary.r $scan_var >> $out_file
    [[ $? == 0 ]] || { echo "An error occured"; exit 2; } # terminate script if it fails at this step

    # Run GMMA
    Rscript gmma01_structure.r gmma_synthetic.txt > gmma01_structure.out
    [[ $? == 0 ]] || { echo "An error occured"; exit 2; }

    Rscript gmma02_init_zero.r > gmma02_init.out
    [[ $? == 0 ]] || { echo "An error occured"; exit 2; }

    Rscript gmma03_graph.r gmma_structured.rda > gmma03_graph.out
    [[ $? == 0 ]] || { echo "An error occured"; exit 2; }

    Rscript gmma04_global_estimation.r > gmma04_global_estimation.out
    [[ $? == 0 ]] || { echo "An error occured"; exit 2; }    
    
    Rscript gmma05_analysis.r gmma_fit_global.rda > gmma05_analysis.out
    [[ $? == 0 ]] || { echo "An error occured"; exit 2; }
    
    # just overwrite normal output
    Rscript eval_synth.r >> $out_file
    
done

# process output and make overview plot
processed="scan_dat_pruned.txt" 
# make sure no results from previous run exist
rm $processed

# header
cat $out_file | grep "LAB" -m1 | cut -c 19- >> $processed
# data and remove preceeding text on each line
cat $out_file | grep "VAL" | cut -c 19- | rev | cut -c 2- | rev >> $processed  

# plot processed output
Rscript plot_synth.r $processed


# remove plots from final scan var in top level directory (clean directory)
rm eval_synth.png
rm global_fit.png
rm lib_compo.png
rm signal_distributions.png
rm stability_distributions.png

