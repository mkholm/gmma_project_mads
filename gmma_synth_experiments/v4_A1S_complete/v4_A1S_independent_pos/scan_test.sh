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

# remove possible output plots from possible previous rurm -r eval_synth
rm -r global_fit
rm -r lib_compo
rm -r signal_distributions
rm -r stability_distributions

# c(rep(seq(-4,-1,0.5), each = 5), rep(seq(-0.9,1,0.1), each = 5), rep(seq(1.5,8,0.5), each = 5))
for scan_var in -4.0 -4.0 -4.0 -4.0 -4.0 -3.5 -3.5 -3.5 -3.5 -3.5 -3.0 -3.0 -3.0 -3.0 -3.0 -2.5 -2.5 -2.5 -2.5 -2.5 -2.0 -2.0 -2.0 -2.0 -2.0 -1.5 -1.5 -1.5 -1.5 -1.5 -1.0 -1.0 -1.0 -1.0 -1.0 -0.9 -0.9 -0.9 -0.9 -0.9 -0.8 -0.8 -0.8 -0.8 -0.8 -0.7 -0.7 -0.7 -0.7 -0.7 -0.6 -0.6 -0.6 -0.6 -0.6 -0.5 -0.5 -0.5 -0.5 -0.5 -0.4 -0.4 -0.4 -0.4 -0.4 -0.3 -0.3 -0.3 -0.3 -0.3 -0.2 -0.2 -0.2 -0.2 -0.2 -0.1 -0.1 -0.1 -0.1 -0.1  0.0  0.0  0.0  0.0  0.0  0.1  0.1  0.1  0.1  0.1  0.2  0.2  0.2  0.2  0.2  0.3  0.3  0.3  0.3  0.3  0.4  0.4  0.4  0.4  0.4  0.5  0.5  0.5  0.5  0.5  0.6  0.6  0.6  0.6  0.6  0.7  0.7  0.7 0.7  0.7  0.8  0.8  0.8  0.8  0.8  0.9  0.9  0.9  0.9  0.9  1.0  1.0  1.0  1.0  1.0  1.5  1.5  1.5  1.5  1.5  2.0  2.0  2.0  2.0  2.0  2.5  2.5  2.5  2.5  2.5  3.0  3.0  3.0  3.0  3.0  3.5  3.5  3.5  3.5  3.5  4.0  4.0  4.0  4.0  4.0  4.5  4.5  4.5  4.5  4.5  5.0  5.0  5.0  5.0  5.0  5.5  5.5 5.5  5.5  5.5  6.0  6.0  6.0  6.0  6.0  6.5  6.5  6.5  6.5  6.5  7.0  7.0  7.0  7.0  7.0  7.5  7.5  7.5  7.5  7.5  8.0  8.0  8.0  8.0  8.0; do
    f=$(printf "Scanning variable %d" $scan_var)

    # Make synthetic data
    echo "make_synthetic_binary.r $scan_var >> $out_file"
    [[ $? == 0 ]] || { echo "An error occured"; exit 2; } # terminate script if it fails at this step

    # Run GMMA
    echo "gmma01_structure.r gmma_synthetic.txt > gmma01_structure.out"
    [[ $? == 0 ]] || { echo "An error occured"; exit 2; }

    echo "gmma02_init_zero.r > gmma02_init.out"
    [[ $? == 0 ]] || { echo "An error occured"; exit 2; }

    echo "gmma03_graph.r gmma_structured.rda > gmma03_graph.out"
    [[ $? == 0 ]] || { echo "An error occured"; exit 2; }

    echo "gmma04_global_estimation.r > gmma04_global_estimation.out"
    	# [[ $? == 0 ]] || { echo "An error occured"; exit 2; }    
	# just continue with next loop iteration if failed gmma05 error analysis
	[[ $? == 0 ]] || { echo "An error occured in gmma04"; continue; }
    
    echo "gmma05_analysis.r gmma_fit_global.rda > gmma05_analysis.out"
	# [[ $? == 0 ]] || { echo "An error occured"; exit 2; }
	# just continue with next loop iteration if failed gmma05 error analysis
	[[ $? == 0 ]] || { echo "An error occured in gmma05"; continue; }
    
    # just overwrite normal output
    echo "eval_synth.r >> $out_file"
    
done

# process output and make overview plot
# processed="scan_dat_pruned.txt" 
# make sure no results from previous run exist
#rm $processed

# header
#cat $out_file | grep "LAB" -m1 | cut -c 19- >> $processed
# data and remove preceeding text on each line
#cat $out_file | grep "VAL" | cut -c 19- | rev | cut -c 2- | rev >> $processed  

# plot processed output
#Rscript plot_synth.r $processed

