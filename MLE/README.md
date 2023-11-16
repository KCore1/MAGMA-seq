## Maximum Likelihood Estimation ##

Updated: October 10, 2022

Sequencing data must be formated as in `test.csv`. For generation of fake sorting data for testing purposes, use the `new_generate_data_sort.ipynb` notebook. 
Use the `calc_sigma_from_data.ipynb` notebook to get sigma estimates from labeled populations.  
Submission scripts for alpine and blanca: `run_mle_alpine.sh` and `run_mle_blanca.sh`

Parameters defined:

CSV file parameters:
i: variant
j: bin
Label/k: label concentration
Count/rijk: number of reads
Hj: max fluorescence limit of bin
Gj: min fluorescence limit of bin
nk: total number of cells sorted at label concentration
fi: variant frequency in reference population
njk: cells sorted into given bin
A: variant-dependent maximum fluorescence

Global parameters:
B: variant-independent cell autofluorescence
sigma: standard deviation of population
Kd_pred: Kd guess value
bounds: solver limits

For brightness adjustment:
	1. Perform haplotyping and scanning as normal
	2. Add column with name "brightness_adjust" with values to be multiplied by Fmax (can use the adjust_brightness.py script to automate this)
	3. Run MLE with brightness_adjust=True 
