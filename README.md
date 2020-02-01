# FPCA_GWAS
code for the sorghum time-series GWAS project

## 1.nonparametric_regression
`nonparametric_fitting.R`: The R script for nonparametric fitting on the sorghum growth curves.

`example_data.csv`: An example dataset of plant height on different date points. The first column is the sorghum line name. The second column represents time points. The last column is the real plant heights. 

`fitting_results.png`: An image generated from the script showing the fitted growth curve. 


## 2.FPCA

R script to conduct functional principal component analysis.  
`fpca.R`: The main FPCA function. Users only running the main function will be enough. 

`pca_fun.R`, `pca_score.R`, `tuning_nointer.R` : Other three functions used in the main function. 

`data.Rdata`: An example of the input dataset. It includes two lists. The first is a list of time points (x-axis) for each sample. The second is a list of observations/measurements (y-axis) corresponding to the time points for each sample.