# Supplementary materials: Using radial frequency patterns to investigate the effect of shape features on the detection of camouflage targets 
### Shared data for: Samuel P. Smithers, Wei Hau Lew, Yulong Shao, Kevin Travieso, Daniel R. Coates, Peter J. Bex. (under review) *Using radial frequency patterns to investigate the effect of shape features on the detection of camouflage targets.*

## Repository contents
- Accumulated raw date from Experiments 1-4 set up for conducting survival analysis within R (files ending "*_accumulated_rawdata.csv*").
- ```Calculate RMST from accumulated raw data.R```: R code used to conduct the survival analysis and output the Restricted Mean Survival Time (RMST) as "*_RMST_data.csv*".
- Restricted Mean Survival Time (RMST) data from Experiments 1-4 generated by ```Calculate RMST from accumulated raw data.R``` (files ending "*_RMST_data.csv*").
- ```Plots and statistical analysis for RMST data.R```: R code used to generate the figures for, and perform the statistical analysis on, the RMST data from Experiments 1-4 reported in the main manuscript.

## Running the R code 
### Requirements & dependencies
- R version >= 4
- The following packages must be installed: 
  - RISCA (Note that the package RISCA includes the survival package, and several other packages, as dependences) 
  - ggplot2
  - scales
  - gridExtra (optional to arrange/combine plots for saving)
  - svglite (optional to save plots)
  - Matrix
  - lme4
  - performance
  - bannerCommenter (optional)

#### Running instructions
1. Manually download ```Calculate RMST from accumulated raw data.R``` and ```Plots and statistical analysis for RMST data.R``` from this repository. You will also need to download the "*_accumulated_rawdata.csv*" and "*_RMST_data.csv*" data files. Alternatively, you can clone this repository ([instructions for Rstudio here](https://datacarpentry.org/rr-version-control/03-git-in-rstudio/index.html)). 
2. Install any uninstalled dependencies, above.
3. Open ```Calculate RMST from accumulated raw data.R``` and/or ```Plots and statistical analysis for RMST data.R``` and follow instructions within the script. 
