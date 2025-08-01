# Supplementary materials: Using radial frequency patterns to investigate the effect of shape features on the detection of camouflage targets 
### Shared data for: Samuel P. Smithers, Wei Hau Lew, Yulong Shao, Kevin Travieso, Daniel R. Coates, Peter J. Bex. (under review) *Using radial frequency patterns to investigate the effect of shape features on the detection of camouflage targets.*

## Repository contents
- Accumulated raw date from all experiments set up for conducting survival analysis within R (files ending "*_accumulated_rawdata.csv*").
- ```Calculate RMST from accumulated raw data.R```: R code used to conduct the survival analysis and output the Restricted Mean Survival Time (RMST) as "*_RMST_data.csv*".
- Restricted Mean Survival Time (RMST) data from all experiments generated by ```Calculate RMST from accumulated raw data.R``` (files ending "*_RMST_data.csv*").
- ```Plots and statistical analysis for RMST data.R```: R code used to generate the figures for, and perform the statistical analysis on, the RMST data reported in the main manuscript.
- Accumulated mouse click data from all experiment (files ending "*_accumulated mouse click data.csv*").
-  ```Graphs and stats for false positives.R```: R code used to generate the figures for, and perform the statistical analysis of the number of, false positives reported in the main manuscript.
-  ```Effect of target position on detection time.R```: R code used to generate supplimentary figures S2-S5.
-  ```Mouse click position relative to target pos.R```: R code used to generate supplimentary figures S6-S9.

## Running the R code 
### Requirements & dependencies
- R version >= 4
- The following packages must be installed: 
  - RISCA (Note that the package RISCA includes the survival package, and several other packages, as dependences) 
  - ggplot2
  - viridis
  - scales
  - grid
  - png
  - gridExtra (optional to arrange/combine plots for saving)
  - svglite (optional to save plots)
  - Matrix
  - lme4
  - performance
  - dplyr
  - gridExtra
  - glmmTMB
  - DHARMa
  - bannerCommenter (optional)
 
#### Running instructions
1. Manually download all files from this repository. Alternatively, you can clone this repository ([instructions for Rstudio here](https://datacarpentry.org/rr-version-control/03-git-in-rstudio/index.html)). 
2. Install any uninstalled dependencies.
3. Open R scripts and follow instructions within the script. 
