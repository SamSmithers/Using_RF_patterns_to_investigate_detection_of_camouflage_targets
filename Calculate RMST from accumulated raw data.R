# This is a copy of the R script used to conduct the survival analysis and output the Restricted Mean Survival
# Time that is reported in Smithers et al. (under review). For an explanation of the statistical analysis
# conducted, please refer to the 'Analysis of Search Time Data' section of the methods in the main manuscript.
# The script plots a KM survival curve for each subject and then outputs the Restricted Mean Survival Time (RMST)
# for each stimulus condition by measuring the area under the curve.

# Script written by Dr Samuel P. Smithers, Northeastern University, 2023-2024
# Last edited January 2024

# Corresponding authors SPS (s.smithers@northeastern.edu), WHL (whlew@meei.harvard.edu) and PJB (p.bex@northeastern.edu)

# Session info:
# R version 4.3.2 (2023-10-31 ucrt)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 19045)
# 
# Matrix products: default
# 
# 
# locale:
#   [1] LC_COLLATE=English_United States.utf8  LC_CTYPE=English_United States.utf8    LC_MONETARY=English_United States.utf8 LC_NUMERIC=C                          
# [5] LC_TIME=English_United States.utf8    
# 
# time zone: America/New_York
# tzcode source: internal
# 
# attached base packages:
#   [1] splines   stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] RISCA_1.0.4       mosaicData_0.20.4 ggformula_0.12.0  dplyr_1.1.4       Matrix_1.6-1.1    ggplot2_3.4.4     lattice_0.21-9    tune_1.1.2       
# [9] reticulate_1.34.0 relsurv_2.2-9     date_1.2-42       survival_3.5-7   
# 
# loaded via a namespace (and not attached):
#   [1] RColorBrewer_1.1-3    rstudioapi_0.15.0     jsonlite_1.8.8        shape_1.4.6           magrittr_2.0.3        vctrs_0.6.5           memoise_2.0.1        
# [8] penalized_0.9-52      polynom_1.4-1         htmltools_0.5.7       forcats_1.0.0         dials_1.2.0           haven_2.5.4           deSolve_1.40         
# [15] pROC_1.18.5           caret_6.0-94          parallelly_1.36.0     htmlwidgets_1.6.4     plyr_1.8.9            plotly_4.10.3         lubridate_1.9.3      
# [22] cachem_1.0.8          gam_1.22-3            mosaicCalc_0.6.0      lifecycle_1.0.4       iterators_1.0.14      pkgconfig_2.0.3       R6_2.5.1             
# [29] fastmap_1.1.1         future_1.33.1         glmnetUtils_1.1.9     digest_0.6.33         numDeriv_2016.8-1.1   colorspace_2.1-0      furrr_0.3.1          
# [36] metR_0.14.1           mosaic_1.9.0          pkgload_1.3.3         ncvreg_3.14.1         fansi_1.0.6           yardstick_1.2.0       timechange_0.2.0     
# [43] nnls_1.5              httr_1.4.7            compiler_4.3.2        withr_2.5.2           doParallel_1.0.17     backports_1.4.1       SuperLearner_2.0-28.1
# [50] MASS_7.3-60           lava_1.7.3            ModelMetrics_1.2.2.2  tools_4.3.2           mstate_0.3.2          survivalmodels_0.1.13 future.apply_1.11.1  
# [57] nnet_7.3-19           glue_1.6.2            quadprog_1.5-8        DiagrammeR_1.0.10     nlme_3.1-163          grid_4.3.2            checkmate_2.3.1      
# [64] reshape2_1.4.4        generics_0.1.3        orthopolynom_1.0-6.1  recipes_1.0.9         gtable_0.3.4          Ryacas_1.1.5          labelled_2.12.0      
# [71] class_7.3-22          tidyr_1.3.0           data.table_1.14.10    hms_1.1.3             rsample_1.2.0         sp_2.1-2              Deriv_4.1.3          
# [78] utf8_1.2.4            foreach_1.5.2         pillar_1.9.0          stringr_1.5.1         lhs_1.1.6             tidyselect_1.2.0      muhaz_1.2.6.4        
# [85] gridExtra_2.3         hdnom_6.0.2           stats4_4.3.2          statmod_1.5.0         hardhat_1.3.0         mosaicCore_0.9.4.0    timeDate_4032.109    
# [92] visNetwork_2.1.2      stringi_1.8.3         lazyeval_0.2.2        DiceDesign_1.10       pec_2023.04.12        workflows_1.1.3       codetools_0.2-19     
# [99] kernlab_0.9-32        data.tree_1.1.0       tibble_3.2.1          cli_3.6.2             rpart_4.1.21          randomForestSRC_3.2.3 munsell_0.5.0        
# [106] Rcpp_1.0.11           globals_0.16.2        png_0.1-8             parallel_4.3.2        gower_1.0.1           parsnip_1.1.1         timeROC_0.4          
# [113] cubature_2.1.0        GPfit_1.0-8           listenv_0.9.0         glmnet_4.1-8          timereg_2.0.5         viridisLite_0.4.2     mvtnorm_1.2-4        
# [120] ipred_0.9-14          flexsurv_2.2.2        scales_1.3.0          prodlim_2023.08.28    ggridges_0.5.5        purrr_1.0.2           rlang_1.1.2          


## Required dependencies/packages. 
library(RISCA) # Note that the package RISCA includes the survival package, along with several over packages, as dependences.
#chrome-extension://efaidnbmnnnibpcajpcglclefindmkaj/https://cran.r-project.org/web/packages/RISCA/RISCA.pdf
#https://www.rdocumentation.org/packages/RISCA/versions/0.9/topics/rmst


rm(list = ls(all = TRUE))


#Set up Working directory 
setwd("***Pathway to folder containing this R script and all '_accumulated_rawdata.csv' files***")

#/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
# Select the name and data file of the experiment you want to process the data for.

ExpName <- "Exp1"
data <- read.csv("Exp1_accumulated_rawdata.csv", header = T)
 
# ExpName <- "Exp2"
# data <-read.csv("Exp2_accumulated_rawdata.csv", header=T)
 
# ExpName <- "Exp3"
# data <-read.csv("Exp3_accumulated_rawdata.csv", header=T)

# ExpName <- "Exp4"
# data <-read.csv("Exp4_accumulated_rawdata.csv", header=T)

#\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/

# Get list of all subjects

SubjectList <- unique(data["Subject_ID"])
ConditionList <- unique(data["Condition_name"])
CumulativeData <- data.frame()

for (ii in 1:nrow(SubjectList)) {
  subData <- subset(data, Subject_ID == SubjectList[ii, 1])
  Check_n_trials <- subset(subData, Block == 1)
  print(paste("Number of trials for", SubjectList[ii, 1],"in block 1 =", 
              count(unique(Check_n_trials["Order"])))) # Check that all subjects did all 26 trials and/or all 26 trials saved correctly (should be 13 trial per block)
  Check_n_trials <- subset(subData, Block == 2)
  print(paste("Number of trials for",SubjectList[ii, 1],"in block 2 =",
              count(unique(Check_n_trials["Order"]))))
  
  # Compute KM survival curve
  fit <-
    summary(survfit(Surv(TimeOfEvent, Found_or_Censored) ~ Condition_name, data =
                      subData))
  
  # For each stimulus condition extract condition info and save this along side the calculated RMST
  temp <- data.frame()
  for (cc in 1:nrow(ConditionList)) {
    tempConditionData <- subset(subData, Condition_name == ConditionList[cc, 1])
    temp[cc, "Subject_ID"] <- tempConditionData[1, "Subject_ID"]
    temp[cc, "Condition_name"] <- tempConditionData[1, "Condition_name"]
    temp[cc, "Radial_frequency"] <- tempConditionData[1, "Radial_frequency"]
    temp[cc, "Amplitude"] <- tempConditionData[1, "Amplitude"]
    temp[cc, "perimeter_cm"] <- tempConditionData[1, "perimeter_cm"]
    temp[cc, "Condition_number"] <- tempConditionData[1, "Condition_number"]
    
    # Get RMST (area under the curve) from the KM curve
    temp[cc, "RMST"] <- 
      rmst(times = fit$time[as.character(fit$strata) == paste("Condition_name=", ConditionList[cc, 1], sep = "")], 
           surv.rates = fit$surv[as.character(fit$strata) == paste("Condition_name=", ConditionList[cc, 1], sep = "")],
           max.time = 180,type = "s")
  }
  
  temp <- temp[order(temp$Condition_number), ]
  CumulativeData <- rbind(CumulativeData, temp)
  rm(temp)
}


#Save
today <- Sys.Date()
today <- format(today, "%Y%m%d")
FileName <-paste(today, "_RF_", ExpName, "_RMST data.csv", sep = "")
write.csv(CumulativeData, FileName, row.names = FALSE)
