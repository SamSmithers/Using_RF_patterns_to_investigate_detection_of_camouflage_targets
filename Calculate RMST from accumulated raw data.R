# This is a copy of the R script used to conduct the survival analysis and output the Restricted Mean Survival
# Time that is reported in Smithers et al. (under review). For an explanation of the statistical analysis
# conducted, please refer to the 'Analysis of Search Time Data' section of the methods in the main manuscript.
# The script plots a KM survival curve for each subject and then outputs the Restricted Mean Survival Time (RMST)
# for each stimulus condition by measuring the area under the curve.

# Script written by Dr Samuel P. Smithers, Northeastern University, 2023-2025
# Last edited April 2025

# Corresponding authors SPS (s.smithers@northeastern.edu), WHL (whlew@meei.harvard.edu) and PJB (p.bex@northeastern.edu)

# R version 4.2.2 (2022-10-31 ucrt)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 19045)
# 
# Matrix products: default
# 
# locale:
#   [1] LC_COLLATE=English_United States.utf8  LC_CTYPE=English_United States.utf8    LC_MONETARY=English_United States.utf8 LC_NUMERIC=C                          
# [5] LC_TIME=English_United States.utf8    
# 
# attached base packages:
#   [1] splines   stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] RISCA_1.0.4       mosaicData_0.20.3 ggformula_0.10.2  dplyr_1.1.1       Matrix_1.5-1      ggplot2_3.4.3     lattice_0.20-45   tune_1.1.0        reticulate_1.28  
# [10] relsurv_2.2-9     date_1.2-42       survival_3.4-0   
# 
# loaded via a namespace (and not attached):
# [1] backports_1.4.1       workflows_1.1.3       plyr_1.8.8            mstate_0.3.2          lazyeval_0.2.2        sp_1.6-0              hdnom_6.0.1           listenv_0.9.0        
# [9] digest_0.6.30         foreach_1.5.2         htmltools_0.5.3       yardstick_1.1.0       parsnip_1.0.4         fansi_1.0.3           magrittr_2.0.3        checkmate_2.1.0      
# [17] memoise_2.0.1         doParallel_1.0.17     mosaicCore_0.9.2.1    recipes_1.0.5         globals_0.16.2        flexsurv_2.2.2        gower_1.0.1           hardhat_1.3.0        
# [25] timechange_0.2.0      rsample_1.1.1         dials_1.2.0           colorspace_2.0-3      haven_2.5.2           jsonlite_1.8.3        iterators_1.0.14      glue_1.6.2           
# [33] polyclip_1.10-4       gtable_0.3.1          nnls_1.4              ipred_0.9-14          kernlab_0.9-32        future.apply_1.10.0   shape_1.4.6           scales_1.2.1         
# [41] penalized_0.9-52      mvtnorm_1.1-3         data.tree_1.0.0       mosaic_1.8.4.2        Rcpp_1.0.9            viridisLite_0.4.1     ggstance_0.3.6        GPfit_1.0-8          
# [49] Ryacas_1.1.5          deSolve_1.35          randomForestSRC_3.2.1 stats4_4.2.2          lava_1.7.2.1          prodlim_2019.11.13    glmnet_4.1-7          SuperLearner_2.0-28  
# [57] htmlwidgets_1.5.4     httr_1.4.4            DiagrammeR_1.0.9      RColorBrewer_1.1-3    ellipsis_0.3.2        pkgconfig_2.0.3       farver_2.1.1          nnet_7.3-18          
# [65] utf8_1.2.2            caret_6.0-94          mosaicCalc_0.6.0      tidyselect_1.2.0      rlang_1.1.0           DiceDesign_1.9        reshape2_1.4.4        polynom_1.4-1        
# [73] munsell_0.5.0         tools_4.2.2           visNetwork_2.1.2      cachem_1.0.6          cli_3.6.1             generics_0.1.3        glmnetUtils_1.1.8     ggridges_0.5.4       
# [81] stringr_1.4.1         fastmap_1.1.0         ncvreg_3.13.0         ModelMetrics_1.2.2.2  timereg_2.0.5         purrr_1.0.1           pec_2022.05.04        future_1.32.0        
# [89] nlme_3.1-160          gam_1.22-2            survivalmodels_0.1.13 compiler_4.2.2        rstudioapi_0.14       plotly_4.10.1         png_0.1-7             tibble_3.2.1         
# [97] lhs_1.1.6             statmod_1.5.0         tweenr_2.0.2          stringi_1.7.8         forcats_0.5.2         cubature_2.0.4.6      vctrs_0.6.1           pillar_1.8.1         
# [105] lifecycle_1.0.3       furrr_0.3.1           data.table_1.14.8     orthopolynom_1.0-6.1  R6_2.5.1              muhaz_1.2.6.4         gridExtra_2.3         timeROC_0.4          
# [113] parallelly_1.35.0     codetools_0.2-18      MASS_7.3-58.1         pkgload_1.3.1         withr_2.5.0           Deriv_4.1.3           parallel_4.2.2        hms_1.1.2            
# [121] quadprog_1.5-8        grid_4.2.2            rpart_4.1.19          labelled_2.10.0       metR_0.14.0           timeDate_4022.108     tidyr_1.2.1           class_7.3-20         
# [129] ggforce_0.4.1         pROC_1.18.0           numDeriv_2016.8-1.1   lubridate_1.9.2 

## Required dependencies/packages. 
library(RISCA) # Note that the package RISCA includes the survival package, along with several over packages, as dependences.
#chrome-extension://efaidnbmnnnibpcajpcglclefindmkaj/https://cran.r-project.org/web/packages/RISCA/RISCA.pdf
#https://www.rdocumentation.org/packages/RISCA/versions/0.9/topics/rmst


rm(list = ls(all = TRUE))


#Set up Working directory 
setwd("***Pathway to folder containing this R script and all '_accumulated_rawdata.csv' files***")

#/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
# Select the name and data file of the experiment you want to process the data for.

ExpName <- "Exp1A"
data <- read.csv("Exp1A_accumulated_rawdata.csv", header = T)

#ExpName <- "Exp1B"
#data <-read.csv("Exp1B_accumulated_rawdata.csv", header=T)

# ExpName <- "Exp2A"
# data <-read.csv("Exp2A_accumulated_rawdata.csv", header=T)

# ExpName <- "Exp2B"
# data <-read.csv("Exp2B_accumulated_rawdata.csv", header=T)

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
FileName <-paste(today, '_', ExpName, "_RMST_data.csv", sep = "")
write.csv(CumulativeData, FileName, row.names = FALSE)


##================================================================
##        Plot individual survival curves for each subject.      =
##================================================================
library(survminer)
library(ggplot2)
CustomPalette = c("#DC267F", "#FE6100", "#FFB000") # Colour blind friendly. Matches colour scheme of main results figures

if (ExpName == "Exp1A" || ExpName == "Exp2A" ){ #Area matched experiments 
  data$Condition_name <- factor(data$Condition_name, levels =c(
    'RF0A0','RF3A2','RF3A5','RF3A10','RF4A2','RF4A5','RF4A10','RF5A2','RF5A5','RF5A10','RF8A2','RF8A5','RF8A10')) 
} else if (ExpName == "Exp1B" || ExpName == "Exp2B" ){ #Size matched experiments
  data$Condition_name <- factor(data$Condition_name, levels =c(
    'RF0A0','RF8A2','RF8A5','RF8A10','RF10A2','RF10A5','RF10A10','RF14A2','RF14A5','RF14A10','RF20A2','RF20A5','RF20A10')) 
}

SubjectList <- unique(data["Subject_ID"])  


for (aa in 1:nrow(SubjectList)){
  PlotData <- subset(data, Subject_ID == SubjectList[aa,1])
  PlotData <- subset(PlotData, Radial_frequency==8) # You can use this row to change the RF being plotted or comment it out to plot all RFs on the same figure. 
  fit<- survfit(Surv(TimeOfEvent,Found_or_Censored) ~ Condition_name, data=PlotData) #Fit KM model for plotting
  Surv_all <- ggsurvplot(fit, data = PlotData, title = paste("KM plot for",SubjectList[aa,1]), legend ="right", palette = CustomPalette) 
  print(Surv_all)
}

# Example from supplementary material. 
#Exp 1B s11 RF8
PlotData <- subset(data, Subject_ID == SubjectList[10,1])
PlotData <- subset(PlotData, Radial_frequency==8)
fit<- survfit(Surv(TimeOfEvent,Found_or_Censored) ~ Condition_name, data=PlotData) #Fit KM model for plotting
Surv_all <- ggsurvplot(fit, data = PlotData, title = paste("KM plot for",SubjectList[aa,1]), legend ="right", palette = CustomPalette) 
print(Surv_all)

