# This is a copy of the R script used to generate the figures for, and perform the statistical analysis on, 
# the data for the number of target misses (false positives) from the experiments reported in Smithers et al. 
# (under review). For an explanation of the statistical analysis conducted, please refer to the 'Analysis of 
# Target Misses (false positives)' section of the methods in the main manuscript.

# Script written by Dr Samuel P. Smithers, Northeastern University, 2023-2025
# Last edited July 2025

# Corresponding author: PJB (p.bex@northeastern.edu)  

# R version 4.3.2 (2023-10-31 ucrt)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 19045)
# 
# Matrix products: default
# 
# 
# locale:
#   [1] LC_COLLATE=English_United States.utf8  LC_CTYPE=English_United States.utf8   
# [3] LC_MONETARY=English_United States.utf8 LC_NUMERIC=C                          
# [5] LC_TIME=English_United States.utf8    
# 
# time zone: America/New_York
# tzcode source: internal
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] DHARMa_0.4.6          glmmTMB_1.1.9         performance_0.10.8    lme4_1.1-35.1        
# [5] Matrix_1.6-5          gridExtra_2.3         bannerCommenter_1.0.0 dplyr_1.1.4          
# [9] ggplot2_3.4.4        
# 
# loaded via a namespace (and not attached):
#   [1] gtable_0.3.4        compiler_4.3.2      tidyselect_1.2.0    Rcpp_1.0.11         splines_4.3.2      
# [6] scales_1.3.0        boot_1.3-28.1       lattice_0.21-9      R6_2.5.1            generics_0.1.3     
# [11] MASS_7.3-60         tibble_3.2.1        nloptr_2.0.3        insight_0.19.7      munsell_0.5.0      
# [16] minqa_1.2.6         pillar_1.9.0        TMB_1.9.14          rlang_1.1.2         utf8_1.2.4         
# [21] cli_3.6.2           mgcv_1.9-0          withr_2.5.2         magrittr_2.0.3      grid_4.3.2         
# [26] rstudioapi_0.15.0   lifecycle_1.0.4     nlme_3.1-163        vctrs_0.6.5         glue_1.6.2         
# [31] numDeriv_2016.8-1.1 fansi_1.0.6         colorspace_2.1-0    tools_4.3.2         pkgconfig_2.0.3  

library(ggplot2) #for graphs
library(dplyr)
library(bannerCommenter) # For the section banners (optional)
library(gridExtra)
library(Matrix) # Matrix_1.6-5 or later (lme4 does not work with Matrix_1.6-1.1 or earlier)
library(lme4) # For stats
library(performance) # For modeling checking
library(glmmTMB)
library(DHARMa)

rm(list=ls(all=TRUE))
#Set up Working directory 
setwd("***Pathway to folder containing this R script and all '_accumulated mouse click data.csv' files***")

#Create graph theme
#Create graph theme
Graph.theme<-theme_bw()+
  theme(axis.text=element_text(size=30,colour="black"),
        axis.title=element_text(size=35, colour="black"),
        axis.title.y=element_text(vjust=1),axis.title.x=element_text(vjust=-3),
        plot.title=element_text(size=20),
        plot.margin = unit(c(1,1,1,1), "cm"),#This adjusts the margines of the plot so that bits don't get cut off.
        legend.text=element_text(size=20),
        legend.title=element_text(size=20),
        legend.key.size=unit(1, "cm"),
        legend.background=element_rect(),
        legend.key = element_blank(), #this gets ride of the annoying box around each of
        legend.position= "right",
        #legend.position= "none",
        axis.line = element_line(colour = "black", linewidth =1),
        axis.ticks = element_line(colour = "black"),
        panel.border = element_blank(),
        panel.background = element_rect(fill="white"),
        #plot.background = element_rect(fill = "black", colour="white"),
        panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        strip.background =element_rect(fill="white", colour= "black", linewidth=1),
        strip.text=element_text(face="bold", size=25,colour = "black"))
#axis.text.x = element_text(angle = 45, hjust = 1))

### Function to calculate mean misses per second ###

Cal_miss_per_sec <- function(data) {
  processed_data <- data %>%
    group_by(Subject_ID, Condition_name, Radial_frequency, Amplitude, Condition_repeat_number,Experiment) %>%
    summarise(
      total_misses = sum(TargetClicked == "miss"),
      unique_targets = n_distinct(TargetClicked[TargetClicked != "miss"]),
      trial_duration = ifelse(unique_targets == 5, max(ClickTime),180),
      misses_per_second = total_misses / trial_duration
    ) %>%
    ungroup()
  
  return(processed_data)
}

CustomPalette = c("#648FFF", "#DC267F", "#FE6100", "#FFB000") # Colour blind friendly

############################################################################
############################################################################
###                                                                      ###
###         MEAN NUMBER OF MISSES PER SECOND FOR EACH CONDITION.         ###
###                                                                      ###
############################################################################
############################################################################
banner("Mean number of misses per second for each condition.", emph = TRUE)


##===============================================================
##    Boxplot: Area matched targets (Experiments 1A and 1B)     =
##===============================================================
boxup("Boxplot: Area matched targets (Experiments 1A and 1B)", bandChar = "=")

Exp1A_data <-read.csv("Exp1A_accumulated mouse click data.csv", header=T)
Exp1A_data["Experiment"] <- "Exp 1A"
Exp1A_data <- subset(Exp1A_data, TargetID=="T1") 
processed_Exp1A_data <- Cal_miss_per_sec(Exp1A_data)
Exp1B_data <-read.csv("Exp1B_accumulated mouse click data.csv", header=T)
Exp1B_data["Experiment"] <- "Exp 1B"
Exp1B_data <- subset(Exp1B_data, TargetID=="T1")
processed_Exp1B_data <- Cal_miss_per_sec(Exp1B_data)
AreaMatched_data <- rbind (processed_Exp1A_data, processed_Exp1B_data)

AreaMatched_missesPerS_Box <- ggplot(AreaMatched_data, aes(x=factor(Radial_frequency), y= misses_per_second, fill = factor(Amplitude))) + 
  geom_boxplot(lwd=0.75, fatten =2, outlier.shape=NA) + xlab("Radial frequency") + Graph.theme + scale_fill_manual(values=CustomPalette) +
  geom_point(position=position_jitterdodge(jitter.width = 0.25), shape = 21, size=2, alpha=0.6) +
  scale_y_continuous( name="Number of misses averaged over trial duration in seconds")+#, breaks= seq(0,5,1), limits=c(0,5))+#, expand=c(0,0)) + 
  guides(fill=guide_legend(title="Amplitude")) + ggtitle("Exp 1 and Exp 3")
AreaMatched_missesPerS_Box


##===============================================================
##    Boxplot: Size matched targets (Experiments 2A and 2B)     =
##===============================================================
boxup("Boxplot: Size matched targets (Experiments 2A and 2B)", bandChar = "=")

Exp2A_data <-read.csv("Exp2A_accumulated mouse click data.csv", header=T)
Exp2A_data["Experiment"] <- "Exp 2A"
Exp2A_data <- subset(Exp2A_data, TargetID=="T1")
processed_Exp2A_data <- Cal_miss_per_sec(Exp2A_data)
Exp2B_data <-read.csv("Exp2B_accumulated mouse click data.csv", header=T)
Exp2B_data["Experiment"] <- "Exp 2B"
Exp2B_data <- subset(Exp2B_data, TargetID=="T1")
processed_Exp2B_data <- Cal_miss_per_sec(Exp2B_data)
SizeMatched_data <- rbind (processed_Exp2A_data, processed_Exp2B_data)

SizeMatched_missesPerS_Box <- ggplot(SizeMatched_data, aes(x=factor(Radial_frequency), y= misses_per_second, fill = factor(Amplitude))) + 
  geom_boxplot(lwd=0.75, fatten =2, outlier.shape=NA) + xlab("Radial frequency") + Graph.theme + scale_fill_manual(values=CustomPalette) +
  geom_point(position=position_jitterdodge(jitter.width = 0.25), shape = 21, size=2, alpha=0.6) +
  scale_y_continuous(name="Number of misses averaged over trial duration in seconds") + #, breaks= seq(0,5,1), limits=c(0,5), expand=c(0,0)) + 
  guides(fill=guide_legend(title="Amplitude")) + ggtitle("Exp 2 and Exp 4") 
SizeMatched_missesPerS_Box

##---------------------------------------------------------------
##                          Save graphs                         -
##---------------------------------------------------------------
boxup("Save graphs", bandChar = "-")
today <- Sys.Date()
today <- format(today, "%Y%m%d")
All_box <-grid.arrange(AreaMatched_missesPerS_Box, SizeMatched_missesPerS_Box, ncol = 1)
ggsave(file =paste(today, "Misses per second against RF & Amp boxplots.svg"), device = 'svg', plot = All_box, width = 20, height = 16)


###########################################################################
###########################################################################
###                                                                     ###
###                STATS FOR NUMBER OF MISSES PER SECOND                ###
###                                                                     ###
###########################################################################
###########################################################################
banner("Stats for number of misses per second", emph = TRUE)

#################################################################
##         Stats for RF and amplitude as fixed effects         ##
#################################################################
banner("Stats for RF and amplitude as fixed effects", emph = FALSE)

##===============================================================
##      Stats: Area matched targets (Experiments 1A and 1B)     =
##===============================================================
boxup("Stats: Area matched targets (Experiments 1A and 1B)", bandChar = "=")

Exp1A_data <-read.csv("Exp1A_accumulated mouse click data.csv", header=T)
Exp1A_data["Experiment"] <- "Exp 1A"
Exp1A_data <- subset(Exp1A_data, TargetID=="T1" & Radial_frequency!=0) #*1 
processed_Exp1A_data <- Cal_miss_per_sec(Exp1A_data)
Exp1B_data <-read.csv("Exp1B_accumulated mouse click data.csv", header=T)
Exp1B_data["Experiment"] <- "Exp 1B"
Exp1B_data <- subset(Exp1B_data, TargetID=="T1" & Radial_frequency!=0) #*1
processed_Exp1B_data <- Cal_miss_per_sec(Exp1B_data)
AreaMatched_data <- rbind (processed_Exp1A_data, processed_Exp1B_data)
# *1: The circle arguably makes it so that RF and amp aren't independent of 
# each other because one cannot be zero without the other also being zero.

# Fit the full Zero-Inflated mixed effects Gamma Model using glmmTMB.
zig_model <- glmmTMB(misses_per_second ~ Radial_frequency * Amplitude + (1|Experiment/Subject_ID),
                     zi = ~ Radial_frequency * Amplitude,
                     family = ziGamma(link = "log"),
                     data = AreaMatched_data)

simulationOutput<- simulateResiduals(fittedModel = zig_model)
plot(simulationOutput)
# Some of the automatic tests come back as significant. However, this does not necessarily mean there is a problem or that
# the model is not usable. As explained in the vignette for the DHARMa package (https://cran.r-project.org/web/packages/DHARMa/vignettes/DHARMa.html),
# significance in hypothesis tests depends on at least 2 factors, the strength of the signal, and the number of data points. 
# If you have a lot of data points (which we do), residual diagnostics will nearly inevitably become significant, because 
# having a perfectly fitting model is very unlikely. Therefore, it is important to inspect the plots and use the magnitude to 
# decide if there is a problem. Based on visual inspection, the the model is acceptable and the best fit to the data that is 
# achievable. This also applies to the models below that form part of the LRT. 

#***************************************************************
# Test significance of RF and amplitude interaction 
zig_model_2 <- glmmTMB(misses_per_second ~ Radial_frequency + Amplitude + (1|Experiment/Subject_ID),
                     zi = ~ Radial_frequency + Amplitude,
                     family = ziGamma(link = "log"),
                     data = AreaMatched_data)
simulationOutput<- simulateResiduals(fittedModel = zig_model_2)
plot(simulationOutput)
anova(zig_model,zig_model_2) # Chi(2) = 34.412, p= 3.369e-08 ***
#***************************************************************


##===============================================================
##      Stats: Size matched targets (Experiments 2A and 2B)     =
##===============================================================
boxup("Stats: Size matched targets (Experiments 2A and 2B)", bandChar = "=")
# For the sized matched experiment I go straight into using the zero-inflated mixed effect gamma model

Exp2A_data <-read.csv("Exp2A_accumulated mouse click data.csv", header=T)
Exp2A_data["Experiment"] <- "Exp 2A"
Exp2A_data <- subset(Exp2A_data, TargetID=="T1"& Radial_frequency!=0) #*1
processed_Exp2A_data <- Cal_miss_per_sec(Exp2A_data)
Exp2B_data <-read.csv("Exp2B_accumulated mouse click data.csv", header=T)
Exp2B_data["Experiment"] <- "Exp 2B"
Exp2B_data <- subset(Exp2B_data, TargetID=="T1" & Radial_frequency!=0) #*1
processed_Exp2B_data <- Cal_miss_per_sec(Exp2B_data)
SizeMatched_data <- rbind (processed_Exp2A_data, processed_Exp2B_data)

# Fit the full Zero-Inflated mixed effects Gamma Model using glmmTMB.
zig_model <- glmmTMB(misses_per_second ~ Radial_frequency * Amplitude + (1|Experiment/Subject_ID),
                     zi = ~ Radial_frequency * Amplitude,
                     family = ziGamma(link = "log"),
                     data = SizeMatched_data)

simulationOutput<- simulateResiduals(fittedModel = zig_model)
plot(simulationOutput)
# Some of the automatic tests come back as significant. However, this does not necessarily mean there is a problem or that
# the model is not usable. As explained in the vignette for the DHARMa package (https://cran.r-project.org/web/packages/DHARMa/vignettes/DHARMa.html),
# significance in hypothesis tests depends on at least 2 factors, the strength of the signal, and the number of data points. 
# If you have a lot of data points (which we do), residual diagnostics will nearly inevitably become significant, because 
# having a perfectly fitting model is very unlikely. Therefore, it is important to inspect the plots and use the magnitude to 
# decide if there is a problem. Based on visual inspection, although not perfect, the the model is acceptable enough and the 
# best fit to the data that is achievable. This also applies to the models below that form part of the LRT. 

#***************************************************************
# Test significance of RF and amplitude interaction 
zig_model_2 <- glmmTMB(misses_per_second ~ Radial_frequency + Amplitude + (1|Experiment/Subject_ID),
                       zi = ~ Radial_frequency + Amplitude,
                       family = ziGamma(link = "log"),
                       data = SizeMatched_data)
simulationOutput<- simulateResiduals(fittedModel = zig_model_2)
plot(simulationOutput)
anova(zig_model,zig_model_2) # Chi(2) = 68.552, p = 1.3e-15 ***
#***************************************************************
