
# This is a copy of the R script used to generate the figures for, and perform the statistical analysis on, 
# the Restricted Mean Survival Time (RMST) data from the experiments reported in Smithers et al. 
# (under review). For an explanation of the statistical analysis conducted, please refer to the 'Analysis 
# of Search Time Data' section of the methods in the main manuscript.

# Script written by Dr Samuel P. Smithers, Northeastern University, 2023-2025
# Last edited March 2025

# Corresponding author: PJB (p.bex@northeastern.edu)  

# R version 4.3.2 (2023-10-31 ucrt)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 19045)
# 
# Matrix products: default
# 
# 
# locale:
#   [1] LC_COLLATE=English_United States.utf8  LC_CTYPE=English_United States.utf8    LC_MONETARY=English_United States.utf8
# [4] LC_NUMERIC=C                           LC_TIME=English_United States.utf8    
# 
# time zone: America/New_York
# tzcode source: internal
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] bannerCommenter_1.0.0 performance_0.10.8    lme4_1.1-35.1         Matrix_1.6-5          svglite_2.1.3        
# [6] gridExtra_2.3         scales_1.3.0          ggplot2_3.4.4        
# 
# loaded via a namespace (and not attached):
#   [1] gtable_0.3.4      dplyr_1.1.4       compiler_4.3.2    tidyselect_1.2.0  Rcpp_1.0.11       splines_4.3.2    
# [7] systemfonts_1.0.5 boot_1.3-28.1     lattice_0.21-9    R6_2.5.1          generics_0.1.3    MASS_7.3-60      
# [13] tibble_3.2.1      nloptr_2.0.3      insight_0.19.7    munsell_0.5.0     minqa_1.2.6       pillar_1.9.0     
# [19] rlang_1.1.2       utf8_1.2.4        datawizard_0.9.1  cli_3.6.2         withr_2.5.2       magrittr_2.0.3   
# [25] grid_4.3.2        rstudioapi_0.15.0 lifecycle_1.0.4   nlme_3.1-163      vctrs_0.6.5       glue_1.6.2       
# [31] farver_2.1.1      fansi_1.0.6       colorspace_2.1-0  tools_4.3.2       pkgconfig_2.0.3  


## Required dependencies/packages. 
library(ggplot2) #for graphs
library(scales) # to access break formatting functions
library(gridExtra) # For arranging plots for saving (optional)
library(svglite) # For saving plots (optional)
library(Matrix) # Matrix_1.6-5 or later (lme4 does not work with Matrix_1.6-1.1 or earlier)
library(lme4) # For stats
library(performance) # For modeling checking
library(bannerCommenter) # For the section banners (optional)


#Set up Working directory 
setwd("***Pathway to folder containing this R script and all '_RMST_data.csv' files***")

############################################################################
############################################################################
###                                                                      ###
###                         STATISTICAL ANALYSIS                         ###
###                                                                      ###
############################################################################
############################################################################
banner("Statistical analysis", emph = TRUE)
 
rm(list=ls(all=TRUE))

#################################################################
##         Stats for RF and amplitude as fixed effects         ##
#################################################################
banner("Stats for RF and amplitude as fixed effects", emph = FALSE)


##===============================================================
##      Stats: Area matched targets (Experiments 1A and 1B)     =
##===============================================================
boxup("Stats: Area matched targets (Experiments 1A and 1B)", bandChar = "=")

Exp1A_data <-read.csv("Exp1A_RMST_data.csv", header=T)
Exp1A_data["Experiment"] <- "Exp 1A"
Exp1A_data<- subset(Exp1A_data, Radial_frequency!=0) # Remove circle condition *1
Exp1B_data <-read.csv("Exp1B_RMST_data.csv", header=T)
Exp1B_data["Experiment"] <- "Exp 1B"
Exp1B_data<- subset(Exp1B_data, Radial_frequency!=0) # Remove circle condition *1
AreaMatched_data <- rbind (Exp1A_data,Exp1B_data)
# *1: The circle makes it so that RF and amp aren't independent of each
# other because one cannot be zero without the other also being zero.

# We first use Cleveland dotplots to check for overly influential outliers. 
# Note that we plot log(RMST) here as it is natural log transformed in the model below.
dotchart(log(AreaMatched_data$RMST),
         ylab = "Order of observations", xlab = "RMST",
         main = "Cleveland dotplot")
# Based on this plot none of the values are extreme. 

# First we check that we do not have an issue of multicollinearity of model terms. For this we use
# the variance inflation factor (VIF) which is a measure to analyze the magnitude of multicollinearity 
# of model terms. A VIF less than 5 indicates a low correlation of that predictor with other predictors. 
# A value between 5 and 10 indicates a moderate correlation, while VIF values larger than 10 are a sign 
# of high, not tolerable correlation of model predictors.

# First we fit a model with the effects of interest without the interaction. We do not include the interaction
# because it can result in incorrectly inflated VIF values. 
Model<-lmer(log(RMST)~ Radial_frequency+Amplitude+ (1|Experiment/Subject_ID), data= AreaMatched_data, REML=TRUE)

# We then check for multicollinearity of model terms. 
check_collinearity(Model) # https://easystats.github.io/performance/reference/check_collinearity.html
# The check finds low correlation between predictors indicating that multicollinearity is not a problem.

# We can now fit the full mixed model containing fixed effects and pairwise interaction of the fixed effects 
M0<-lmer(log(RMST)~ Radial_frequency*Amplitude + (1|Experiment/Subject_ID), data= AreaMatched_data, REML=TRUE) 
# The response variable is natural log transformed to correct for positive skew in the distribution of residuals.

#Check normality and homoscedasticity
sresid <- resid(M0)  # Extract the residuals
hist(sresid) # Plot histogram and QQ plot to look at distribution of residuals
qqnorm(sresid)
qqline(sresid)
fitted.glmm <- fitted(M0, level=1) # Extract the fitted (predicted) values
plot(sresid ~ fitted.glmm) # Check for homoscedasticity 
plot(M0)
# The plot of residuals vs fitted suggests that there may be some slight heterogeneity (i.e. a
# weak pattern in the residual), however it is likely not enough to be problematic. To be sure 
# it is not a problem we use the performance package to check for heterogeneity statistically. 
check_heteroscedasticity(M0)
# This confirms that heterogeneity is not a problem. 

#***************************************************************
# Test significance of the interaction using a likelihood ratio rest (LRT) to compare 
# the full model with the same model but with the interaction removed. 
m1 <-update(M0,~.-Radial_frequency:Amplitude)
anova(M0, m1) # Chi(1) = 71.745, p= < 2.2e-16 ***
# The interaction is significant
#***************************************************************
rm(list=ls(all=TRUE))

##===============================================================
##      Stats: Size matched targets (Experiments 2 and 4)       =
##===============================================================
boxup("Stats: Size matched targets (Experiments 2 and 4)", bandChar = "=")

Exp2A_data <-read.csv("Exp2A_RMST_data.csv", header=T)
Exp2A_data["Experiment"] <- "Exp 2A"
Exp2A_data<- subset(Exp2A_data, Radial_frequency!=0) # Remove circle condition *1
Exp2B_data <-read.csv("Exp2B_RMST_data.csv", header=T)
Exp2B_data["Experiment"] <- "Exp 2B"
Exp2B_data<- subset(Exp2B_data, Radial_frequency!=0) # Remove circle condition *1
SizeMatched_data <- rbind (Exp2A_data,Exp2B_data)

# We first use Cleveland dotplots to check for overly influential outliers. 
# Note that we plot log(RMST) here as it is natural log transformed in the model below.
dotchart(log(SizeMatched_data$RMST),
         ylab = "Order of observations", xlab = "RMST",
         main = "Cleveland dotplot")
# Based on this plot none of the values are extreme. 

# First we check that we do not have an issue of multicollinearity of model terms. For this we use
# the variance inflation factor (VIF) which is a measure to analyze the magnitude of multicollinearity 
# of model terms. A VIF less than 5 indicates a low correlation of that predictor with other predictors. 
# A value between 5 and 10 indicates a moderate correlation, while VIF values larger than 10 are a sign 
# of high, not tolerable correlation of model predictors.

# First we fit a model with the effects of interest without the interactions. We do not include the interaction
# because it can result in incorrectly inflated VIF values. 
Model<-lmer(log(RMST)~ Radial_frequency+Amplitude + (1|Experiment/Subject_ID), data= SizeMatched_data, REML=TRUE) 

# We then check for multicollinearity of model terms. 
check_collinearity(Model) # https://easystats.github.io/performance/reference/check_collinearity.html
# The check finds low correlation between predictors indicating that multicollinearity is not a problem.

# We can now fit the full mixed model containing fixed effects and pairwise interaction of the fixed effects 
M0<-lmer(log(RMST)~ Radial_frequency*Amplitude + (1|Experiment/Subject_ID), data= SizeMatched_data, REML=TRUE) 
# The response variable is natural log transformed to correct for positive skew in the distribution of residuals.

#Check normality and homoscedasticity
sresid <- resid(M0)  # Extract the residuals
hist(sresid) # Plot histogram and QQ plot to look at distribution of residuals
qqnorm(sresid)
qqline(sresid)
fitted.glmm <- fitted(M0, level=1) # Extract the fitted (predicted) values
plot(sresid ~ fitted.glmm) # Check for homoscedasticity 
plot(M0)
# The plot of residuals vs fitted suggests that there may be some slight heterogeneity (i.e. a
# weak pattern in the residual), however it is likely not enough to be problematic. To be sure
# it is not a problem we use the performance package to check for heterogeneity statistically. 
check_heteroscedasticity(M0)
# This confirms that heterogeneity is not a problem. 

#***************************************************************
# Test significance of the interaction using a likelihood ratio rest (LRT) to compare 
# the full model with the same model but with the interaction removed. 
m1 <-update(M0,~.-Radial_frequency:Amplitude)
anova(M0, m1) # Chi(1) = 88.122, p= < 2.2e-16 ***
# The interaction is significant
#***************************************************************
rm(list=ls(all=TRUE))


#####################################################################
##  Stats for RF, amplitude and perimeter length as fixed effects  ##
#####################################################################
banner("Stats for RF, amplitude and perimeter length as fixed effects", emph = FALSE)

##===============================================================
##      Stats: Area matched targets (Experiments 1 and 3)       =
##===============================================================
boxup("Stats: Area matched targets (Experiments 1 and 3)", bandChar = "=")

Exp1A_data <-read.csv("Exp1A_RMST_data.csv", header=T)
Exp1A_data["Experiment"] <- "Exp 1A"
Exp1A_data<- subset(Exp1A_data, Radial_frequency!=0) # Remove circle condition *1
Exp1B_data <-read.csv("Exp1B_RMST_data.csv", header=T)
Exp1B_data["Experiment"] <- "Exp 1B"
Exp1B_data<- subset(Exp1B_data, Radial_frequency!=0) # Remove circle condition *1
AreaMatched_data <- rbind (Exp1A_data,Exp1B_data)

# We first use Cleveland dotplots to check for overly influential outliers. 
# Note that we plot log(RMST) here as it is natural log transformed in the model below.
dotchart(log(AreaMatched_data$RMST),
         ylab = "Order of observations", xlab = "RMST",
         main = "Cleveland dotplot")
# Based on this plot none of the values are extreme. 

# First we check that we do not have an issue of multicollinearity of model terms. For this we use
# the variance inflation factor (VIF), which is a measure to analyze the magnitude of multicollinearity 
# of model terms. A VIF less than 5 indicates a low correlation of that predictor with other predictors. 
# A value between 5 and 10 indicates a moderate correlation, while VIF values larger than 10 are a sign 
# of high, not tolerable correlation of model predictors.

# First we fit a model with the effects of interest without any interactions. We do not include the interactions
# because they can result in incorrectly inflated VIF values. 
Model<-lmer(log(RMST)~ Radial_frequency+Amplitude+perimeter_cm+ (1|Experiment/Subject_ID), data= AreaMatched_data, REML=TRUE) 

# We then check for multicollinearity of model terms. 
check_collinearity(Model) # https://easystats.github.io/performance/reference/check_collinearity.html
# The check finds low correlation between predictors indicating that multicollinearity is not a problem.

# We can now fit the full mixed model containing fixed effects and all possible pairwise interactions of the fixed effects 
M0<-lmer(log(RMST)~ Radial_frequency*Amplitude*perimeter_cm+ (1|Experiment/Subject_ID), data= AreaMatched_data, REML=TRUE) 
# The response variable is natural log transformed to correct for positive skew in the distribution of residuals.

#Check normality and homoscedasticity
sresid <- resid(M0)  # Extract the residuals
hist(sresid) # Plot histogram and QQ plot to look at distribution of residuals
qqnorm(sresid)
qqline(sresid)
fitted.glmm <- fitted(M0, level=1) # Extract the fitted (predicted) values
plot(sresid ~ fitted.glmm) # Check for homoscedasticity 
plot(M0)
# The plot of residuals vs fitted suggests that there may be some slight heterogeneity (i.e. a
# weak pattern in the residual), however it is likely not enough to be problematic. To be sure 
# it is not a problem we use the performance package to check for heterogeneity statistically. 
check_heteroscedasticity(M0)
# This confirms that heterogeneity is not a problem. 

# Test significance of the three-way interaction using a likelihood ratio rest (LRT) to compare 
# the full model with the same model but with the interaction removed. 
m1 <-update(M0,~.-Radial_frequency:Amplitude:perimeter_cm)
anova(M0, m1) # Chi(1) = 1.2547, p= 0.2627
# The interaction is not significant so is removed

#***************************************************************
m2a <-update(m1,~.-Radial_frequency:Amplitude)
anova(m2a, m1) # Chi(1) = 42.379, p = 7.52e-11 ***
m2b <-update(m1,~.-Radial_frequency:perimeter_cm)
anova(m2b, m1) # Chi(1) = 23.626, p = 1.17e-06 ***
m2c <-update(m1,~.-Amplitude:perimeter_cm)
anova(m2c, m1) # Chi(1) = 48.577 , p = 3.176e-12 ***
#***************************************************************
rm(list=ls(all=TRUE))

##===============================================================
##      Stats: Size matched targets (Experiments 2 and 4)       =
##===============================================================
boxup("Stats: Size matched targets (Experiments 2 and 4)", bandChar = "=")

Exp2A_data <-read.csv("Exp2A_RMST_data.csv", header=T)
Exp2A_data["Experiment"] <- "Exp 2A"
Exp2A_data<- subset(Exp2A_data, Radial_frequency!=0) # Remove circle condition *1
Exp2B_data <-read.csv("Exp2B_RMST_data.csv", header=T)
Exp2B_data["Experiment"] <- "Exp 2B"
Exp2B_data<- subset(Exp2B_data, Radial_frequency!=0) # Remove circle condition *1
SizeMatched_data <- rbind (Exp2A_data,Exp2B_data)

# We first use Cleveland dotplots to check for overly influential outliers. 
# Note that we plot log(RMST) here as it is natural log transformed in the model below.
dotchart(log(SizeMatched_data$RMST),
         ylab = "Order of observations", xlab = "RMST",
         main = "Cleveland dotplot")
# Based on this plot none of the values are extreme. 

# First we check that we do not have an issue of multicollinearity of model terms. For this we use
# the variance inflation factor (VIF), which is a measure to analyze the magnitude of multicollinearity 
# of model terms. A VIF less than 5 indicates a low correlation of that predictor with other predictors. 
# A value between 5 and 10 indicates a moderate correlation, while VIF values larger than 10 are a sign 
# of high, not tolerable correlation of model predictors.

# First we fit a model with the effects of interest without any interactions. We do not include the interactions
# because they can result in incorrectly inflated VIF values. 
Model<-lmer(log(RMST)~ Radial_frequency+Amplitude+perimeter_cm+ (1|Experiment/Subject_ID), data= SizeMatched_data, REML=TRUE) 

# We then check for multicollinearity of model terms. 
check_collinearity(Model) # https://easystats.github.io/performance/reference/check_collinearity.html
# The check finds low correlation between predictors indicating that multicollinearity is not a problem.

# We can now fit the full mixed model containing fixed effects and all possible pairwise interactions of the fixed effects 
M0<-lmer(log(RMST)~ Radial_frequency*Amplitude*perimeter_cm+ (1|Experiment/Subject_ID), data= SizeMatched_data, REML=TRUE) 
# The response variable is natural log transformed to correct for positive skew in the distribution of residuals.

#Check normality and homoscedasticity
sresid <- resid(M0)  # Extract the residuals
hist(sresid) # Plot histogram and QQ plot to look at distribution of residuals
qqnorm(sresid)
qqline(sresid)
fitted.glmm <- fitted(M0, level=1) # Extract the fitted (predicted) values
plot(sresid ~ fitted.glmm) # Check for homoscedasticity 
plot(M0)
# The plot of residuals vs fitted suggests that there may be some slight heterogeneity (i.e. a
# weak pattern in the residual), however it is likely not enough to be problematic. To be sure
# it is not a problem we use the performance package to check for heterogeneity statistically. 
check_heteroscedasticity(M0)
# This confirms that heterogeneity is not a problem. 

# Test significance of the three-way interaction using a likelihood ratio rest (LRT) to compare 
# the full model with the same model but with the interaction removed. 
m1 <-update(M0,~.-Radial_frequency:Amplitude:perimeter_cm)
anova(M0, m1) # Chi(1) = 3.6841, p=  0.05493
# The interaction is not significant so is removed

#***************************************************************
m2a <-update(m1,~.-Radial_frequency:Amplitude)
anova(m2a, m1) # Chi(1) = 26.008, p = 3.4e-07 ***
m2b <-update(m1,~.-Radial_frequency:perimeter_cm)
anova(m2b, m1) # Chi(1) = 6.6347, p = 0.01 *
m2c <-update(m1,~.-Amplitude:perimeter_cm)
anova(m2c, m1) # Chi(1) = 6.4767, p = 0.01093 *
#***************************************************************
rm(list=ls(all=TRUE))


############################################################################
############################################################################
###                                                                      ###
###                                GRAPHS                                ###
###                                                                      ###
############################################################################
############################################################################
banner("Graphs", emph = TRUE)

#Create graph theme
Graph.theme<-theme_bw()+
  theme(axis.text=element_text(size=30,colour="black"),
        axis.title=element_text(size=35, colour="black"),
        axis.title.y=element_text(vjust=1),axis.title.x=element_text(vjust=-3),
        plot.title=element_text(size=20),
        plot.margin = unit(c(1,1,1,1), "cm"),
        legend.text=element_text(size=20),
        legend.title=element_text(size=20),
        legend.key.size=unit(1, "cm"),
        legend.background=element_rect(),
        legend.key = element_blank(), 
        legend.position= "right",
        #legend.position= "none", # Uncomment to remove the legend/key
        axis.line = element_line(colour = "black", linewidth =1),
        axis.ticks = element_line(colour = "black"),
        panel.border = element_blank(),
        panel.background = element_rect(fill="white"),
        panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        strip.background =element_rect(fill="white", colour= "black", linewidth=1),
        strip.text=element_text(face="bold", size=25,colour = "black"))


######################################################################################
##  Figure 4- Boxplots for RMST against radial frequency with legend for amplitude  ##
######################################################################################
banner("Figure 4- Boxplots for RMST against radial frequency with legend for amplitude", emph = FALSE)
CustomPalette = c("#648FFF", "#DC267F", "#FE6100", "#FFB000") # Colour blind friendly

##====================================================================
##  Boxplot: Area matched targets (Figure 4A- Experiments 1A and 1B)   =
##====================================================================
boxup("Boxplot: Area matched targets (Figure 4A- Experiments 1 and 3)", bandChar = "=")

Exp1A_data <-read.csv("Exp1A_RMST_data.csv", header=T)
Exp1B_data <-read.csv("Exp1B_RMST_data.csv", header=T)
N1 <- nrow(unique(Exp1A_data["Subject_ID"]))
print(N1) # N=23
N3 <- nrow(unique(Exp1B_data["Subject_ID"]))
print(N3) # N=20
AreaMatched_data <- rbind (Exp1A_data,Exp1B_data)

AreaMatched_RMST_Box <- ggplot(AreaMatched_data, aes(x=factor(Radial_frequency), y= RMST, fill = factor(Amplitude))) + 
  geom_boxplot(lwd=0.75, fatten =2, outlier.shape=NA) + xlab("Radial frequency") + Graph.theme + scale_fill_manual(values=CustomPalette) +
  geom_point(position=position_jitterdodge(jitter.width = 0.25), shape = 21, size=2, alpha=0.6) +
  scale_y_continuous(breaks= seq(0,100,20), limits=c(0,105), expand=c(0,0),  name="Restricted Mean Survival Time (s)") + 
  guides(fill=guide_legend(title="Amplitude")) + ggtitle("Exp 1 and Exp 3")
AreaMatched_RMST_Box


##====================================================================
##  Boxplot: Size matched targets (Figure 4B- Experiments 2A and 2B)   =
##====================================================================
boxup("Boxplot: Size matched targets (Figure 4B- Experiments 2A and 2B)", bandChar = "=")

Exp2A_data <-read.csv("Exp2A_RMST_data.csv", header=T)
Exp2B_data <-read.csv("Exp2B_RMST_data.csv", header=T)
N2 <- nrow(unique(Exp2A_data["Subject_ID"]))
print(N2) # N=26
N4 <- nrow(unique(Exp2B_data["Subject_ID"]))
print(N4) # N=20
SizeMatched_data <- rbind (Exp2A_data,Exp2B_data)

SizeMatched_RMST_Box <- ggplot(SizeMatched_data, aes(x=factor(Radial_frequency), y= RMST, fill = factor(Amplitude))) + 
  geom_boxplot(lwd=0.75, fatten =2, outlier.shape=NA) + xlab("Radial frequency") + Graph.theme + scale_fill_manual(values=CustomPalette) +
  geom_point(position=position_jitterdodge(jitter.width = 0.25), shape = 21, size=2, alpha=0.6) +
  scale_y_continuous(breaks= seq(0,100,20), limits=c(0,105), expand=c(0,0),  name="Restricted Mean Survival Time (s)") + 
  guides(fill=guide_legend(title="Amplitude")) + ggtitle("Exp 2 and Exp 4")
SizeMatched_RMST_Box

##---------------------------------------------------------------
##                          Save graphs                         -
##---------------------------------------------------------------
boxup("Save graphs", bandChar = "-")
today <- Sys.Date()
today <- format(today, "%Y%m%d")
All_box <-grid.arrange(AreaMatched_RMST_Box, SizeMatched_RMST_Box, ncol = 1)
ggsave(file =paste(today, "Copy of Figure 4 from Smithers et al.svg"), device = 'svg', plot = All_box, width = 20, height = 15)


#####################################################################################################
##  Figure 5- Scatter graphs for RMST against perimeter with color for RF and shape for amplitude  ##
#####################################################################################################
banner("Figure 5- Scatter graphs for RMST against perimeter with color for RF and shape for amplitude", emph = FALSE)
# Scatter plot with mean +- CI and local polynomial regression curve.

CustomPalette = c("#000000", "#00F1FF", "#1A85FF", "#0A00FF", "#FF20B7", "#FFFA00","#FF9000","#FF0000") # Colour blind friendly

##=========================================================================
##  Scatter plot: Area matched targets (Figure 5A- Experiments 1A and 1B)   =
##=========================================================================
boxup("Scatter plot: Area matched targets (Figure 5A- Experiments 1A and 1B)", bandChar = "=")
                          
Exp1A_data <-read.csv("Exp1A_RMST_data.csv", header=T)
Exp1B_data <-read.csv("Exp1B_RMST_data.csv", header=T)
AreaMatched_data <- rbind (Exp1A_data,Exp1B_data)

AreaMatched_RMST_vs_perimeter_scatter <- ggplot(AreaMatched_data, aes(x=perimeter_cm, y= RMST)) +
  xlab("Perimeter (cm)") +  ylab("Restricted Mean Survival Time (s)") + Graph.theme + 
  geom_point(size=2, alpha=0.6, aes(fill = factor(Radial_frequency), shape = factor(Amplitude))) +
  stat_summary(geom = "point", fun = "mean", size = 4, aes(fill = factor(Radial_frequency), shape = factor(Amplitude))) + 
  stat_summary(geom = "errorbar", fun.data = "mean_se", width=0.02, linewidth=0.5,  aes(color = factor(Radial_frequency))) + 
  scale_shape_manual(values = c(23,21,24,22)) +
  scale_fill_manual(values=CustomPalette) +  scale_color_manual(values=CustomPalette) +
  geom_smooth(method = 'loess', formula = 'y ~ x', color = "black") + # Add local polynomial regression curve
  scale_y_log10(breaks=c(1,5,10,50,100)) + scale_x_log10(breaks=c(10,20,30,40,50)) + 
  annotation_logticks(outside = TRUE, size=0.8) + coord_cartesian(clip = "off") +
  guides(fill=guide_legend(override.aes=list(shape=21), title="Radial frequency"), 
         shape=guide_legend(override.aes=list(shape=c(18,16,17,15)), title="Amplitude"), 
         color = 'none') + ggtitle("Area matched: Exp 1 & 3") 
AreaMatched_RMST_vs_perimeter_scatter

##=========================================================================
##  Scatter plot: Size matched targets (Figure 5B- Experiments 2A and 2B)   =
##=========================================================================
boxup("Scatter plot: Size matched targets (Figure 5B- Experiments 2A and 2B)", bandChar = "=") 

Exp2A_data <-read.csv("Exp2A_RMST_data.csv", header=T)
Exp2B_data <-read.csv("Exp2B_RMST_data.csv", header=T)
SizeMatched_data <- rbind (Exp2A_data,Exp2B_data)

SizeMatched_RMST_vs_perimeter_scatter <- ggplot(SizeMatched_data, aes(x=perimeter_cm, y= RMST)) +
  xlab("Perimeter (cm)") +  ylab("Restricted Mean Survival Time (s)") + Graph.theme +
  geom_point(size=2, alpha=0.6, aes(fill = factor(Radial_frequency), shape = factor(Amplitude))) +
  stat_summary(geom = "point", fun = "mean", size = 4, aes(fill = factor(Radial_frequency), shape = factor(Amplitude))) + 
  stat_summary(geom = "errorbar", fun.data = "mean_se", width=0.02, linewidth=0.5,  aes(color = factor(Radial_frequency))) + 
  scale_shape_manual(values = c(23,21,24,22)) +
  scale_fill_manual(values=CustomPalette) +  scale_color_manual(values=CustomPalette) +
  geom_smooth(method = 'loess', formula = 'y ~ x', color = "black") + # Add local polynomial regression curve
  scale_y_log10(breaks=c(1,5,10,50,100)) + scale_x_log10(breaks=c(10,20,30)) + 
  annotation_logticks(outside = TRUE, size=0.8) + coord_cartesian(clip = "off") +
  guides(fill=guide_legend(override.aes=list(shape=21), title="Radial frequency"), 
         shape=guide_legend(override.aes=list(shape=c(18,16,17,15)), title="Amplitude"), 
         color = 'none') + ggtitle("Angular size matched: Exp 2 & 4") 
SizeMatched_RMST_vs_perimeter_scatter

##---------------------------------------------------------------
##                          Save graphs                         -
##---------------------------------------------------------------
boxup("Save graphs", bandChar = "-")
today <- Sys.Date()
today <- format(today, "%Y%m%d")
All_scatter <-grid.arrange(AreaMatched_RMST_vs_perimeter_scatter,SizeMatched_RMST_vs_perimeter_scatter, ncol = 1)
ggsave(file =paste(today, "Copy of Figure 5 from Smithers et al.svg"), device = 'svg', plot = All_scatter, width = 20, height = 15)
