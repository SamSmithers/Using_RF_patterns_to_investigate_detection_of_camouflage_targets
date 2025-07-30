# This is a copy of the R script used to generate supplementary figures S2-S5 from Smithers et al.
# (under review).

# Script written by Dr Samuel P. Smithers, Northeastern University, 2023-2025
# Last edited July 2025

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
#   [1] bannerCommenter_1.0.0 viridis_0.6.4         viridisLite_0.4.2     ggplot2_3.4.4        
# 
# loaded via a namespace (and not attached):
#   [1] Matrix_1.6-5      gtable_0.3.4      dplyr_1.1.4       compiler_4.3.2    tidyselect_1.2.0 
# [6] Rcpp_1.0.11       DHARMa_0.4.6      gridExtra_2.3     splines_4.3.2     scales_1.3.0     
# [11] boot_1.3-28.1     lattice_0.21-9    R6_2.5.1          generics_0.1.3    MASS_7.3-60      
# [16] tibble_3.2.1      nloptr_2.0.3      munsell_0.5.0     minqa_1.2.6       pillar_1.9.0     
# [21] TMB_1.9.14        rlang_1.1.2       utf8_1.2.4        cli_3.6.2         withr_2.5.2      
# [26] magrittr_2.0.3    grid_4.3.2        rstudioapi_0.15.0 lme4_1.1-35.1     lifecycle_1.0.4  
# [31] nlme_3.1-163      vctrs_0.6.5       glue_1.6.2        fansi_1.0.6       colorspace_2.1-0 
# [36] tools_4.3.2       pkgconfig_2.0.3  

library(ggplot2) #for graphs
library(viridis)
library(bannerCommenter) # For the section banners (optional)

rm(list=ls(all=TRUE))

#Set up Working directory 
setwd("***Pathway to folder containing this R script and all '_accumulated_rawdata.csv' files***")

#Create graph theme
Graph.theme<-theme_bw()+
  theme(axis.text=element_text(size=18,colour="black"),
        axis.title=element_text(size=25, colour="black"),
        axis.title.y=element_text(vjust=1),axis.title.x=element_text(vjust=-3),
        plot.margin = unit(c(1,1,1,1), "cm"),
        legend.text=element_text(size=20),
        legend.title=element_text(size=25),
        legend.key.size=unit(1, "cm"),
        legend.background=element_rect(),
        legend.key = element_blank(),
        legend.position="right",
        axis.line = element_line(colour = "black", linewidth=1),
        axis.ticks = element_line(colour = "black"),
        #panel.border = element_blank(),
        panel.background = element_rect(fill="white"),
        #panel.grid.major.x=element_blank(),
        panel.grid.minor=element_blank(),
        strip.background =element_rect(fill="white", colour= "black", linewidth=1),
        strip.text=element_text(face="bold", size=25,colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1))



#######################################################################################
#######################################################################################
###                                                                                 ###
###  PLOTS OF TARGET SCREEN POSITION (X & Y) AGAINST DETECTION TIME (COLOUR SCALE)  ###
###                                                                                 ###
#######################################################################################
#######################################################################################
banner("Plots of target screen position (x & Y) against detection time (colour scale)", emph = TRUE)

##===============================================================
##                        Experiment 1A                         =
##===============================================================
boxup("Experiment 1A", bandChar = "=") 
  
### Exp 1A data                     
data <-read.csv("Exp1A_accumulated_rawdata.csv", header=T)

temp1<- subset(data, TimeOfEvent < 180)
temp2<- subset(data, TimeOfEvent == 180)

TargetPosVsDetection_exp1A<- ggplot() +
  # Layer for all targets that were found
  geom_point(data=temp1, aes(x=TargetPos_X, y= TargetPos_Y, fill= TimeOfEvent),shape = 21, size=2) + 
  scale_fill_viridis(name="Detection\ntime (s)", option="viridis",trans = "log10", breaks = c(1, 10, 100, 180), labels = c("1", "10", "100", "180")) +
  # Layer for all targets that were NOT found
  geom_point(data=temp2, aes(x=TargetPos_X, y= TargetPos_Y, color= "Target\n   not\n found"), shape = 20, size=3, show.legend = TRUE) +  # Ensure this layer shows in the legend
  # Add a custom legend entry for temp2 points
  scale_color_manual(values = c("Target\n   not\n found" = "red")) +  
  Graph.theme +
  facet_grid(Amplitude ~ Radial_frequency) + xlab("X target position (deg)") + ylab("Y target position (deg)") +
  scale_y_continuous(breaks= seq(-20,20,10), limits=c(-20,20), expand=c(0,0)) +
  scale_x_continuous(breaks= seq(-20,20,10), limits=c(-20,20), expand=c(0,0))
 # coord_fixed(ratio=1) #This ensures that the x and Y axis are the same length
TargetPosVsDetection_exp1A

##===============================================================
##                        Experiment 1B                         =
##===============================================================
boxup("Experiment 1B", bandChar = "=") 

### Exp 1B data     
data <-read.csv("Exp1B_accumulated_rawdata.csv", header=T)

temp1<- subset(data, TimeOfEvent < 180)
temp2<- subset(data, TimeOfEvent == 180)

TargetPosVsDetection_exp1B<- ggplot() +
  # Layer for all targets that were found
  geom_point(data=temp1, aes(x=TargetPos_X, y= TargetPos_Y, fill= TimeOfEvent),shape = 21, size=2) + 
  scale_fill_viridis(name="Detection\ntime (s)", option="viridis",trans = "log10", breaks = c(1, 10, 100, 180), labels = c("1", "10", "100", "180")) +
  # Layer for all targets that were NOT found
  geom_point(data=temp2, aes(x=TargetPos_X, y= TargetPos_Y, color= "Target\n   not\n found"), shape = 20, size=3, show.legend = TRUE) +  # Ensure this layer shows in the legend
  # Add a custom legend entry for temp2 points
  scale_color_manual(values = c("Target\n   not\n found" = "red")) +  
  Graph.theme +
  facet_grid(Amplitude ~ Radial_frequency) + xlab("X target position (deg)") + ylab("Y target position (deg)") +
  scale_y_continuous(breaks= seq(-20,20,10), limits=c(-20,20), expand=c(0,0)) +
  scale_x_continuous(breaks= seq(-20,20,10), limits=c(-20,20), expand=c(0,0))
# coord_fixed(ratio=1) #This ensures that the x and Y axis are the same length
TargetPosVsDetection_exp1B

##===============================================================
##                        Experiment 2A                         =
##===============================================================
boxup("Experiment 2A", bandChar = "=") 

### Exp 2A data     
data <-read.csv("Exp2A_accumulated_rawdata.csv", header=T)

temp1<- subset(data, TimeOfEvent < 180)
temp2<- subset(data, TimeOfEvent == 180)

TargetPosVsDetection_exp2A<- ggplot() +
  # Layer for all targets that were found
  geom_point(data=temp1, aes(x=TargetPos_X, y= TargetPos_Y, fill= TimeOfEvent),shape = 21, size=2) + 
  scale_fill_viridis(name="Detection\ntime (s)", option="viridis",trans = "log10", breaks = c(1, 10, 100, 180), labels = c("1", "10", "100", "180")) +
  # Layer for all targets that were NOT found
  geom_point(data=temp2, aes(x=TargetPos_X, y= TargetPos_Y, color= "Target\n   not\n found"), shape = 20, size=3, show.legend = TRUE) +  # Ensure this layer shows in the legend
  # Add a custom legend entry for temp2 points
  scale_color_manual(values = c("Target\n   not\n found" = "red")) +  
  Graph.theme +
  facet_grid(Amplitude ~ Radial_frequency) + xlab("X target position (deg)") + ylab("Y target position (deg)") +
  scale_y_continuous(breaks= seq(-20,20,10), limits=c(-20,20), expand=c(0,0)) +
  scale_x_continuous(breaks= seq(-20,20,10), limits=c(-20,20), expand=c(0,0))
# coord_fixed(ratio=1) #This ensures that the x and Y axis are the same length
TargetPosVsDetection_exp2A


##===============================================================
##                        Experiment 2B                         =
##===============================================================
boxup("Experiment 2B", bandChar = "=") 

### Exp 2B data     
data <-read.csv("Exp2B_accumulated_rawdata.csv", header=T)

temp1<- subset(data, TimeOfEvent < 180)
temp2<- subset(data, TimeOfEvent == 180)

TargetPosVsDetection_exp2B<- ggplot() +
  # Layer for all targets that were found
  geom_point(data=temp1, aes(x=TargetPos_X, y= TargetPos_Y, fill= TimeOfEvent),shape = 21, size=2) + 
  scale_fill_viridis(name="Detection\ntime (s)", option="viridis",trans = "log10", breaks = c(1, 10, 100, 180), labels = c("1", "10", "100", "180")) +
  # Layer for all targets that were NOT found
  geom_point(data=temp2, aes(x=TargetPos_X, y= TargetPos_Y, color= "Target\n   not\n found"), shape = 20, size=3, show.legend = TRUE) +  # Ensure this layer shows in the legend
  # Add a custom legend entry for temp2 points
  scale_color_manual(values = c("Target\n   not\n found" = "red")) +  
  Graph.theme +
  facet_grid(Amplitude ~ Radial_frequency) + xlab("X target position (deg)") + ylab("Y target position (deg)") +
  scale_y_continuous(breaks= seq(-20,20,10), limits=c(-20,20), expand=c(0,0)) +
  scale_x_continuous(breaks= seq(-20,20,10), limits=c(-20,20), expand=c(0,0))
# coord_fixed(ratio=1) #This ensures that the x and Y axis are the same length
TargetPosVsDetection_exp2B


# Save plots
ggsave(file ="Target Pos vs detection time for Exp1A_v1_1.svg", device = 'svg', plot = TargetPosVsDetection_exp1A, width = 15, height = 11)
ggsave(file ="Target Pos vs detection time for Exp1B_v1_1.svg", device = 'svg', plot = TargetPosVsDetection_exp1B, width = 15, height = 11)
ggsave(file ="Target Pos vs detection time for Exp2A_v1_1.svg", device = 'svg', plot = TargetPosVsDetection_exp2A, width = 15, height = 11)
ggsave(file ="Target Pos vs detection time for Exp2B_v1_1.svg", device = 'svg', plot = TargetPosVsDetection_exp2B, width = 15, height = 11)
