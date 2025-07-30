# This is a copy of the R script used to generate supplementary figures S6-S9 from Smithers et al. (under review).

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
#   [1] LC_COLLATE=English_United States.utf8  LC_CTYPE=English_United States.utf8    LC_MONETARY=English_United States.utf8
# [4] LC_NUMERIC=C                           LC_TIME=English_United States.utf8    
# 
# time zone: America/New_York
# tzcode source: internal
# 
# attached base packages:
#   [1] grid      stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] png_0.1-8             bannerCommenter_1.0.0 viridis_0.6.4         viridisLite_0.4.2     ggplot2_3.4.4        
# 
# loaded via a namespace (and not attached):
#   [1] vctrs_0.6.5       cli_3.6.2         rlang_1.1.2       generics_0.1.3    labeling_0.4.3    glue_1.6.2       
# [7] colorspace_2.1-0  gridExtra_2.3     scales_1.3.0      fansi_1.0.6       munsell_0.5.0     tibble_3.2.1     
# [13] lifecycle_1.0.4   compiler_4.3.2    dplyr_1.1.4       pkgconfig_2.0.3   rstudioapi_0.15.0 farver_2.1.1     
# [19] R6_2.5.1          tidyselect_1.2.0  utf8_1.2.4        pillar_1.9.0      magrittr_2.0.3    tools_4.3.2      
# [25] withr_2.5.2       gtable_0.3.4 

library(ggplot2) #for graphs
library(viridis)
library(bannerCommenter) # For the section banners (optional)
library(grid)
library(png)

rm(list=ls(all=TRUE))

#Set up Working directory 
setwd("***Pathway to folder containing this R script and all '_accumulated mouse click data.csv' files***")

#Create graph theme
Graph.theme<-theme_bw()+
  theme(axis.text=element_text(size=20,colour="black"),
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


############################################################################
############################################################################
###                                                                      ###
###       MOUSE CLICK POSITION HEAT MAPS FOR EACH TARGET CONDITION       ###
###                                                                      ###
############################################################################
############################################################################
banner("Mouse click position heat maps for each target condition", emph = TRUE)
# The graphs show the position of each click and the colour of each point shows the number of neighbouring points within a 0.25deg radius of it.

###### Function to count the number of points within a set radius around each point #######
# The function is based off the function used by "geom_pointdensity" from here: https://github.com/LKremer/ggpointdensity/tree/master
count_neighbors <- function(x, y, r2) { # Note that this assumes that the x and y axis are on the same scale. 
  sapply(1:length(x), function(i) {
    sum(((x[i] - x) ^ 2) + ((y[i] - y) ^ 2) < r2)
  })
}
r2 = 0.25 # Radius around each point within which to count the number of neighboring points
PlotRange = 4 # Range (in deg) over which to plots clicks either side of the center of the target 


##===============================================================
##                        Experiment 1A                         =
##===============================================================
boxup("Experiment 1A", bandChar = "=")  

Exp1A_data <-read.csv("Exp1A_accumulated mouse click data.csv", header=T)
Exp1A_data["HITorMISS"] <- Exp1A_data["TargetClicked"]
Exp1A_data$HITorMISS[Exp1A_data$HITorMISS != "miss"] <- "hit"
N <- nrow(unique(Exp1A_data["Subject_ID"]))
print(N) # N=23

ConditionList <- unique(Exp1A_data["Condition_name"])
ClickDensityData <- data.frame()

# Calculate the density of clicks for each target condition
for (stim in 1:nrow(ConditionList)){
  temp<- subset(Exp1A_data, Condition_name== ConditionList[stim,1])
  temp$density <- count_neighbors(temp$RelativeClickPos_X, temp$RelativeClickPos_Y, r2 = r2)
  ClickDensityData <- rbind(ClickDensityData,temp)
}
  
Subdata<- subset(ClickDensityData, RelativeClickPos_X < PlotRange & RelativeClickPos_X > -PlotRange & RelativeClickPos_Y < PlotRange & RelativeClickPos_Y > -PlotRange) # Make it so that we only have the relative position of the click that was counted as the target being detected 

# Plot graph
Exp1A_ClickPos_ScatterWDensity <- ggplot(Subdata, aes(x=RelativeClickPos_X, y= RelativeClickPos_Y, color=density, shape=HITorMISS))+#,alpha=log10(density))) + 
  geom_point(size=2) + scale_color_viridis(option="viridis",trans = "log10") + Graph.theme +
  scale_shape_manual(values = c("hit" = 19, "miss" = 4)) +
  facet_grid(Amplitude ~ Radial_frequency) + xlab("X click position relative to target center (deg)") + ylab("Y click position relative to target center (deg)") +
  coord_fixed(ratio=1) #This ensures that the x and Y axis are the same length

Exp1A_ClickPos_ScatterWDensity 


##===============================================================
##                        Experiment 1B                         =
##===============================================================
boxup("Experiment 1B", bandChar = "=")     

Exp1B_data <-read.csv("Exp1B_accumulated mouse click data.csv", header=T)
Exp1B_data["HITorMISS"] <- Exp1B_data["TargetClicked"]
Exp1B_data$HITorMISS[Exp1B_data$HITorMISS != "miss"] <- "hit"
N <- nrow(unique(Exp1B_data["Subject_ID"]))
print(N) # N=20

ConditionList <- unique(Exp1B_data["Condition_name"])
ClickDensityData <- data.frame()

# Calculate the density of clicks for each target condition
for (stim in 1:nrow(ConditionList)){
  temp<- subset(Exp1B_data, Condition_name== ConditionList[stim,1])
  temp$density <- count_neighbors(temp$RelativeClickPos_X, temp$RelativeClickPos_Y, r2 = r2)
  ClickDensityData <- rbind(ClickDensityData,temp)
}

Subdata<- subset(ClickDensityData, RelativeClickPos_X < PlotRange & RelativeClickPos_X > -PlotRange & RelativeClickPos_Y < PlotRange & RelativeClickPos_Y > -PlotRange) # Make it so that we only have the relative position of the click that was counted as the target being detected 

# Plot graphs
Exp1B_ClickPos_ScatterWDensity <- ggplot(Subdata, aes(x=RelativeClickPos_X, y= RelativeClickPos_Y, color=density, shape=HITorMISS))+#,alpha=log10(density))) + 
  geom_point(size=2) + scale_color_viridis(option="viridis",trans = "log10") + Graph.theme +
  scale_shape_manual(values = c("hit" = 19, "miss" = 4)) +
  facet_grid(Amplitude ~ Radial_frequency) + xlab("X click position relative to target center (deg)") + ylab("Y click position relative to target center (deg)") +
  coord_fixed(ratio=1) 

Exp1B_ClickPos_ScatterWDensity 


##===============================================================
##                        Experiment 2A                         =
##===============================================================
boxup("Experiment 2A", bandChar = "=")     

Exp2A_data <-read.csv("Exp2A_accumulated mouse click data.csv", header=T)
Exp2A_data["HITorMISS"] <- Exp2A_data["TargetClicked"]
Exp2A_data$HITorMISS[Exp2A_data$HITorMISS != "miss"] <- "hit"
N <- nrow(unique(Exp2A_data["Subject_ID"]))
print(N) # N=26

ConditionList <- unique(Exp2A_data["Condition_name"])
ClickDensityData <- data.frame()

# Calculate the density of clicks for each target condition
for (stim in 1:nrow(ConditionList)){
  temp<- subset(Exp2A_data, Condition_name== ConditionList[stim,1])
  temp$density <- count_neighbors(temp$RelativeClickPos_X, temp$RelativeClickPos_Y, r2 = r2)
  ClickDensityData <- rbind(ClickDensityData,temp)
}

Subdata<- subset(ClickDensityData, RelativeClickPos_X < PlotRange & RelativeClickPos_X > -PlotRange & RelativeClickPos_Y < PlotRange & RelativeClickPos_Y > -PlotRange) # Make it so that we only have the relative position of the click that was counted as the target being detected 

# Plot graphs
Exp2A_ClickPos_ScatterWDensity <- ggplot(Subdata, aes(x=RelativeClickPos_X, y= RelativeClickPos_Y, color=density, shape=HITorMISS))+#,alpha=log10(density))) + 
  geom_point(size=2) + scale_color_viridis(option="viridis",trans = "log10") + Graph.theme +
  scale_shape_manual(values = c("hit" = 19, "miss" = 4)) +
  facet_grid(Amplitude ~ Radial_frequency) + xlab("X click position relative to target center (deg)") + ylab("Y click position relative to target center (deg)") +
  coord_fixed(ratio=1)

Exp2A_ClickPos_ScatterWDensity 


##===============================================================
##                        Experiment 2B                         =
##===============================================================
boxup("Experiment 2B", bandChar = "=")    

Exp2B_data <-read.csv("Exp2B_accumulated mouse click data.csv", header=T)
Exp2B_data["HITorMISS"] <- Exp2B_data["TargetClicked"]
Exp2B_data$HITorMISS[Exp2B_data$HITorMISS != "miss"] <- "hit"
N <- nrow(unique(Exp2B_data["Subject_ID"]))
print(N) # N=20

ConditionList <- unique(Exp2B_data["Condition_name"])
ClickDensityData <- data.frame()

# Calculate the density of clicks for each target condition
for (stim in 1:nrow(ConditionList)){
  temp<- subset(Exp2B_data, Condition_name== ConditionList[stim,1])
  temp$density <- count_neighbors(temp$RelativeClickPos_X, temp$RelativeClickPos_Y, r2 = r2)
  ClickDensityData <- rbind(ClickDensityData,temp)
}

Subdata<- subset(ClickDensityData, RelativeClickPos_X < PlotRange & RelativeClickPos_X > -PlotRange & RelativeClickPos_Y < PlotRange & RelativeClickPos_Y > -PlotRange) # Make it so that we only have the relative position of the click that was counted as the target being detected 

# plot graph
Exp2B_ClickPos_ScatterWDensity <- ggplot(Subdata, aes(x=RelativeClickPos_X, y= RelativeClickPos_Y, color=density, shape=HITorMISS))+#,alpha=log10(density))) + 
  geom_point(size=2) + scale_color_viridis(option="viridis",trans = "log10") + Graph.theme +
  scale_shape_manual(values = c("hit" = 19, "miss" = 4)) +
  facet_grid(Amplitude ~ Radial_frequency) + xlab("X click position relative to target center (deg)") + ylab("Y click position relative to target center (deg)") +
  coord_fixed(ratio=1) 

Exp2B_ClickPos_ScatterWDensity 

# Save plots
ggsave(file ="Click Pos density plot for Exp1A.svg", device = 'svg', plot = Exp1A_ClickPos_ScatterWDensity, width = 14, height = 14)
ggsave(file ="Click Pos density plot for Exp1B.svg", device = 'svg', plot = Exp1B_ClickPos_ScatterWDensity, width = 14, height = 14)
ggsave(file ="Click Pos density plot for Exp2A.svg", device = 'svg', plot = Exp2A_ClickPos_ScatterWDensity, width = 14, height = 14)
ggsave(file ="Click Pos density plot for Exp2B.svg", device = 'svg', plot = Exp2B_ClickPos_ScatterWDensity, width = 14, height = 14)


