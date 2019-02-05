## Code to run the analyses and create figures for Little, Fronhofer, & Altermatt 2019
## "Dispersal syndromes can impact ecosystem functioning in spatially structured freshwater populations"
## Date: 05 February 2019

#### load packages and functions ####

# load packages 

library(knitr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(lme4)
library(lmerTest)
library(multcomp)
library(arm)
library(gridExtra)
library(MuMIn)

# call this summary funtion

summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  require(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

#### load paths ####

# data.path <- PATH TO DOWNLOADED DATA from Dryad

#### data cleanup and formatting for G. fossarum ####

# load data

GF_disp_weights_df <- read.csv(paste(data.path, 
                                     "/Suppdata1_G_fossarum_dispnet_bodyweights.csv",
                                     sep=""), skip=3, sep=";", dec=",")

GF_consumption_data <- read.csv(paste(data.path, 
                                      "/Suppdata2_GF_weight_and_survival.csv",
                                      sep=""), skip=3, sep=",", dec=",")

# correct a typo
GF_disp_weights_df[which(GF_disp_weights_df$Res_or_Dips=="Disperer"),2] <-
  c("Disperser")

# order the resident and disperser levels
GF_disp_weights_df$Res_or_Disp <- factor(GF_disp_weights_df$Res_or_Disp, 
                                         levels=c("Residents", "Dispersers"),
                                         ordered=TRUE)

# extract just the relevant columns
GF_disp_weights_df <- GF_disp_weights_df[,-1]
colnames(GF_disp_weights_df) <- c("Mesocosm_ID", "Res_or_Dips", "Treatment",
                                  "Mass_mg", "Note")  

# take only the measurements after the feeding experiment (see Dryad README)
GF_disp_weights_df <- 
  GF_disp_weights_df[which(GF_disp_weights_df$Note=="post-feeding-experiment"),]
GF_disp_weights_df <- GF_disp_weights_df[,-5]

# convert g to mg
str(GF_disp_weights_df)
GF_disp_weights_df$Mass_mg <- as.numeric(as.character(GF_disp_weights_df$Mass_mg))
GF_disp_weights_df$Mass_mg <- GF_disp_weights_df$Mass_mg*1000

#### calculate the "weight days" for G. fossarum mesocosms ####

# calculate average weight of surviving amphipods for each mesocosm/replicate

GF_dispnet_average_weights <- as.data.frame(GF_disp_weights_df %>% 
                                              group_by(Mesocosm_ID, Res_or_Dips, Treatment) %>% 
                                              dplyr::summarise(mean_weight = mean(Mass_mg),
                                                               n_weighed = n()))

# make a version that removes res/disp codes
# keep treatment codes though! 

GF_all_averages <- GF_dispnet_average_weights[,-2]

# calculate global mean and sd weight 

GF_meanweight <- mean(GF_disp_weights_df$Mass_mg)
GF_sdweight <- sd(GF_disp_weights_df$Mass_mg)

# join the consumption data with the weights from density experiments
GF_consumption_data_joined = 
  merge(GF_consumption_data, GF_all_averages,
        by.y="Mesocosm_ID", by.x="Mesocosm", all.x=TRUE)

# order the residents and dispersers
GF_consumption_data_joined$Treatment.x <-
  factor(GF_consumption_data_joined$Treatment.x,
         levels=c("Residents", "Dispersers"), ordered=TRUE)

# rename the treatment columns
names(GF_consumption_data_joined)[names(GF_consumption_data_joined) ==
                                    'Treatment.x'] <- 'ResDisp'
names(GF_consumption_data_joined)[names(GF_consumption_data_joined) ==
                                    'Treatment.y'] <- 'Condition'

# data reformat
GF_consumption_data_joined$n_weighed <- 
  as.numeric(as.character(GF_consumption_data_joined$n_weighed))
GF_consumption_data_joined$Amphipod_Days <- 
  as.numeric(as.character(GF_consumption_data_joined$Amphipod_Days))

# calculate biomass-days for weighed amphipods
GF_consumption_data_joined$weight_days_measured <- 
  GF_consumption_data_joined$mean_weight*19

# calculate biomass-days for non-surviving amphipods based on global average
GF_consumption_data_joined$missing_amphi_days <- 
  GF_consumption_data_joined$Amphipod_Days -
  (GF_consumption_data_joined$n_weighed*19)

GF_consumption_data_joined$weight_days_estimated <- 
  GF_consumption_data_joined$missing_amphi_days * GF_meanweight

# for mesocosms w/ "missing" amphipod-days, sum the two to get total biomass-days,
# or else just multiply the mean weight by the number of amphipod-days
GF_consumption_data_joined$weight_days_total <- 
  ifelse(GF_consumption_data_joined$missing_amphi_days > 0,
         rowSums(GF_consumption_data_joined[,c("weight_days_estimated",
                                               "weight_days_measured")], na.rm=TRUE),
         GF_consumption_data_joined$Amphipod_Days*GF_consumption_data_joined$mean_weight)

# need to also do this for ones with NO amphipods measured (NA)
GF_consumption_data_joined[is.na(GF_consumption_data_joined$n_weighed),]$weight_days_total <- 
  GF_consumption_data_joined[is.na(GF_consumption_data_joined$n_weighed),]$Amphipod_Days*GF_meanweight

#### adjust mass loss by bodyweight for G. fossarum ####

GF_consumption_data_joined$Mass_lost_g <- 
  as.numeric(as.character(GF_consumption_data_joined$Mass_lost_g))

GF_consumption_data_joined$Mass_loss_per_biomass_day_mg <- 
  (GF_consumption_data_joined$Mass_lost_g*1000)/GF_consumption_data_joined$weight_days_total

# replace negative values with zero
GF_consumption_data_joined$Mass_loss_per_biomass_day_mg <-
  ifelse(GF_consumption_data_joined$Mass_loss_per_biomass_day_mg < 0, 
         0,
         GF_consumption_data_joined$Mass_loss_per_biomass_day_mg)

#### load and cleanup D. villosus data ####

DV_disp_weights_df <- read.csv(paste(data.path,
                                     "/Suppdata4_D_villosus_dispnet_bodyweights.csv",
                                     sep=""), skip=3, sep=";", dec=",")

DV_consumption_data <- read.csv(paste(data.path, 
                                      "/Suppdata3_DV_weight_and_survival.csv",
                                      sep=""), skip=3, sep=",", dec=",")

# order the resident and disperser levels
DV_disp_weights_df$Res_or_Disp <- factor(DV_disp_weights_df$Res_or_Disp, 
                                         levels=c("Residents", "Dispersers"),
                                         ordered=TRUE)

# extract relevant parts of dataframe
DV_disp_weights_df <- DV_disp_weights_df[,-1]
colnames(DV_disp_weights_df) <- c("Mesocosm_ID", "Res_or_Dips", "Treatment",
                                  "Mass_mg", "Note")  

# correct a typo
DV_disp_weights_df[which(DV_disp_weights_df$Res_or_Dips=="Disperer"),2] <-
  c("Disperser")

# take only the measurements after the feeding experiment
DV_disp_weights_df <- 
  DV_disp_weights_df[which(DV_disp_weights_df$Note=="post-feeding-experiment"),]
DV_disp_weights_df <- DV_disp_weights_df[,-5]

# convert g to mg
str(DV_disp_weights_df)
DV_disp_weights_df$Mass_mg <- as.numeric(as.character(DV_disp_weights_df$Mass_mg))
DV_disp_weights_df$Mass_mg <- DV_disp_weights_df$Mass_mg*1000

#### calculate the "weight days" for D. villosus mesocosms ####

# calculate average weight of surviving amphipods for each mesocosm/replicate

DV_dispnet_average_weights <- as.data.frame(DV_disp_weights_df %>% 
                                              group_by(Mesocosm_ID, Res_or_Dips, Treatment) %>% 
                                              dplyr::summarise(mean_weight = mean(Mass_mg),
                                                               n_weighed = n()))

# make a version that removes treatment codes
DV_all_averages <- DV_dispnet_average_weights[,-2]

# calculate global mean and sd weight for D. villosus
DV_meanweight <- mean(DV_disp_weights_df$Mass_mg)
DV_sdweight <- sd(DV_disp_weights_df$Mass_mg)

# join the consumption data with the weights from density experiments
DV_consumption_data_joined = 
  merge(DV_consumption_data, DV_all_averages,
        by.y="Mesocosm_ID", by.x="Mesocosm", all.x=TRUE)

# order the residents and dispersers
DV_consumption_data_joined$Treatment.x <-
  factor(DV_consumption_data_joined$Treatment.x,
         levels=c("Residents", "Dispersers"), ordered=TRUE)

# rename the treatment columns
names(DV_consumption_data_joined)[names(DV_consumption_data_joined) ==
                                    'Treatment.x'] <- 'ResDisp'
names(DV_consumption_data_joined)[names(DV_consumption_data_joined) ==
                                    'Treatment.y'] <- 'Condition'

# data reformat
DV_consumption_data_joined$n_weighed <- 
  as.numeric(as.character(DV_consumption_data_joined$n_weighed))
DV_consumption_data_joined$Amphipod_Days <- 
  as.numeric(as.character(DV_consumption_data_joined$Amphipod_Days))

# calculate biomass-days for weighed amphipods
DV_consumption_data_joined$weight_days_measured <- 
  DV_consumption_data_joined$mean_weight*19

# calculate biomass-days for non-surviving amphipods based on global average
DV_consumption_data_joined$missing_amphi_days <- 
  DV_consumption_data_joined$Amphipod_Days -
  (DV_consumption_data_joined$n_weighed*19)

DV_consumption_data_joined$weight_days_estimated <- 
  DV_consumption_data_joined$missing_amphi_days * DV_meanweight

# for mesocosms w/ "missing" amphipod-days, sum the two to get total biomass-days,
# or else just multiply the mean weight by the number of amphipod-days
DV_consumption_data_joined$weight_days_total <- 
  ifelse(DV_consumption_data_joined$missing_amphi_days > 0,
         rowSums(DV_consumption_data_joined[,c("weight_days_estimated",
                                               "weight_days_measured")], na.rm=TRUE),
         DV_consumption_data_joined$Amphipod_Days*DV_consumption_data_joined$mean_weight)

# need to also do this for ones with NO amphipods measured (NA)
DV_consumption_data_joined[is.na(DV_consumption_data_joined$n_weighed),]$weight_days_total <- 
  DV_consumption_data_joined[is.na(DV_consumption_data_joined$n_weighed),]$Amphipod_Days*DV_meanweight

#### adjust the mass loss by bodyweight days for D. villosus ####

DV_consumption_data_joined$Mass_lost_g <- 
  as.numeric(as.character(DV_consumption_data_joined$Mass_lost_g))

DV_consumption_data_joined$Mass_loss_per_biomass_day_mg <- 
  (DV_consumption_data_joined$Mass_lost_g*1000)/DV_consumption_data_joined$weight_days_total

#### add some data in manually: set up replicate blocks ####

# need to pair the resident and disperser mesocosms from each two-patch mesocosm setup

GF_consumption_data_joined$replicate_block <- NA

GF_consumption_data_joined$replicate_block[which(GF_consumption_data_joined$Mesocosm==6|
                                                   GF_consumption_data_joined$Mesocosm==2)] <- 
  1
GF_consumption_data_joined$replicate_block[which(GF_consumption_data_joined$Mesocosm==9|
                                                   GF_consumption_data_joined$Mesocosm==8)] <- 
  2
GF_consumption_data_joined$replicate_block[which(GF_consumption_data_joined$Mesocosm==10|
                                                   GF_consumption_data_joined$Mesocosm==13)] <- 
  3
GF_consumption_data_joined$replicate_block[which(GF_consumption_data_joined$Mesocosm==14|
                                                   GF_consumption_data_joined$Mesocosm==21)] <- 
  4
GF_consumption_data_joined$replicate_block[which(GF_consumption_data_joined$Mesocosm==16|
                                                   GF_consumption_data_joined$Mesocosm==25)] <- 
  5
GF_consumption_data_joined$replicate_block[which(GF_consumption_data_joined$Mesocosm==18|
                                                   GF_consumption_data_joined$Mesocosm==26)] <- 
  6
GF_consumption_data_joined$replicate_block[which(GF_consumption_data_joined$Mesocosm==33|
                                                   GF_consumption_data_joined$Mesocosm==28)] <- 
  7
GF_consumption_data_joined$replicate_block[which(GF_consumption_data_joined$Mesocosm==34|
                                                   GF_consumption_data_joined$Mesocosm==29)] <- 
  8
GF_consumption_data_joined$replicate_block[which(GF_consumption_data_joined$Mesocosm==38|
                                                   GF_consumption_data_joined$Mesocosm==30)] <- 
  9
GF_consumption_data_joined$replicate_block[which(GF_consumption_data_joined$Mesocosm==43|
                                                   GF_consumption_data_joined$Mesocosm==31)] <- 
  10
GF_consumption_data_joined$replicate_block[which(GF_consumption_data_joined$Mesocosm==47|
                                                   GF_consumption_data_joined$Mesocosm==32)] <- 
  11
GF_consumption_data_joined$replicate_block[which(GF_consumption_data_joined$Mesocosm==48|
                                                   GF_consumption_data_joined$Mesocosm==41)] <- 
  12
GF_consumption_data_joined$replicate_block[which(GF_consumption_data_joined$Mesocosm==55|
                                                   GF_consumption_data_joined$Mesocosm==53)] <- 
  13
GF_consumption_data_joined$replicate_block[which(GF_consumption_data_joined$Mesocosm==58|
                                                   GF_consumption_data_joined$Mesocosm==56)] <- 
  14
GF_consumption_data_joined$replicate_block[which(GF_consumption_data_joined$Mesocosm==60|
                                                   GF_consumption_data_joined$Mesocosm==65)] <- 
  15
GF_consumption_data_joined$replicate_block[which(GF_consumption_data_joined$Mesocosm==61|
                                                   GF_consumption_data_joined$Mesocosm==71)] <- 
  16
GF_consumption_data_joined$replicate_block[which(GF_consumption_data_joined$Mesocosm==64|
                                                   GF_consumption_data_joined$Mesocosm==73)] <- 
  17
GF_consumption_data_joined$replicate_block[which(GF_consumption_data_joined$Mesocosm==68|
                                                   GF_consumption_data_joined$Mesocosm==78)] <- 
  18
GF_consumption_data_joined$replicate_block[which(GF_consumption_data_joined$Mesocosm==82|
                                                   GF_consumption_data_joined$Mesocosm==81)] <- 
  19
GF_consumption_data_joined$replicate_block[which(GF_consumption_data_joined$Mesocosm==91|
                                                   GF_consumption_data_joined$Mesocosm==86)] <- 
  20
GF_consumption_data_joined$replicate_block[which(GF_consumption_data_joined$Mesocosm==94|
                                                   GF_consumption_data_joined$Mesocosm==90)] <- 
  21
GF_consumption_data_joined$replicate_block[which(GF_consumption_data_joined$Mesocosm==98|
                                                   GF_consumption_data_joined$Mesocosm==97)] <- 
  22
GF_consumption_data_joined$replicate_block[which(GF_consumption_data_joined$Mesocosm==103|
                                                   GF_consumption_data_joined$Mesocosm==101)] <- 
  23
GF_consumption_data_joined$replicate_block[which(GF_consumption_data_joined$Mesocosm==105|
                                                   GF_consumption_data_joined$Mesocosm==106)] <- 
  24
GF_consumption_data_joined$replicate_block[which(GF_consumption_data_joined$Mesocosm==110|
                                                   GF_consumption_data_joined$Mesocosm==116)] <- 
  25
GF_consumption_data_joined$replicate_block[which(GF_consumption_data_joined$Mesocosm==111|
                                                   GF_consumption_data_joined$Mesocosm==121)] <- 
  26
GF_consumption_data_joined$replicate_block[which(GF_consumption_data_joined$Mesocosm==112|
                                                   GF_consumption_data_joined$Mesocosm==124)] <- 
  27
GF_consumption_data_joined$replicate_block[which(GF_consumption_data_joined$Mesocosm==117|
                                                   GF_consumption_data_joined$Mesocosm==127)] <- 
  28
GF_consumption_data_joined$replicate_block[which(GF_consumption_data_joined$Mesocosm==125|
                                                   GF_consumption_data_joined$Mesocosm==129)] <- 
  29
GF_consumption_data_joined$replicate_block[which(GF_consumption_data_joined$Mesocosm==135|
                                                   GF_consumption_data_joined$Mesocosm==138)] <- 
  30
GF_consumption_data_joined$replicate_block[which(GF_consumption_data_joined$Mesocosm==137|
                                                   GF_consumption_data_joined$Mesocosm==139)] <- 
  31
GF_consumption_data_joined$replicate_block[which(GF_consumption_data_joined$Mesocosm==148|
                                                   GF_consumption_data_joined$Mesocosm==140)] <- 
  32
GF_consumption_data_joined$replicate_block[which(GF_consumption_data_joined$Mesocosm==152|
                                                   GF_consumption_data_joined$Mesocosm==142)] <- 
  33
GF_consumption_data_joined$replicate_block[which(GF_consumption_data_joined$Mesocosm==153|
                                                   GF_consumption_data_joined$Mesocosm==143)] <- 
  34
GF_consumption_data_joined$replicate_block[which(GF_consumption_data_joined$Mesocosm==155|
                                                   GF_consumption_data_joined$Mesocosm==145)] <- 
  35
GF_consumption_data_joined$replicate_block[which(GF_consumption_data_joined$Mesocosm==160|
                                                   GF_consumption_data_joined$Mesocosm==168)] <- 
  36
GF_consumption_data_joined$replicate_block[which(GF_consumption_data_joined$Mesocosm==165|
                                                   GF_consumption_data_joined$Mesocosm==173)] <- 
  37
GF_consumption_data_joined$replicate_block[which(GF_consumption_data_joined$Mesocosm==174|
                                                   GF_consumption_data_joined$Mesocosm==176)] <- 
  38
GF_consumption_data_joined$replicate_block[which(GF_consumption_data_joined$Mesocosm==175|
                                                   GF_consumption_data_joined$Mesocosm==177)] <- 
  39
GF_consumption_data_joined$replicate_block[which(GF_consumption_data_joined$Mesocosm==182|
                                                   GF_consumption_data_joined$Mesocosm==178)] <- 
  40

DV_consumption_data_joined$replicate_block <- NA

DV_consumption_data_joined$replicate_block[which(DV_consumption_data_joined$Mesocosm==6|
                                                   DV_consumption_data_joined$Mesocosm==2)] <- 
  1
DV_consumption_data_joined$replicate_block[which(DV_consumption_data_joined$Mesocosm==9|
                                                   DV_consumption_data_joined$Mesocosm==8)] <- 
  2
DV_consumption_data_joined$replicate_block[which(DV_consumption_data_joined$Mesocosm==10|
                                                   DV_consumption_data_joined$Mesocosm==13)] <- 
  3
DV_consumption_data_joined$replicate_block[which(DV_consumption_data_joined$Mesocosm==14|
                                                   DV_consumption_data_joined$Mesocosm==21)] <- 
  4
DV_consumption_data_joined$replicate_block[which(DV_consumption_data_joined$Mesocosm==16|
                                                   DV_consumption_data_joined$Mesocosm==25)] <- 
  5
DV_consumption_data_joined$replicate_block[which(DV_consumption_data_joined$Mesocosm==18|
                                                   DV_consumption_data_joined$Mesocosm==26)] <- 
  6
DV_consumption_data_joined$replicate_block[which(DV_consumption_data_joined$Mesocosm==33|
                                                   DV_consumption_data_joined$Mesocosm==28)] <- 
  7
DV_consumption_data_joined$replicate_block[which(DV_consumption_data_joined$Mesocosm==34|
                                                   DV_consumption_data_joined$Mesocosm==29)] <- 
  8
DV_consumption_data_joined$replicate_block[which(DV_consumption_data_joined$Mesocosm==38|
                                                   DV_consumption_data_joined$Mesocosm==30)] <- 
  9
DV_consumption_data_joined$replicate_block[which(DV_consumption_data_joined$Mesocosm==43|
                                                   DV_consumption_data_joined$Mesocosm==31)] <- 
  10
DV_consumption_data_joined$replicate_block[which(DV_consumption_data_joined$Mesocosm==47|
                                                   DV_consumption_data_joined$Mesocosm==32)] <- 
  11
DV_consumption_data_joined$replicate_block[which(DV_consumption_data_joined$Mesocosm==48|
                                                   DV_consumption_data_joined$Mesocosm==41)] <- 
  12
DV_consumption_data_joined$replicate_block[which(DV_consumption_data_joined$Mesocosm==55|
                                                   DV_consumption_data_joined$Mesocosm==53)] <- 
  13
DV_consumption_data_joined$replicate_block[which(DV_consumption_data_joined$Mesocosm==58|
                                                   DV_consumption_data_joined$Mesocosm==56)] <- 
  14
DV_consumption_data_joined$replicate_block[which(DV_consumption_data_joined$Mesocosm==60|
                                                   DV_consumption_data_joined$Mesocosm==65)] <- 
  15
DV_consumption_data_joined$replicate_block[which(DV_consumption_data_joined$Mesocosm==61|
                                                   DV_consumption_data_joined$Mesocosm==71)] <- 
  16
DV_consumption_data_joined$replicate_block[which(DV_consumption_data_joined$Mesocosm==64|
                                                   DV_consumption_data_joined$Mesocosm==73)] <- 
  17
DV_consumption_data_joined$replicate_block[which(DV_consumption_data_joined$Mesocosm==68|
                                                   DV_consumption_data_joined$Mesocosm==78)] <- 
  18
DV_consumption_data_joined$replicate_block[which(DV_consumption_data_joined$Mesocosm==82|
                                                   DV_consumption_data_joined$Mesocosm==81)] <- 
  19
DV_consumption_data_joined$replicate_block[which(DV_consumption_data_joined$Mesocosm==91|
                                                   DV_consumption_data_joined$Mesocosm==86)] <- 
  20
DV_consumption_data_joined$replicate_block[which(DV_consumption_data_joined$Mesocosm==94|
                                                   DV_consumption_data_joined$Mesocosm==90)] <- 
  21
DV_consumption_data_joined$replicate_block[which(DV_consumption_data_joined$Mesocosm==98|
                                                   DV_consumption_data_joined$Mesocosm==97)] <- 
  22
DV_consumption_data_joined$replicate_block[which(DV_consumption_data_joined$Mesocosm==103|
                                                   DV_consumption_data_joined$Mesocosm==101)] <- 
  23
DV_consumption_data_joined$replicate_block[which(DV_consumption_data_joined$Mesocosm==105|
                                                   DV_consumption_data_joined$Mesocosm==106)] <- 
  24
DV_consumption_data_joined$replicate_block[which(DV_consumption_data_joined$Mesocosm==110|
                                                   DV_consumption_data_joined$Mesocosm==116)] <- 
  25
DV_consumption_data_joined$replicate_block[which(DV_consumption_data_joined$Mesocosm==111|
                                                   DV_consumption_data_joined$Mesocosm==121)] <- 
  26
DV_consumption_data_joined$replicate_block[which(DV_consumption_data_joined$Mesocosm==112|
                                                   DV_consumption_data_joined$Mesocosm==124)] <- 
  27
DV_consumption_data_joined$replicate_block[which(DV_consumption_data_joined$Mesocosm==117|
                                                   DV_consumption_data_joined$Mesocosm==127)] <- 
  28
DV_consumption_data_joined$replicate_block[which(DV_consumption_data_joined$Mesocosm==125|
                                                   DV_consumption_data_joined$Mesocosm==129)] <- 
  29
DV_consumption_data_joined$replicate_block[which(DV_consumption_data_joined$Mesocosm==135|
                                                   DV_consumption_data_joined$Mesocosm==138)] <- 
  30
DV_consumption_data_joined$replicate_block[which(DV_consumption_data_joined$Mesocosm==137|
                                                   DV_consumption_data_joined$Mesocosm==139)] <- 
  31
DV_consumption_data_joined$replicate_block[which(DV_consumption_data_joined$Mesocosm==148|
                                                   DV_consumption_data_joined$Mesocosm==140)] <- 
  32
DV_consumption_data_joined$replicate_block[which(DV_consumption_data_joined$Mesocosm==152|
                                                   DV_consumption_data_joined$Mesocosm==142)] <- 
  33
DV_consumption_data_joined$replicate_block[which(DV_consumption_data_joined$Mesocosm==153|
                                                   DV_consumption_data_joined$Mesocosm==143)] <- 
  34
DV_consumption_data_joined$replicate_block[which(DV_consumption_data_joined$Mesocosm==155|
                                                   DV_consumption_data_joined$Mesocosm==145)] <- 
  35
DV_consumption_data_joined$replicate_block[which(DV_consumption_data_joined$Mesocosm==160|
                                                   DV_consumption_data_joined$Mesocosm==168)] <- 
  36
DV_consumption_data_joined$replicate_block[which(DV_consumption_data_joined$Mesocosm==165|
                                                   DV_consumption_data_joined$Mesocosm==173)] <- 
  37
DV_consumption_data_joined$replicate_block[which(DV_consumption_data_joined$Mesocosm==174|
                                                   DV_consumption_data_joined$Mesocosm==176)] <- 
  38
DV_consumption_data_joined$replicate_block[which(DV_consumption_data_joined$Mesocosm==175|
                                                   DV_consumption_data_joined$Mesocosm==177)] <- 
  39
DV_consumption_data_joined$replicate_block[which(DV_consumption_data_joined$Mesocosm==182|
                                                   DV_consumption_data_joined$Mesocosm==178)] <- 
  40

#### GF model specification and summary ####

# random intercept for each replicate block
# dispersal status and density as fixed factors
GFdispnetlmer_randomintercept <- 
  lmer(Mass_loss_per_biomass_day_mg ~ ResDisp + Amphipod_Days +
         (1|replicate_block), data=GF_consumption_data_joined)
shapiro.test(resid(GFdispnetlmer_randomintercept))

# because the residuals are really not normal, try square-root transforming
GFdispnetlmer_randomintercept_sqrt <- 
  lmer(sqrt(Mass_loss_per_biomass_day_mg) ~ ResDisp + Amphipod_Days +
         (1|replicate_block), data=GF_consumption_data_joined)
shapiro.test(resid(GFdispnetlmer_randomintercept_sqrt))
# this is not normal, but it is better
plot(GFdispnetlmer_randomintercept_sqrt)
qqnorm(resid(GFdispnetlmer_randomintercept_sqrt))
hist(resid(GFdispnetlmer_randomintercept_sqrt))

# summary and results from model
summary(GFdispnetlmer_randomintercept_sqrt)

r.squaredGLMM(GFdispnetlmer_randomintercept_sqrt)

summary(glht(GFdispnetlmer_randomintercept_sqrt, 
             linfct = mcp(ResDisp = "Tukey")))

#### DV model specification and summary ####

# random intercept for each replicate block
# dispersal status and density as fixed factors
DVdispnetlmer_randomintercept <- 
  lmer(Mass_loss_per_biomass_day_mg ~ ResDisp + Amphipod_Days +
         (1|replicate_block), data=DV_consumption_data_joined)

# square root transforming does improve the resiuals, etc., as in GF
DVdispnetlmer_randomintercept_sqrt <- 
  lmer(sqrt(Mass_loss_per_biomass_day_mg) ~ ResDisp + Amphipod_Days +
         (1|replicate_block), data=DV_consumption_data_joined)
shapiro.test(resid(DVdispnetlmer_randomintercept_sqrt))
plot(DVdispnetlmer_randomintercept_sqrt)
qqnorm(resid(DVdispnetlmer_randomintercept_sqrt))
hist(resid(DVdispnetlmer_randomintercept_sqrt))

# model results and summary
summary(DVdispnetlmer_randomintercept_sqrt)

r.squaredGLMM(DVdispnetlmer_randomintercept_sqrt)

summary(glht(DVdispnetlmer_randomintercept_sqrt, 
             linfct = mcp(ResDisp = "Tukey")))

#### Figure 1 ####

sumGF <-
  summarySE(data=GF_consumption_data_joined, 
            measurevar = c("Mass_loss_per_biomass_day_mg"), 
            groupvars=c("ResDisp"))
sumGF$ResDisp <- factor(sumGF$ResDisp,
                        levels=c("Residents", "Dispersers"),
                        ordered=TRUE)

# colors
dispcols <- c("Dispersers" = "#b2182b", "Residents" = "#67a9cf")

(GFdispnetplot <- 
    ggplot(sumGF, 
           aes(x = ResDisp, y = Mass_loss_per_biomass_day_mg,
               fill = ResDisp)) + 
    geom_bar(stat='identity', colour="black") +
    geom_point(data = GF_consumption_data_joined, 
               mapping = aes(x = ResDisp, 
                             y = Mass_loss_per_biomass_day_mg),
               fill="darkgray", alpha=0.5, size=3)+
    ylab(bquote('Consumption (mg mg'^-1*' amphipod day'^-1*')')) + 
    ylim(c(0,15)) + ggtitle("Gammarus fossarum")+
    scale_fill_manual("Dispersal Type", values=dispcols)+
    theme_classic() +
    theme(axis.line.x = element_line(colour = 'black', 
                                     size=0.5, linetype='solid'),
          axis.line.y = element_line(colour = 'black', size=0.5,
                                     linetype='solid'),
          axis.title.x = element_blank(),
          axis.text.x  = element_text(size=12),
          axis.title.y = element_text(size=14), 
          axis.text.y  = element_text(size=12),
          plot.title = element_text(face="italic"),
          legend.position = "none") +
    geom_errorbar(aes(ymin=Mass_loss_per_biomass_day_mg-se,
                      ymax=Mass_loss_per_biomass_day_mg+se),width=.2,
                  position=position_dodge(.9)))

sumDV <-
  summarySE(data=DV_consumption_data_joined, 
            measurevar = c("Mass_loss_per_biomass_day_mg"), 
            groupvars=c("ResDisp"))
sumDV$ResDisp <- factor(sumDV$ResDisp,
                        levels=c("Residents", "Dispersers"),
                        ordered=TRUE)

(DVdispnetplot <- 
    ggplot(sumDV, 
           aes(x = ResDisp, y = Mass_loss_per_biomass_day_mg,
               fill = ResDisp)) + 
    geom_bar(stat='identity', colour="black") +
    geom_point(data = DV_consumption_data_joined, 
               mapping = aes(x = ResDisp, 
                             y = Mass_loss_per_biomass_day_mg),
               fill="darkgray", alpha=0.5, size=3)+
    ylab(bquote('Consumption (mg mg'^-1*' amphipod day'^-1*')')) + 
    ylim(c(0,15)) + ggtitle("Dikerogammarus villosus")+
    scale_fill_manual("Dispersal Type", values=dispcols)+
    annotate('text', x = 2, y = 7, label='*', size=7, fontface="bold")+
    theme_classic() +
    theme(axis.line.x = element_line(colour = 'black', 
                                     size=0.5, linetype='solid'),
          axis.line.y = element_line(colour = 'black', size=0.5,
                                     linetype='solid'),
          axis.title.x = element_blank(),
          axis.text.x  = element_text(size=12),
          axis.title.y = element_text(size=14), 
          axis.text.y  = element_text(size=12),
          plot.title = element_text(face="italic"),
          legend.title = element_text(size=14),
          legend.text = element_text(size=12),
          legend.position = "none") +
    geom_errorbar(aes(ymin=Mass_loss_per_biomass_day_mg-se,
                      ymax=Mass_loss_per_biomass_day_mg+se),width=.2,
                  position=position_dodge(.9)))

grid.arrange(GFdispnetplot, DVdispnetplot, ncol=2)

#### Correlation of density between residents and dispersers ####

# subset the data to things of interest for this comparison
GF_data_wide_density_compare_AD <- GF_consumption_data_joined[,c(2,19,29)]
DV_data_wide_density_compare_AD <- DV_consumption_data_joined[,c(2,13,23)]

# reshape long to wide
GF_density_compare_wide_AD <- 
  spread(GF_data_wide_density_compare_AD, ResDisp, Amphipod_Days)
DV_density_compare_wide_AD <- 
  spread(DV_data_wide_density_compare_AD, ResDisp, Amphipod_Days)

cor.test(GF_density_compare_wide_AD$Dispersers,
         GF_density_compare_wide_AD$Residents)

cor.test(DV_density_compare_wide_AD$Dispersers,
         DV_density_compare_wide_AD$Residents)

#### Checking if residuals associated with DispNet condition types ####

GF_resid_frame <- as.data.frame(resid(GFdispnetlmer_randomintercept_sqrt))

GF_resid_frame <- cbind(GF_consumption_data_joined, GF_resid_frame)

colnames(GF_resid_frame)[31] <- "model_resid"

DV_resid_frame <- as.data.frame(resid(DVdispnetlmer_randomintercept_sqrt))

DV_resid_frame <- cbind(DV_consumption_data_joined, DV_resid_frame)

colnames(DV_resid_frame)[24] <- "model_resid"

# change my labels to the ones from DispNet
# CTR-CTR = RA std., PRED no
levels(GF_resid_frame$Condition)[levels(GF_resid_frame$Condition)=="CTR-CTR"] <- 
  "RA std., PRED no"
levels(DV_resid_frame$Condition)[levels(DV_resid_frame$Condition)=="CTR-CTR"] <- 
  "RA std., PRED no"

# FISH-CTR = RA std., PRED yes
levels(GF_resid_frame$Condition)[levels(GF_resid_frame$Condition)=="FISH-CTR"] <- 
  "RA std., PRED yes"
levels(DV_resid_frame$Condition)[levels(DV_resid_frame$Condition)=="FISH-CTR"] <- 
  "RA std., PRED yes"

# CTR-LOW = RA low, PRED no
levels(GF_resid_frame$Condition)[levels(GF_resid_frame$Condition)=="CTR-LOW"] <- 
  "RA low, PRED no"
levels(DV_resid_frame$Condition)[levels(DV_resid_frame$Condition)=="CTR-LOW"] <- 
  "RA low, PRED no"

# FISH-LOW = RA low, PRED yes
levels(GF_resid_frame$Condition)[levels(GF_resid_frame$Condition)=="FISH-LOW"] <- 
  "RA low, PRED yes"
levels(DV_resid_frame$Condition)[levels(DV_resid_frame$Condition)=="FISH-LOW"] <- 
  "RA low, PRED yes"

# there are some rows which do not have a "Condition" label
# because this came from the weights dataframe
# and these replicates had no surviving amphipods to weigh


GF_resid_frame$Condition[GF_resid_frame$Mesocosm == 10] <- "RA low, PRED yes"
GF_resid_frame$Condition[GF_resid_frame$Mesocosm == 16] <- "RA std., PRED yes"
GF_resid_frame$Condition[GF_resid_frame$Mesocosm == 68] <- "RA low, PRED no"
GF_resid_frame$Condition[GF_resid_frame$Mesocosm == 112] <- "RA std., PRED yes"
GF_resid_frame$Condition[GF_resid_frame$Mesocosm == 125] <- "RA std., PRED yes"
GF_resid_frame$Condition[GF_resid_frame$Mesocosm == 127] <- "RA low, PRED yes"
GF_resid_frame$Condition[GF_resid_frame$Mesocosm == 139] <- "RA low, PRED yes"
GF_resid_frame$Condition[GF_resid_frame$Mesocosm == 152] <- "RA std., PRED no"
GF_resid_frame$Condition[GF_resid_frame$Mesocosm == 155] <- "RA low, PRED yes"
GF_resid_frame$Condition[GF_resid_frame$Mesocosm == 174] <- "RA std., PRED no"
GF_resid_frame$Condition[GF_resid_frame$Mesocosm == 175] <- "RA low, PRED no"

DV_resid_frame$Condition[DV_resid_frame$Mesocosm == 6] <- "RA low, PRED no"
DV_resid_frame$Condition[DV_resid_frame$Mesocosm == 9] <- "RA std., PRED yes"
DV_resid_frame$Condition[DV_resid_frame$Mesocosm == 16] <- "RA std., PRED yes"
DV_resid_frame$Condition[DV_resid_frame$Mesocosm == 18] <- "RA low, PRED yes"
DV_resid_frame$Condition[DV_resid_frame$Mesocosm == 33] <- "RA std., PRED no"
DV_resid_frame$Condition[DV_resid_frame$Mesocosm == 34] <- "RA low, PRED no"
DV_resid_frame$Condition[DV_resid_frame$Mesocosm == 38] <- "RA std., PRED yes"
DV_resid_frame$Condition[DV_resid_frame$Mesocosm == 43] <- "RA std., PRED no"
DV_resid_frame$Condition[DV_resid_frame$Mesocosm == 47] <- "RA low, PRED no"
DV_resid_frame$Condition[DV_resid_frame$Mesocosm == 55] <- "RA std., PRED yes"
DV_resid_frame$Condition[DV_resid_frame$Mesocosm == 58] <- "RA std., PRED no"
DV_resid_frame$Condition[DV_resid_frame$Mesocosm == 60] <- "RA low, PRED yes"
DV_resid_frame$Condition[DV_resid_frame$Mesocosm == 82] <- "RA std., PRED yes"
DV_resid_frame$Condition[DV_resid_frame$Mesocosm == 98] <- "RA std., PRED no"
DV_resid_frame$Condition[DV_resid_frame$Mesocosm == 105] <- "RA low, PRED yes"
DV_resid_frame$Condition[DV_resid_frame$Mesocosm == 110] <- "RA low, PRED no"
DV_resid_frame$Condition[DV_resid_frame$Mesocosm == 111] <- "RA std., PRED no"
DV_resid_frame$Condition[DV_resid_frame$Mesocosm == 117] <- "RA low, PRED yes"
DV_resid_frame$Condition[DV_resid_frame$Mesocosm == 125] <- "RA std., PRED yes"
DV_resid_frame$Condition[DV_resid_frame$Mesocosm == 135] <- "RA low, PRED no"
DV_resid_frame$Condition[DV_resid_frame$Mesocosm == 138] <- "RA low, PRED no"
DV_resid_frame$Condition[DV_resid_frame$Mesocosm == 148] <- "RA std., PRED yes"
DV_resid_frame$Condition[DV_resid_frame$Mesocosm == 153] <- "RA std., PRED yes"
DV_resid_frame$Condition[DV_resid_frame$Mesocosm == 155] <- "RA low, PRED yes"
DV_resid_frame$Condition[DV_resid_frame$Mesocosm == 160] <- "RA low, PRED no"
DV_resid_frame$Condition[DV_resid_frame$Mesocosm == 165] <- "RA low, PRED yes"
DV_resid_frame$Condition[DV_resid_frame$Mesocosm == 174] <- "RA std., PRED no"
DV_resid_frame$Condition[DV_resid_frame$Mesocosm == 175] <- "RA low, PRED no"
DV_resid_frame$Condition[DV_resid_frame$Mesocosm == 177] <- "RA low, PRED no"

# test if residuals differ with dispersal condition
summary(lm(model_resid ~ Condition, data=GF_resid_frame))

summary(lm(model_resid ~ Condition, data=DV_resid_frame))

#### Figure S1 ####

dispnet_colors <- c("RA std., PRED no" = "darkblue",
                    "RA low, PRED no" = "blue",
                    "RA std., PRED yes" = "darkred",
                    "RA low, PRED yes" = "red")

sumGF_resid <-
  summarySE(data=GF_resid_frame, 
            measurevar = c("model_resid"), 
            groupvars=c("Condition"))

(GFdispnet_resid_treatment <- 
    ggplot(sumGF_resid, 
           aes(x = Condition, y = model_resid,
               fill = Condition)) + 
    geom_hline(yintercept=0)+
    geom_bar(stat='identity', color="black", alpha=0.5) +
    ylab("Residual from model") + 
    xlab("DispNet Treatment") + ggtitle("G. fossarum")+
    scale_color_manual("DispNet Treatment", values=dispnet_colors)+
    scale_fill_manual("DispNet Treatment", values=dispnet_colors)+
    theme_classic() +
    theme(axis.line.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.line.y = element_line(colour = 'black', size=0.5,
                                     linetype='solid'),
          axis.title.x = element_text(size=12),
          axis.text.x  = element_text(size=10, color="black"),
          axis.title.y = element_text(size=12), 
          axis.text.y  = element_text(size=10, color="black"),
          plot.title = element_text(face="italic"),
          legend.position = "none") +
    geom_point(data = GF_resid_frame, 
               mapping = aes(x = Condition, y = model_resid),
               fill="gray", alpha=0.3, size=3)+
    geom_errorbar(aes(ymin=model_resid-se,
                      ymax=model_resid+se),width=.2,
                  position=position_dodge(.9), color="black"))

sumDV_resid <-
  summarySE(data=DV_resid_frame, 
            measurevar = c("model_resid"), 
            groupvars=c("Condition"))

(DVdispnet_resid_treatment <- 
    ggplot(sumDV_resid, 
           aes(x = Condition, y = model_resid,
               fill = Condition)) + 
    geom_hline(yintercept=0)+
    geom_bar(stat='identity', color="black", alpha=0.5) +
    ylab("Residual from model") + 
    xlab("DispNet Treatment") + ggtitle("D. villosus")+
    scale_color_manual("DispNet Treatment", values=dispnet_colors)+
    scale_fill_manual("DispNet Treatment", values=dispnet_colors)+
    theme_classic() +
    theme(axis.line.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.line.y = element_line(colour = 'black', size=0.5,
                                     linetype='solid'),
          axis.title.x = element_text(size=12),
          axis.text.x  = element_text(size=10, color="black"),
          axis.title.y = element_text(size=12), 
          axis.text.y  = element_text(size=10, color="black"),
          plot.title = element_text(face="italic"),
          legend.position = "none") +
    geom_point(data = DV_resid_frame, 
               mapping = aes(x = Condition, y = model_resid),
               fill="gray", alpha=0.3, size=3)+
    geom_errorbar(aes(ymin=model_resid-se,
                      ymax=model_resid+se),width=.2,
                  position=position_dodge(.9), color="black"))

grid.arrange(GFdispnet_resid_treatment, DVdispnet_resid_treatment,
             nrow=2)

#### Survival rates #####

GF_consumption_data_joined$final_alive <-
  ifelse(is.na(GF_consumption_data_joined$Final_Alive_28_Nov) == FALSE, 
         GF_consumption_data_joined$Final_Alive_28_Nov,
         GF_consumption_data_joined$Final_Alive_29_Nov)

GF_surv <- sum(GF_consumption_data_joined$final_alive)/
  sum(GF_consumption_data_joined$Alive_11_Nov)

DV_surv <- sum(DV_consumption_data_joined$Final_Alive_8_Feb)/
  sum(DV_consumption_data_joined$Start_Alive_26_Jan)
