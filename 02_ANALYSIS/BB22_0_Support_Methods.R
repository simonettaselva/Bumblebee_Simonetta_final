################################################
# Support Methods
# by Simonetta Selva
#
# Created: May 2nd, 2023
# Project: Bumblebee 2022
################################################

# information: 
# 1) every subsection works in itself
# 2) all the plots not used in the main work are in the script but # are used to not print them

# SITES OVERVIEW

## preparation ----
# clear environment
rm(list=ls())

#load libraries
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)

# set working directory to main repository
input <- "~/Library/CloudStorage/GoogleDrive-simo1996s@gmail.com/My Drive/ETH/Master Thesis/bumblebee_try/01_DATA"
output <- "~/Library/CloudStorage/GoogleDrive-simo1996s@gmail.com/My Drive/ETH/Master Thesis/bumblebee_try/03_OUTPUT"

setwd(input)

# load data
BB22.habitat <- read_csv("habitat maps/BB22_habitat_buffer.csv") %>%
  mutate(type = substring(as.character(TypoCH_NUM), 1, 1))
BB22.habitat$type <- as.factor(BB22.habitat$type)

## analysis ----
# rename levels according to Price et al., 2021
levels(BB22.habitat$type) <- c("Water body", "Riparian and wetland", "Glacier, rock, rubble and scree", 
                               "Greenland", "Bush and shrubbery", "Forest", "Turf and ruderal areas", 
                               "Tree and field crops", "Building")
# rename sites to match
BB22.habitat$site <- factor(BB22.habitat$site, levels = c("ZHUA", "ZHUB", "ZHUC", "ZHRD", "ZHRE", "ZHRF", "BEUA", "BEUB", "BEUC", "BERD", "BSUA", "BSUB", "BSUC", "BSRD", "BSRE", "BSRF"))

# calculate sum of each landcover type per site
BB22.habitat.area <- BB22.habitat %>%
  group_by(site, type)%>%
  summarise(location = Location,
            landscape = Habitat,
            replicate = Code,
            area = sum(area)) %>%
  distinct()

### FIGURE S1c ----
# plotting
ggplot(BB22.habitat.area, aes(fill=type, y=area, x=site)) + 
  geom_bar(position="stack", stat="identity") +
  ylab("area [m2]") + xlab("")+
  scale_fill_manual(values = c("#3470f3", "#d4f7fe", "#d8e7f1", "#f8e3cd", "#f0f04f", "#b5e674", "#c09853", "#d779de", "#a0a0a0"))+
  theme_classic(base_size = 20) +  
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
setwd(output)
ggsave(paste("./preparation/habitat map/sites_habitat_map.png", sep = ""), width = 16, height = 10)
setwd(input)

# SUGAR DATA

## preparation ----
# clear environment
rm(list = ls())

# load required library
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)

# set working directory to main repository
input <- "~/Library/CloudStorage/GoogleDrive-simo1996s@gmail.com/My Drive/ETH/Master Thesis/bumblebee_try/01_DATA"
output <- "~/Library/CloudStorage/GoogleDrive-simo1996s@gmail.com/My Drive/ETH/Master Thesis/bumblebee_try/03_OUTPUT"

# load data
setwd(input) # set working directory
BB22.full <- read_csv("BB22_full.csv") # import data set with

# most abundant plants species to justify use of sugar data
BB22_plants_abund <- BB22.full %>%
  filter(binom.abund != 0) 
BB22_plants_abund  <- BB22_plants_abund %>%
  group_by(plant.species)%>%
  summarise(cum.rel.abundance = sum(Abundance),
            sum.occurance = sum(binom.abund))

# see  on which of the most abundant species there is sugar data
# ABUNDANCE
BB22_plants_abund[,"abund.sum"] <- NA
BB22_plants_abund[,"abund.bin"] <- NA

BB22_plants_abund <- BB22_plants_abund%>%                                     
  arrange(desc(cum.rel.abundance))

#find the 95% of the distribution
for (i in 1:228){
  if (i == 1) {
    BB22_plants_abund$abund.sum[i] <- BB22_plants_abund$cum.rel.abundance[i]
  } else {
    BB22_plants_abund$abund.sum[i] <- BB22_plants_abund$cum.rel.abundance[i] + BB22_plants_abund$abund.sum[i-1]
  }}

for (i in 1:228){
  if (BB22_plants_abund$abund.sum[i] < sum(BB22_plants_abund$cum.rel.abundance)*0.95) {
    BB22_plants_abund$abund.bin[i] <- 1
  } else {
    BB22_plants_abund$abund.bin[i] <- 0
  }}

BB22_plants_abund$abund.bin <- as.factor(BB22_plants_abund$abund.bin)

# setwd(output)
# ggplot(BB22_plants_abund, aes(x=reorder(plant.species, -cum.rel.abundance), y=cum.rel.abundance, fill=abund.bin)) +
#   geom_bar(stat="identity") +
#   labs(title="Abundance",x ="plant species", y = "sum of relative abundance") +  
#   theme_bw() + coord_flip() + scale_fill_manual(values=c("red", "black"))
# ggsave("plant_species_sum_abundance.png", width = 8, height = 22)
# setwd(input)

# refine plot: show only the most abundant 29 species and summarize the rest into one
BB22_plants_abund_cut <- BB22_plants_abund %>%
  arrange(desc(cum.rel.abundance)) %>% 
  summarise(plant.species = plant.species,
            cum.rel.abundance = cum.rel.abundance) %>% 
  slice(1:29)
BB22_plants_abund_cut <- rbind(BB22_plants_abund_cut, c("other (199 species)", 1))
BB22_plants_abund_cut$cum.rel.abundance <- as.numeric(BB22_plants_abund_cut$cum.rel.abundance)

# check if 29 most abundant have sugar data
BB22_plant_sugar <- BB22.full%>%  
  summarise(plant.species = plant.species,
            sugar.concentration = sugar.concentration) %>%  
  filter(plant.species %in% BB22_plants_abund_cut$plant.species)

BB22_plants_abund_cut[,"sugar.data"] <- NA
BB22_plants_abund_cut <- BB22_plants_abund_cut%>%                                     
  arrange(plant.species)

BB22_plant_sugar <- BB22_plant_sugar%>%                                     
  arrange(plant.species)

# add binary variable (0 - no sugar data available, 1 - sugar data available)
for (i in 1:30){
  if (is.na(BB22_plant_sugar$sugar.concentration[i])) {
    BB22_plants_abund_cut$sugar.data[i] <- 0
  } else {
    BB22_plants_abund_cut$sugar.data[i] <- 1
  }}
BB22_plants_abund_cut$sugar.data <- as.factor(BB22_plants_abund_cut$sugar.data)

### FIGURE S2 ----
setwd(output)
ggplot(BB22_plants_abund_cut, aes(x = reorder(plant.species, -cum.rel.abundance), y=cum.rel.abundance, fill = sugar.data)) +
  geom_bar(stat="identity") +
  labs(title="",x ="plant species", y = "sum of relative abundance") +  
  theme_classic(base_size = 20) + 
  theme(axis.text.x = element_text(angle = 90))+ 
  scale_fill_discrete(name = "sugar data", labels=c('absent', 'present'))
ggsave("./preparation/sugar/plant_species_abundance_sugar.png", width = 15, height =8)
setwd(input)

# percentages of the most abundant species
total.abund <- sum(BB22_plants_abund_cut$cum.rel.abundance)
most.abund <- sum(BB22_plants_abund_cut$cum.rel.abundance[order(-BB22_plants_abund_cut$cum.rel.abundance)][1:3])

perc <- most.abund/total.abund #the 3 most abundant species make up ~79% of all abundances



