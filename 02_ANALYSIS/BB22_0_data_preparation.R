################################################
# R data preparation
# by Simonetta Selva
#
# Created: May 2nd, 2023
# Project: Bumblebee 2022
################################################

# PREPARE FULL DATASET ----

# Here, the data set calles BB22_full is prepared. It is needed for all other scripts

# clear work environment
rm(list = ls())

#load libraries
library(dplyr)
library(tidyverse)

# set working directory to main repository
input <- "YOUR PATH/Bumblebee_Simonetta_final/01_DATA" #change „YOUR PATH“ to where the repo is
output <- "YOUR PATH/Bumblebee_Simonetta_final/03_OUTPUT" #change „YOUR PATH“ to where the repo is

# load data
setwd(input)
BB22 <- read.csv("BB22_data_tres_0.01.csv")

#create binom. abundance variable
BB22 <- BB22 %>%
  mutate(binom.abund = if_else(Abundance == 0, 0, 1),
         ID = substring(Sample, 2),
         site = paste(location, landscape, replicate, sep = ""),
         bbspecies = as_factor(bbspecies)) %>%
  select(-Sample, -X, -project, -xID)
levels(BB22$bbspecies) <- c("B.pascuorum", "B.lapidarius") # rename BB species
BB22 <- BB22[, c(15, 16, 1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 2, 14)] # reorder columns
BB22$OTU <- recode(BB22$OTU, # match plant names
                   "Phedimus spurius" = "Sedum spurium",
                   "Phedimus aizoon" = "Sedum aizoon",
                   "Clematis sp. YX-2018" = "Clematis sp.",
                   "Salvia amplexicaulis" = "Salvia nemorosa",
                   "Phedimus kamtschaticus" = "Sedum aizoon",
                   "Tilia americana x Tilia x moltkei" = "Tilia americana",
                   "Hylotelephium telephium" = "Sedum telephium",
                   "Petunia sp. LR-2018" = "Petunia sp.",
                   "Onopordum illyricum" = "Onopordum acanthium",
                   "Fabaceae spc" = "Fabaceae sp.",
                   "Potentilla glabra" = "Dasiphora fruticosa",
                   "Linaria spc" = "Linaria sp.",
                   "Begonia spc" = "Begonia sp.",
                   "Lamiaceae spc" = "Lamiaceae sp.",
                   "Allium spc" = "Allium sp.",
                   "Asteraceae spc" = "Asteraceae sp.",
                   "Boraginaceae spc" = "Boraginaceae sp.",
                   "Betonica officinalis" = "Stachys officinalis",
                   "Cyclamen spc" = "Cyclamen sp.",
                   "Crepis spc" = "Crepis sp.",
                   "Eruca pinnatifida" = "Eruca vesicaria",
                   "Sedum montanum" = "Sempervivum montanum",
                   "Trifolium spc" = "Trifolium sp.",
                   "Lotus spc" = "Lotus sp.",
                   "Hypochaeris spc" = "Hypochaeris sp.",
                   "Nepeta spc" = "Nepeta sp.",
                   "x Chitalpa tashkentensis" = "xChitalpa tashkentensis")

# 230 species
# 47 families

#remove entries with abundance = 0
setwd(input)
BB22.abund <- BB22 %>% 
  dplyr::filter(binom.abund == 1)

# add information on plants species (file BB22_data_plants.R)
BB22_plant_traits <- read_csv("BB22_plant_traits.csv")
BB22 <- BB22 %>% rename(plant.species = OTU)

# join the two dataframes into one -> BB22_full
BB22_full <- full_join(BB22, BB22_plant_traits, by = "plant.species") %>%
  dplyr::select(-Family, -Genus, -species) %>%
  dplyr::rename(species = Sp2) %>%
  dplyr::filter(binom.abund == 1) %>%
  dplyr::mutate(native_exotic = as_factor(native_exotic),
         pollination_mode = as_factor(pollination_mode),
         growth_form_category = as_factor(growth_form_category),
         inflorescence = as_factor(inflorescence),
         structural_blossom_class = as_factor(structural_blossom_class),
         symmetry = as_factor(symmetry),
         Flowering_months_duration = as.numeric(Flowering_months_duration))
BB22_full <- BB22_full[, c(1, 2, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 16, 3, 14, 15, 
                           17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30)] # reorder columns

# create numerical values for factors
BB22_full$growth_form_numeric <- BB22_full$growth_form_category
levels(BB22_full$growth_form_numeric) <- 1:5
BB22_full$growth_form_numeric <- as.numeric(BB22_full$growth_form_numeric)
# herb = 1, shrub = 2, climber = 3, tree = 4 , shrub/tree = 5

BB22_full$structural_blossom_numeric <- BB22_full$structural_blossom_class
levels(BB22_full$structural_blossom_numeric) <- 1:7
BB22_full$structural_blossom_numeric <- as.numeric(BB22_full$structural_blossom_numeric)
# flag = 1, gullet = 2, dish_bowl = 3, stalk_disk = 4, bell_trumpet = 5, tube = 6, brush = 7

BB22_full$symmetry_numeric <- BB22_full$symmetry
levels(BB22_full$symmetry_numeric) <- 1:3
BB22_full$symmetry_numeric <- as.numeric(BB22_full$symmetry_numeric)
# sigomorph = 1, actinomorph = 2, no_symmetry = 3

# export as csv
write.csv(BB22_full, "BB22_full.csv")


# PREPARE SITE PLANT OCCURANCES FROM GBIF AND INFOFLORA ----

# Here, combined data sets form the GBIF and InfoFlora data per site is generated and exported.

# note: the code above needs to have run already!!

# clear work environment
rm(list=ls())

#load libraries
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)

# set working directory to main repository
input <- "YOUR PATH/Bumblebee_Simonetta_final/01_DATA" #change „YOUR PATH“ to where the repo is
output <- "YOUR PATH/Bumblebee_Simonetta_final/03_OUTPUT" #change „YOUR PATH“ to where the repo is
setwd(input)

## species list per site from pollen data ----
# bumblebee data: create species lists for each site and save in list()
sitenames <- c("ZHUA", "ZHUB", "ZHUC", "ZHRD", "ZHRE", "ZHRF", "BEUA", "BEUB", 
               "BEUC", "BERD", "BSUA", "BSUB", "BSUC", "BSRD", "BSRE", "BSRF")
# region <- c("ZHU", "ZHR", "BSU", "BSR", "BEU", "BER")

BB22.full <- read_csv("BB22_full.csv") %>%
  mutate(site = as_factor(site),
         plant.species = as_factor(plant.species),
         region = substr(site, 1, 3)) 

# create species lists of pollen per site
site.list <- list()
for (i in sitenames) {
  site.list[[i]]  <- unique(BB22.full$plant.species[BB22.full$site == i]) %>%
    droplevels()
}

save(site.list, file = "site_list_bb.RData")


## GBIF and InfoFlora per site in 1500m Buffer ----
sitenames <- c("ZHUA", "ZHUB", "ZHUC", "ZHRD", "ZHRE", "ZHRF", "BEUA", "BEUB", 
               "BEUC", "BERD", "BSUA", "BSUB", "BSUC", "BSRD", "BSRE", "BSRF")

### InfoFlora ----
# prepare list to store
site.list.InfoFlora.1500 <- list()

# create species list per site for InfoFlora
for (i in sitenames) {
  # import InfoFlora data and clean it
  site.data.InfoFlora.1500  <- read_csv(paste("./GBIF and InfoFlora/InfoFlora_",i , "_1500.csv", sep = "")) 
  site.data.InfoFlora.1500 <- site.data.InfoFlora.1500 %>%
    summarise(Species = Taxon,
              site = rep(i, nrow(site.data.InfoFlora.1500)),
              region = rep(substr(i, 1, 3), nrow(site.data.InfoFlora.1500)))
  j <- 0
  site.data.InfoFlora.1500$species <- c()
  for (j in 1:nrow(site.data.InfoFlora.1500)){
    # adjust the plant species names since they dont match between GBIF and InfoFlora
    site.data.InfoFlora.1500$species[j] <- paste(unlist(strsplit(site.data.InfoFlora.1500$Species[j], 
                                                                 split=' ', fixed=TRUE))[1],
                                                 unlist(strsplit(site.data.InfoFlora.1500$Species[j], 
                                                                 split=' ', fixed=TRUE))[2])
  }
  # store in list
  site.list.InfoFlora.1500[[i]] <- unique(site.data.InfoFlora.1500$species)
}

### GBIF ----
# prepare list to store
site.list.gbif.1500 <- list()

# create species list per site for GBIF
for (i in sitenames) {
  # import GBIF data and clean it
  site.data.gbif.1500  <- read_csv(paste("./GBIF and InfoFlora/GBIF_", i, "_1500.csv", sep = "")) %>%
    mutate(species = as_factor(species)) %>%
    filter(class == "Magnoliopsida" | class == "Liliopsida") # filter for classes that are potentially visited from bumblebees
  # store in list
  site.list.gbif.1500[[i]] <- unique(site.data.gbif.1500$species) 
}

# combine species list GBIF and InfoFlora (site level)
site.list.1500 <- list()
for (i in sitenames) {
  x <- as.factor(site.list.InfoFlora.1500 [[i]])
  y <- as.factor(site.list.gbif.1500 [[i]])
  site.list.1500[[i]] <- unique(c(x,y)) %>%
    droplevels()
}
# save list
save(site.list.1500, file="sp_list_gbif_infoflora.RData")


## plot of overview of percentage of visited plant by the bumblebees ----
# prepare list
ratios.1500 <- c()

# create data frame with ratios (visited/non-visited)
for (i in sitenames) {
  shared <- length(intersect(site.list[[i]], site.list.1500[[i]]))
  occured <- length(site.list.1500[[i]])
  first <- c(i, shared/occured, "percentage of plants visited by bumblebees")
  second<- c(i, (occured-shared)/occured, "GBIF and InfoFlora (1500m)")
  ratios.1500 <- rbind(ratios.1500, first, second)
  colnames(ratios.1500)<- c("site", "percentage", "group")
  ratios.1500 <- as.data.frame(ratios.1500)%>% 
    mutate(percentage = as.numeric(percentage))
}

# plot and export
setwd(output)
ggplot(ratios.1500, aes(fill=group, y=percentage, x=site)) + 
  geom_bar(position="fill", stat="identity") + xlab("sites") +
  ggtitle("Comparison of avaible plants species and plants visited by bumblebees") + theme_classic()+
  theme(axis.text.x = element_text(angle = 90))
ggsave("./preparation/GBIF and InfoFlora/Comp_1500_visited_per_site.jpeg", width = 8, height = 8)
setwd(input)





