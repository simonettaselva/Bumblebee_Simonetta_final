################################################
# Compute Functional Diversity and Phylogentic Species Metrics
# by Simonetta Selva
#
# Created: May 2nd, 2023
# Project: Bumblebee 2022
################################################

# Here, all Functional Diversity Metrics and the Phylogenetic Species Variability is calculated 
# per individual bumblebee or if wanted per site. They are needed in all other script. 

# information: 
# 1) every subsection works in itself
# 2) all the plots not used in the main work are in the script but # are used to not print them

# FUNCTIONAL DIVERSITY METRICS ----

## preparation ----
# clear environment
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

# load data
BB22_full <- read_csv("BB22_full.csv") %>% 
  mutate(ID = as.character(ID))
BB22_full$ID.short = as.factor(substring(BB22_full$ID,1, nchar(BB22_full$ID)-1)) #remove information on leg and body from ID

# choose only to numeric data
BB22_full.numeric <- BB22_full %>% 
  summarise(ID.short = ID.short, #remove information on leg and body from ID
            ID = as.factor(ID),
            location = as.factor(location),
            landscape = as.factor(landscape),
            replicate = as.factor(replicate),
            bbspecies = as.factor(bbspecies),
            bborgan = as.factor(bborgan),
            site = as.factor(paste(location, landscape, replicate, sep="")),
            plant.species = plant.species,
            binom.abund = binom.abund,
            abundance = Abundance,
            Flowering_duration = Flowering_months_duration,
            Flowering_start = start_flowering,
            growth_form_numeric = growth_form_numeric,
            structural_blossom_numeric = structural_blossom_numeric,
            sugar.concentration = sugar.concentration,
            symmetry_numeric = symmetry_numeric,
            plant_height_m = plant_height_m)
BB22_full.numeric$Flowering_duration <- as.numeric(BB22_full.numeric$Flowering_duration)

## Correlation analysis ----
# correlation analysis of the variables
library(Hmisc)
library(corrplot)
res <- rcorr(as.matrix(BB22_full.numeric[,c(12:18)]),
             type="pearson")
M <- cor(BB22_full.numeric[,c(12:18)], 
         use = "complete.obs")

### FIGURE S3 ----
corrplot::corrplot(M, type="upper", order="hclust", 
                   p.mat = res$P, sig.level = 0.01, tl.col = "black",
                   col = COL2('RdBu', 10)) # plot correlation with p-values

# remove plant height, flowering start and symmetry and look at correlation
res2 <- rcorr(as.matrix(BB22_full.numeric[,c(12, 14, 15, 16)]), type="pearson") # all p values are <0.05
M <- cor(BB22_full.numeric[,c(12, 14, 15, 16)], use = "complete.obs")
corrplot::corrplot(M, type="upper", order="hclust", 
                   p.mat = res2$P, sig.level = 0.01, tl.col = "black",
                   col = COL2('RdBu', 10)) # plot correlation with p-values

# remove plant height, flowering start and symmetry from data set
BB22_full.red <- BB22_full.numeric%>% 
  dplyr::select(-Flowering_start,
                -symmetry_numeric,
                -plant_height_m) 

## compute FDs ----
library(caret)
library(vegan)
library(reshape2)
library(FD)
library(GGally)

# the following code has to be run 4 times with each time different input variables i and j 
# to calculate FD for each bumblebee species and on a site or ID level
# i = "B.lapidarius" or "B.pascuorum"
# j = "site" or "ID.short"
# It can also be performed in a loop, but due to computational time it is easier to them individually. 

i <- "B.lapidarius"
j <- "ID.short"

# filter for bumblebee species, remove plant species entries that do not have traits and prepare data
BB22_full.loop <- BB22_full.red[BB22_full.red$bbspecies == i,]%>%
  filter(plant.species!= "Fabaceae sp.", 
         plant.species!= "Cyclamen sp.",
         plant.species!= "Petunia sp.",
         plant.species!= "Mandevilla atroviolacea" )%>%
  mutate(Flowering_duration = as.numeric(Flowering_duration))%>% 
  droplevels() # remove duplicates

# select columns used in computing FDs
BB22_full.loop.species <- BB22_full.loop %>% 
  dplyr::select(plant.species, Flowering_duration, structural_blossom_numeric, sugar.concentration) %>% 
  distinct() # remove duplicates

# impute missing data
trt.mis.pred <- preProcess(as.data.frame(BB22_full.loop.species[,-c(1)]), "knnImpute")
traits <- predict(trt.mis.pred, BB22_full.loop.species[,-c(1)]); head(trt.mis.pred)
traits <- as.data.frame(traits)
rownames(traits)  <- BB22_full.loop.species$plant.species
traits <- traits[order(row.names(traits)), ] # reorder traits into alphabetical order (needed for FD package)

# relative abundance data
j.unquoted <- rlang::sym(j)

BB22_full.ab <- BB22_full.loop%>%
  group_by(!!j.unquoted, plant.species)%>%
  summarise(abundance = sum(abundance))%>% 
  distinct() # remove duplicates

# create new data frame that recalculates relative abundance data with sum of leg and body pollen
BB22_full.ab.new <- c()
for (h in unique(BB22_full.ab[[j]])) {
  temp <- BB22_full.ab[BB22_full.ab[[j]]==h,]
  perc <- sum(temp$abundance)
  for (k in 1:nrow(temp)) {
    temp$ab.new[k] <- 100/perc*temp$abundance[k]/100
  }
  BB22_full.ab.new <- rbind(BB22_full.ab.new, temp[, c(1,2,4)])
}

# bring into wide format
wide <- dcast(BB22_full.ab.new, BB22_full.ab.new[[j]] ~ plant.species, value.var="ab.new")[,-1]
rownames(wide)  <- levels(BB22_full.ab.new[[j]])
wide[is.na(wide)] <- 0

# create matrix (abundance and presence/absence)
# abundance
sp.ab <- as.matrix(wide)
write.table(sp.ab, file = paste("./FD/community_matrix_", i, "_", j, ".txt", sep=""), sep = ",")

# presence/absence
sp.pa <- decostand(wide, "pa")
sp.pa <- as.matrix(sp.pa) #turn into matrix

# compute FD
fd.weig <- FD::dbFD(x = traits , a = sp.ab, w.abun = T) # weighted 
fd.bino <- FD::dbFD(x = traits , a = sp.pa) # not weighted 

df.FD <- data.frame(nbsp.w = fd.weig$nbsp,
                    nbsp = fd.bino$nbsp,
                    FRic.w = fd.weig$FRic, 
                    FRic = fd.bino$FRic, 
                    FEve.w = fd.weig$FEve, 
                    FEve = fd.bino$FEve, 
                    FDiv.w = fd.weig$FDiv, 
                    FDiv = fd.bino$FDiv)

# look at possible correlation between weighted and non.weighted FDs
ggpairs(df.FD)

cor(df.FD$FRic, df.FD$FRic.w, method=c("pearson"), use = "complete.obs")
cor(df.FD$FEve, df.FD$FEve.w, method=c("pearson"), use = "complete.obs")
cor(df.FD$FDiv, df.FD$FDiv.w, method=c("pearson"), use = "complete.obs")

# assign and export data frame
assign(paste("df.FD", i, sep="_"), df.FD[, c(1,3,5,7)])
write.csv(assign(paste("df.FD",i,j, sep = "_"), df.FD[, c(1,3,5,7)]), 
          file = paste("./FD/FD_", i, "_", j, ".csv", sep = ""))


# PHYLOGENETIC SPECIES METRICS ----

# clear environment
rm(list=ls())

#load libraries
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(V.PhyloMaker2)

# set working directory to main repository
input <- "YOUR PATH/Bumblebee_Simonetta_final/01_DATA" #change „YOUR PATH“ to where the repo is
output <- "YOUR PATH/Bumblebee_Simonetta_final/03_OUTPUT" #change „YOUR PATH“ to where the repo is
setwd(input)


# load data
BB22.full <- read_csv("BB22_full.csv") # import data set with

# create a data frame with taxa (species, genus, family)
phylo <- data.frame(species = BB22.full$plant.species, genus = BB22.full$genus, family = BB22.full$family)

# phylogenetic hypotheses under three scenarios based on a backbone phylogeny
tree.result <- phylo.maker(phylo)
write.tree(tree.result$scenario.3, "BB22_plant_tree.tre") # save for later purpose

# Tree result needs to be re-imported
tree <- ape::read.tree("BB22_plant_tree.tre")

## compute PD ----
# loop for two species and two resolution levels (site and ID)
species <- c("B.lapidarius", "B.pascuorum")
level <- c("ID.short") # or level <- c("site", ID.short") if also site level is needed

for (i in species) {
  for (j in level) {
    sp.ab <- read.delim(paste("./FD/community_matrix_", i, "_", j, ".txt", sep=""), sep = ",", 
                        check.names = F)
    colnames(sp.ab) <- sub(" ", "_", colnames(sp.ab)) # match plant species names
    
    library(picante)
    PD_var <- psv(sp.ab,tree,compute.var=TRUE,scale.vcv=TRUE)
    PD_ric <- psr(sp.ab,tree,compute.var=TRUE,scale.vcv=TRUE)
    PD_eve <- pse(sp.ab,tree,scale.vcv=TRUE)
    PD_clu <- psc(sp.ab,tree,scale.vcv=TRUE)
    
    write.csv(assign(paste("PD", i, j, sep = "_"), 
                     data.frame(ID = rownames(PD_var),
                                nbsp = PD_var$SR,
                                PVar = PD_var$PSVs,
                                PRic = PD_ric$PSR,
                                PEve = PD_eve$PSEs,
                                PClu = PD_clu$PSCs)), 
              file = paste("./PD/PD_", i, "_", j, ".csv", sep = ""))
  }
}



