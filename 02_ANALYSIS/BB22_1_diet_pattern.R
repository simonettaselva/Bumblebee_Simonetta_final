################################################
# Foraging patterns of the two bumblebee species
# by Simonetta Selva
#
# Created: May 2nd, 2023
# Project: Bumblebee 2022
################################################

# HYPOTHESIS: The diet of bumblebees in urban areas are more diverse compared 
# to rural areas and consequent-ly bumblebees in cities will be less dependent 
# on the presence of certain plant species

# information:
# 1) every subsection works in itself
# 2) all the plots not used in the main work are in the script but # are used to not print them


# SPECIES ABBUNDANCES IN POLLEN AND PHYLOGENETIC TREE ----
## preparation ----
# clear environment
rm(list = ls())

# load required library
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(pals)
library(V.PhyloMaker2)
library(ape)

# set working directory to main repository
input <- "~/Library/CloudStorage/GoogleDrive-simo1996s@gmail.com/My Drive/ETH/Master Thesis/bumblebee_try/01_DATA"
output <- "~/Library/CloudStorage/GoogleDrive-simo1996s@gmail.com/My Drive/ETH/Master Thesis/bumblebee_try/03_OUTPUT"

# load data
setwd(input) # set working directory
BB22.full <- read_csv("BB22_full.csv") # import data set with

# initialize region as column
BB22.full$region <- c()

# add site and region as columns
for (i in 1:nrow(BB22.full)) {
  BB22.full$site[i] <- paste(BB22.full$location[i], BB22.full$landscape[i], BB22.full$replicate[i], sep = "")
  BB22.full$region[i] <- paste(BB22.full$location[i], BB22.full$landscape[i], sep = "")
}

## body and corbicula pollen ----

# 1. create a data frame with taxa (species, genus, family)
phylo <- data.frame(species = BB22.full$plant.species, genus = BB22.full$genus, family = BB22.full$family)

# 2. phylogentic tree
# phylogenetic hypotheses under three scenarios based on a backbone phylogeny
tree.result <- phylo.maker(phylo) # new version of package
write.tree(tree.result$scenario.3, "BB22_plant_tree.tre") # save for later purpose

# plot the phylogenies with node ages displayed with sceanario 3
tree <- plot.phylo(tree.result$scenario.3,
                   cex = 0.5,
                   main = "Phylogenetic tree of species in pollen")
tree # print the tree

# get order of species in tree
phylo.order <- data.frame(sps = tree.result$scenario.3$tip.label)
phylo.order$order <- seq(1, length(phylo.order$sps))
colnames(phylo.order) <- c("plant.species", "order")
phylo.order$plant.species <- sub("_", " ", phylo.order$plant.species)

# create new data frame with needed variables
BB22.full.bubble <- BB22.full %>%
  dplyr::group_by(site, plant.species, bbspecies) %>%
  dplyr::summarise(Abundance = sum(Abundance),
                   region = region,
                   landscape = landscape)

# find how many taa were visited per region
BB22.full.taxa <- BB22.full.bubble %>%
  dplyr::group_by(region, bbspecies) %>%
  dplyr::summarise(plant.species = plant.species) %>%
  distinct()

BB22.full.taxa <- BB22.full.taxa %>%
  dplyr::summarise(Nr.taxa = n())

write_csv(BB22.full.taxa, "BB22_NrTaxa_region.csv")


# 3. plotting

# merge two data frames and order along plant species in phylo-tree
BB22.full.bubble_ordered <- merge(x = BB22.full.bubble, y = phylo.order, by.x = "plant.species") %>%
  mutate(plant.species = as_factor(plant.species))
BB22.full.bubble_ordered <- BB22.full.bubble_ordered[order(BB22.full.bubble_ordered$order), ]
BB22.full.bubble_ordered$plant.species <- factor(BB22.full.bubble_ordered$plant.species,
                                                 levels = unique(BB22.full.bubble_ordered$plant.species[order(BB22.full.bubble_ordered$order)]))

# plot bubble plot with relative abundances along order of species in tree
setwd(output)

# on site level
### FIGURE 1 ----
# region level
palette.landscape <- c("#E69F00", "#56B4E9") #create color palette for landscape
ggplot(BB22.full.bubble_ordered, aes(x = region, y = BB22.full.bubble_ordered$plant.species)) +
  geom_point(aes(size = Abundance, color = landscape), alpha = 0.5) +
  facet_wrap(~bbspecies) +
  labs(y = "plant species") +
  theme_classic(base_size = 20) +
  guides(alpha = "none") +
  scale_color_manual(values = palette.landscape, labels = c("rural", "urban"), name = "Landscape") +
  guides(color = guide_legend(override.aes = list(alpha = 1)))
ggsave(paste("./diet pattern/species abundance/Phylo_Bubble_Region.png", sep = ""), width = 16, height = 16)

# POLLEN COMPOSITION ----
## preparation ----

# reset environment
rm(list = ls())

#load libraries
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)

# set working directory to main repository
input <- "~/Library/CloudStorage/GoogleDrive-simo1996s@gmail.com/My Drive/ETH/Master Thesis/bumblebee_try/01_DATA"
output <- "~/Library/CloudStorage/GoogleDrive-simo1996s@gmail.com/My Drive/ETH/Master Thesis/bumblebee_try/03_OUTPUT"

# load data
setwd(input) # set working directory
BB22.full <- read_csv("BB22_full.csv") # import data set

# initialize region as column
BB22.full$region <- c()

# add site and region as columns
for (i in 1:nrow(BB22.full)) {
  BB22.full$site[i] <- paste(BB22.full$location[i], BB22.full$landscape[i], BB22.full$replicate[i], sep = "")
  BB22.full$region[i] <- paste(BB22.full$location[i], BB22.full$landscape[i], sep = "")
}

## analysis ----
# divide data set into leg and body pollen
BB22.full.leg <- BB22.full[BB22.full$bborgan == "L", ]
BB22.full.body <- BB22.full[BB22.full$bborgan == "B", ]

## body and corbicula pollen ----
### plant families per site, region and species ----

# 1. find the 30 most abundant families
families.overview <- BB22.full %>%
  dplyr::group_by(family) %>%
  dplyr::summarise(cum.abund = sum(Abundance))
tot <- sum(families.overview$cum.abund)
families.overview <- families.overview %>%
  dplyr::summarise(family = family,
                   cum.abund = cum.abund,
                   rel.abund = cum.abund / tot)
# plot family abundance
ggplot(families.overview, aes(y = reorder(family, -cum.abund), x = cum.abund)) +
  geom_bar(stat = "identity") + theme_classic()

# use cumulative abundance to select 30 most abundant species
rare.families <- families.overview %>%
  top_n(nrow(families.overview) - 30, -cum.abund)
BB22.full$family.agg <- BB22.full$family

# group all rare families into "other families"
for(h in rare.families$family) {
  for(i in 1:nrow(BB22.full)) {
    if(BB22.full$family.agg[i] == h) {
      BB22.full$family.agg[i] <- "Other families"
    } 
  } # end loop i
} # end loop h

# export data set with new family categories
write_csv(BB22.full, "BB22.full.family")

# 2. produce data frame with families' abundances per site with leg pollen
families.site <- BB22.full %>%
  dplyr::group_by(family.agg, site, bbspecies) %>%
  dplyr::summarise(cum.abund = sum(Abundance))
families.site$family.agg <- as.factor(families.site$family.agg)

# 3. plot families per site
# color palette for families
palette.fams = c("#07575B", "#5D535E", "#C4DFE6", "#336B87", "#FAAF08", "#DFE166", "#1995AD", 
                        "#4897D8", "#80BD9E", "#FA812F", "#66A5AD", "#F34A4A", "#375E97", "#258039",  
                        "#73605B", "#DDBC95", "#A3A599", "#A1D6E2", "#88A550", "#75B1A9", "#D9B44A", 
                        "#4F6457", "#ACD0C0", "#0F1B07", "#F7EFE2", "#D09683", "#F62A00", "#A1BE95", 
                        "#20948B", "#9B4F0F", "#CB0000")

# 4. produce data frame with families' abundances per region
families.region <- BB22.full %>%
  dplyr::group_by(family.agg, region, bbspecies) %>%
  dplyr::summarise(cum.abund = sum(Abundance))
families.region$family.agg <- as.factor(families.region$family.agg)

### FIGURE S5 ####
# 5. plot families per region
setwd(output) 
ggplot(families.region, aes(fill = family.agg, y = cum.abund, x = region)) + 
  geom_bar(position="fill", stat="identity") + 
  theme_classic(base_size = 20) +
  facet_wrap(~bbspecies) + 
  ggtitle("Plant families per region") +
  labs(fill = 'Plant families', 
       x = "regions", 
       y = "relative abundance") +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values = palette.fams, 
                    name = "Plant families")
ggsave(paste("./diet pattern/pollen composition/PlantFamilies_per_Region.png", sep = ""), 
       width = 16, height = 8, device = "png", )
setwd(input) 

### origin status per site, region and species ----
# 1. produce data frame with exotic/native per site
ex.nat.site <- BB22.full %>%
  dplyr::group_by(native_exotic, site, bbspecies) %>%
  dplyr::summarise(abund = sum(binom.abund))
ex.nat.site$native_exotic <- as.factor(ex.nat.site$native_exotic)

# 2. plot native/exotic per site
# create color palette
palette.ex =c("#518A45", "#BDA509")

# 3. produce data frame with exotic/native per region
ex.nat.region <- BB22.full %>%
  dplyr::group_by(native_exotic, region, bbspecies) %>%
  dplyr::summarise(abund = sum(binom.abund))
ex.nat.region$native_exotic <- as.factor(ex.nat.region$native_exotic)

### FIGURE S6c ####
# 4. plot native/exotic per region    
setwd(output)
ggplot(ex.nat.region, aes(fill = native_exotic, y = abund, x = region)) + 
  geom_bar(position="fill", stat="identity") + 
  theme_classic(base_size = 20) +
  facet_wrap(~bbspecies) + 
  ggtitle("Origin status per region LEG") + 
  ylab("Proportion of species in the pollen") +
  labs(fill='Origin status') +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values = palette.ex, 
                    labels = c('Exotic', 'Native'), 
                    name = "")
ggsave(paste("./diet pattern/pollen composition/OriginStatus_per_Region.png", sep = ""), 
       width = 16, height = 8)
setwd(input)

### growth form per site, region and species ----
# 1. produce data frame with growth form per site
growth.site <- BB22.full %>%
  dplyr::group_by(growth_form_category, site, bbspecies) %>%
  dplyr::summarise(abund = sum(binom.abund))
growth.site$growth_form_category <- as.factor(growth.site$growth_form_category)

# 2. plot growth form per site
# create color palette
palette.growth =c("#375E97", "#518A45", "#FA812F", "#F34A4A", "#5A99AD", "#A3A3A3")

# 3. produce data frame with growth form per region
growth.region <- BB22.full %>%
  dplyr::group_by(growth_form_category, region, bbspecies) %>%
  dplyr::summarise(abund = sum(binom.abund))
growth.region$growth_form_category <- as.factor(growth.region$growth_form_category)

### FIGURE S6b ####
# 4. plot growth form per region   
setwd(output)
ggplot(growth.region, aes(fill = growth_form_category, y = abund, x = region)) + 
  geom_bar(position = "fill", stat = "identity") + 
  theme_classic(base_size = 20) +
  facet_wrap(~bbspecies) + 
  ggtitle("Growth Form LEG") + 
  ylab("Proportion of species in the pollen") + 
  labs(fill='Growth form') +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values = palette.growth, name = "")
ggsave(paste("./diet pattern/pollen composition/GrowthForm_per_Region.png", sep = ""), 
       width = 16, height = 8)
setwd(input)

### blossom class per site, region and species ----
# 1. produce data frame with blossom class per site
blossom.site <- BB22.full %>%
  dplyr::group_by(structural_blossom_class, site, bbspecies) %>%
  dplyr::summarise(abund = sum(binom.abund))
blossom.site$structural_blossom_class <- as.factor(blossom.site$structural_blossom_class)

# 2. plot blossom class per site
# create color palette
palette.bloss=c("#375E97", "#80BD9E", "#FA812F", "#758A30", "#07575B", "#D95F56", "#C4DFE6")

# 3. produce data frame with blossom class per region
blossom.region <- BB22.full %>%
  dplyr::group_by(structural_blossom_class, region, bbspecies) %>%
  dplyr::summarise(abund = sum(binom.abund))
blossom.region$structural_blossom_class <- as.factor(blossom.region$structural_blossom_class)

### FIGURE S6a ####
# 4. plot blossom class per region
setwd(output)
ggplot(blossom.region, aes(fill = structural_blossom_class, y = abund, x = region)) + 
  geom_bar(position = "fill", stat = "identity") + 
  theme_classic(base_size = 20) +
  facet_wrap(~bbspecies) + 
  ggtitle("Blossom Class per site LEG") + 
  ylab("Proportion of species in the pollen") +
  labs(fill='Blossom Class') +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values = palette.bloss, 
                    labels = c('Bell Trumpet', 'Brush', "Dish Bowl",
                             "Flag", "Gullet", "Stalk Disk", "Tube"),
                    name = "")
ggsave(paste("./diet pattern/pollen composition/BlossomClass_per_Region.png", sep = ""), 
       width = 16, height = 8)
setwd(input)

# TEST FAMILY, FLOWER, GRWOTH FORM, EXOTIC/NATIVE PROPORTIONS URBAN vs. RURAL ----
## preparation ----
# clear work environment
rm(list = ls()) 

#load libraries
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(insight)

# set working directory to main repository
input <- "~/Library/CloudStorage/GoogleDrive-simo1996s@gmail.com/My Drive/ETH/Master Thesis/bumblebee_try/01_DATA"
output <- "~/Library/CloudStorage/GoogleDrive-simo1996s@gmail.com/My Drive/ETH/Master Thesis/bumblebee_try/03_OUTPUT"

setwd(input)

# import data (added row with information on family)
BB22.full.family <- read_csv("BB22.full.family")

## analysis ----

BB.pasc <- BB22.full.family %>%
  filter(bbspecies == "B.pascuorum")
BB.lapi <- BB22.full.family %>%
  filter(bbspecies == "B.lapidarius")

# define vectors to loop over
resolution <- c("landscape", "region", "site")
fun.trait <- c("family.agg", "native_exotic", "growth_form_category", "structural_blossom_class")

### B.pascuorum ----
chisq.summary <- NULL
for (i in resolution) {
  for (j in fun.trait) {
    chisq <- chisq.test(BB.pasc[[i]], BB.pasc[[j]])
    chisq
    summary <- c(paste(i, j, sep = "~"), round(chisq$statistic, 3), chisq$parameter, chisq$p.value)
    chisq.summary <- rbind(chisq.summary, summary)
  } # end loop j
} # end loop i

chisq.summary <- as.data.frame(chisq.summary)
colnames(chisq.summary) <- c("formula", "X-squared", "df", "p")
rownames(chisq.summary) <- NULL
chisq.summary <- chisq.summary %>%
  mutate(p = format_p(p, stars = TRUE)) %>%
  format_table()

# export the summary of the tests
setwd(output)
write.table(chisq.summary, file = "diet pattern/chisq tests/chisq_Bpascuorum.txt", 
            sep = "\t", row.names = TRUE, col.names = NA)
setwd(input)

### B.lapidarius ----
chisq.summary.1 <- NULL
for (i in resolution) {
  for (j in fun.trait) {
    chisq <- chisq.test(BB.lapi[[i]], BB.lapi[[j]])
    chisq
    summary <- c(paste(i, j, sep ="~"), round(chisq$statistic, 3),chisq$parameter,chisq$p.value)
    chisq.summary.1 <- rbind(chisq.summary.1, summary)
  } # end loop j
} # end loop i

chisq.summary.1 <- as.data.frame(chisq.summary.1)
colnames(chisq.summary.1) <- c("formula", "X-squared", "df", "p")
rownames(chisq.summary.1) <- NULL
chisq.summary.1 <- chisq.summary.1 %>%
  mutate(p = format_p(p, stars = TRUE)) %>%
  format_table()

# export the summary of the tests
setwd(output)
write.table(chisq.summary.1, file = "diet pattern/chisq tests/chisq_Blapidarius.txt", 
            sep = "\t", row.names = TRUE, col.names = NA)
setwd(input)

# DOCUMENTED PLANT SPECIES PER SITE AND FOUND IN POLLEN ----

## preparation ----
# clear work environment
rm(list = ls()) 

#load libraries
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)

# set working directory to main repository
input <- "~/Library/CloudStorage/GoogleDrive-simo1996s@gmail.com/My Drive/ETH/Master Thesis/bumblebee_try/01_DATA"
output <- "~/Library/CloudStorage/GoogleDrive-simo1996s@gmail.com/My Drive/ETH/Master Thesis/bumblebee_try/03_OUTPUT"

setwd(input)

## analysis ----

# import species lists form GBIF and InfoFlora (file data preparation)
site.list.occ <- readRDS("sp_list_gbif_infoflora.RData")

# load data on pollen (file file data preparation)
site.list.bb <- readRDS("site_list_bb.RData")

# combine occurrence data with species found in pollen per site
site.list <- list()
sitenames <- c("ZHUA", "ZHUB", "ZHUC", "ZHRD", "ZHRE", "ZHRF", "BEUA", "BEUB", 
               "BEUC", "BERD", "BSUA", "BSUB", "BSUC", "BSRD", "BSRE", "BSRF")

for (i in sitenames) {
  x <- as.factor(site.list.occ [[i]])
  y <- as.factor(site.list.bb [[i]])
  site.list[[i]] <- unique(c(x,y)) %>%
    droplevels()
}

# calculate mean of number of species per site
mean.landsacpe <- c()
for (i in sitenames) {
  temp <- c(substring(i, 3, 3), nlevels(site.list[[i]]))
  mean.landsacpe <- rbind(mean.landsacpe, temp)
}

# adapt dataframe to be used with t-test
colnames(mean.landsacpe) <- c("landscape", "species_richness")
mean.landsacpe <- as.data.frame(mean.landsacpe)
mean.landsacpe$species_richness <- as.numeric(mean.landsacpe$species_richness)
str(mean.landsacpe)

# perform t-test and display it in a boxplot
library(rstatix)
w.test <- wilcox_test(mean.landsacpe, species_richness ~ landscape, paired = F)
palette.landscape <- c("#E69F00", "#56B4E9") #create color palette for landscape
a <- ggplot(mean.landsacpe, aes(x = landscape, y = species_richness,  fill = landscape)) +
  geom_boxplot(notch = T) +
  ylab("plant diversity of sites") +
  xlab("") +
  scale_x_discrete(labels = c('rural', 'urban')) +
  theme_classic(base_size = 20) +
  theme(aspect.ratio = 1) +
  scale_fill_manual(values = palette.landscape, 
                    guide = "none") +
  labs(subtitle = paste("W = ", w.test$statistic, ", p = ", w.test$p, sep = ""))
a # print the plot

# compile species richness per site in dataframe
rich.occ <- c()
for (i in sitenames) {
  x <- nlevels(site.list.occ[[i]])
  rich.occ <- c(rich.occ, x)
}

rich.bb <- c()
for (i in sitenames) {
  x <- nlevels(site.list.bb[[i]])
  rich.bb <- c(rich.bb, x)
}

df.site <- data_frame(site = sitenames,
                      landscape = rep_len(c(rep("U", 3), rep("R", 3)), length(sitenames)),
                      rich.occ = rich.occ,
                      rich.bb = rich.bb)

# plot the relationship of species richness of data bases and species collected by the bumblebees
ggplot(df.site, aes(x = rich.occ, y = rich.bb)) +
  geom_point()

# test the relationship with a model
fit <- lm(rich.occ ~ rich.bb, df.site)
summary(fit)

# plot it with information on the model
b <- ggplot(df.site, aes(x = rich.bb , y = rich.occ)) + 
  geom_point(aes(color = landscape), size = 3) + 
  labs(x ="species richness in pollen" , 
       y = "plant diversity of sites") +
  theme_classic(base_size = 20) + 
  theme(aspect.ratio = 1) + 
  geom_smooth(method = "lm", 
              se = FALSE, 
              col = "black") +
  scale_color_manual(values = palette.landscape) +
  stat_cor(
    aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), 
    label.x = 3,
    size = 7); b

# arrange them into one file to export
### FIGURE 3 ----
setwd(output)
ggarrange(a, b, ncol = 2, nrow = 1,
          labels = c("A", "B"))
ggsave("./diet pattern/plant richness landscape/sp_rich_occ_bb.png", 
       width = 20, height = 10)
setwd(input)


# DIET SIMILARITIES ----

## preparation ----
# clear work environment
rm(list = ls()) 

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
BB22.full <- read_csv("BB22_full.csv") %>%
  mutate(site = as_factor(site),
         species = as_factor(species),
         region = substr(site, 1, 3)) 

#remove information on leg and body from ID
BB22.full$ID.short = as.factor(substring(BB22.full$ID,1, nchar(BB22.full$ID)-1)) 

# sites' names for looping
sitenames <- c("ZHUA", "ZHUB", "ZHUC", "ZHRD", "ZHRE", "ZHRF", "BEUA", "BEUB", 
               "BEUC", "BERD", "BSUA", "BSUB", "BSUC", "BSRD", "BSRE", "BSRF")

## B.pascuorum ----
# only use B.pascuorum
BB22.full.species <- BB22.full[BB22.full$bbspecies == "B.pascuorum",]

# prepare list element to store correlation matrices for species, genus and family
pasc.ID <- list()
pasc.ID.mean <- list()

### calculate distance ----
for (i in sitenames) {
  # select a site
  # i <- "BSRE"
  BB22.full.species.site <- BB22.full.species[BB22.full.species$site == i,] %>%
    droplevels()
  
  BB22_full.ab <- BB22.full.species.site %>%
    group_by(ID.short, family) %>%
    summarise(abundance = sum(Abundance)) %>%
    distinct() # remove duplicates
  
  BB22_full.ab.new <- c()
  for (h in unique(BB22_full.ab$ID.short)) {
    temp <- BB22_full.ab[BB22_full.ab$ID.short == h,]
    perc <- sum(temp$abundance)
    for (k in 1:nrow(temp)) {
      temp$ab.new[k] <- 100 / perc * temp$abundance[k]/100
    }
    BB22_full.ab.new <- rbind(BB22_full.ab.new, temp[, c(1,2,4)])
  }
  
  # ON FAMILY
  # convert into matrix
  library(reshape2)
  BB22.full.table <- dcast(BB22_full.ab.new, ID.short ~ family, value.var = "ab.new")
  rownames(BB22.full.table) <- BB22.full.table$ID.short
  BB22.full.table[is.na(BB22.full.table)] <- 0
  BB22.full.table <- BB22.full.table[, -1]
  
  # compute measures of distance (or resemblance) and store them in list
  library(vegan)
  dist <- vegdist(BB22.full.table,method = "bray")
  pasc.ID[[i]] <- as.numeric(dist[1:1000])
  pasc.ID.mean[[i]] <- mean(dist)
}


### plotting ----
# transform data
pasc.ID.df <- do.call(cbind, pasc.ID)
pasc.ID.df <- melt(pasc.ID.df)[, -1]
colnames(pasc.ID.df) <- c("site", "dist")
pasc.ID.df$landscape <- substr(pasc.ID.df$site, 3, 3)

palette.landscape <- c("#E69F00", "#56B4E9") #create color palette for landscape

#### FIGURE S7 ----
ggplot(pasc.ID.df, aes(x = site, y = dist, fill = landscape)) +
  geom_boxplot(notch = T) +
  xlab("") +
  ylab("distance") +
  ggtitle("B.pascuorum: distance between individuals per site") +
  theme_classic(base_size = 20) +
  theme(aspect.ratio = 1) +
  scale_fill_manual(values = palette.landscape, guide = "none") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

setwd(output)
ggsave("./diet pattern/distances/pasc_dist_family.png", 
       width = 10, height = 10)
setwd(input)

# add mean and SD of distance to a data frame and export it
pasc.ID.df.site <- pasc.ID.df %>%
  group_by(site) %>%
  summarise(mean.pasc = mean(dist, na.rm = TRUE), 
            sd.pasc = sd(dist, na.rm = TRUE))

setwd(output)
write.csv(pasc.ID.df.site, "./diet pattern/distances/statistics_pasc.csv")
setwd(input)

## B.lapidarius ----
# only use B.B.lapidarius
BB22.full.species <- BB22.full[BB22.full$bbspecies == "B.lapidarius",]

# prepare list element to store correlation matrices for species, genus and family
lapi.ID <- list()
lapi.ID.mean <- list()

### calculate distance ----
for (i in sitenames) {
  BB22.full.species.site <- BB22.full.species[BB22.full.species$site == i,] %>%
    droplevels()
  
  BB22_full.ab <- BB22.full.species.site %>%
    group_by(ID.short, family) %>%
    summarise(abundance = sum(Abundance)) %>%
    distinct() # remove duplicates
  
  BB22_full.ab.new <- c()
  for (h in unique(BB22_full.ab$ID.short)) {
    temp <- BB22_full.ab[BB22_full.ab$ID.short == h,]
    perc <- sum(temp$abundance)
    for (k in 1:nrow(temp)) {
      temp$ab.new[k] <- 100 / perc * temp$abundance[k] / 100
    }
    BB22_full.ab.new <- rbind(BB22_full.ab.new, temp[, c(1,2,4)])
  }
  
  # ON FAMILY
  # convert into matrix
  library(reshape2)
  BB22.full.table <- dcast(BB22_full.ab.new, ID.short ~ family, value.var = "ab.new")
  rownames(BB22.full.table) <- BB22.full.table$ID.short
  BB22.full.table[is.na(BB22.full.table)] <- 0
  BB22.full.table <- BB22.full.table[, -1]
  
  # compute measures of distance (or resemblance) and store them in list
  library(vegan)
  dist <- vegdist(BB22.full.table,method = "bray")
  lapi.ID[[i]] <- as.numeric(dist[1:1000])
  lapi.ID.mean[[i]] <- mean(dist)
}

### plotting ----
# transform data
lapi.ID.df <- do.call(cbind, lapi.ID)
lapi.ID.df <- melt(lapi.ID.df)[, -1]
colnames(lapi.ID.df) <- c("site", "dist")
lapi.ID.df$landscape <- substr(lapi.ID.df$site, 3, 3)

palette.landscape <- c("#E69F00", "#56B4E9") #create color palette for landscape

#### FIGURE S7 ----
ggplot(lapi.ID.df, aes(x = site, y = dist, fill = landscape)) +
  geom_boxplot(notch = T) +
  xlab("") +
  ylab("distance") +
  ggtitle("B.lapidarius: distance between individuals per site") +
  theme_classic(base_size = 20) +
  theme(aspect.ratio = 1) +
  scale_fill_manual(values = palette.landscape, 
                    guide = "none") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

setwd(output)
ggsave("./diet pattern/distances/lapi_dist_family.png", 
       width = 10, height = 10)
setwd(input)

# add mean and SD of distance to a dataframe and export it
lapi.ID.df.site <- lapi.ID.df %>%
  group_by(site) %>%
  summarise(mean.lapi = mean(dist, na.rm = TRUE), 
            sd.lapi = sd(dist, na.rm = TRUE))

setwd(output)
write.csv(lapi.ID.df.site, "./diet pattern/distances/statistics_lapi.csv")
setwd(input)

# compare species summary statistics ----
library(stats)
library(tidyverse)
library(rstatix)
library(ggpubr)

# compile sd and mean for both species in one dataframe and export it
site.summary <- data.frame(site = rep(lapi.ID.df.site$site, 2))
site.summary <- site.summary %>%
  mutate(bbspecies = c(rep("B.lapidarius", 16), rep("B.pascuorum", 16)),
         landscape = as_factor(substring(site, 3, 3)),
         mean = c(lapi.ID.df.site$mean.lapi, pasc.ID.df.site$mean.pasc),
         sd = c(lapi.ID.df.site$sd.lapi, pasc.ID.df.site$sd.pasc)) %>%
  mutate(landscape = fct_relevel(landscape, c("R", "U")))

setwd(output)
write.csv(site.summary, "./diet pattern/distances/distance_summary_statistics.csv")
setwd(input)

# subset for different comparisons
site.summary.lapi <- site.summary %>% filter(bbspecies == "B.lapidarius")
site.summary.pasc <- site.summary %>% filter(bbspecies == "B.pascuorum")
site.summary.urban <- site.summary %>% filter(landscape == "U")
site.summary.rural <- site.summary %>% filter(landscape == "R")

# prepare list for plots
plot.list <- list()

### compare the mean of landscape for B.lapidarius ----
res.aov1 <- site.summary.lapi %>% anova_test(mean ~ landscape)
res.aov1
pwc1 <- site.summary.lapi %>%
  pairwise_t_test(mean ~ landscape, p.adjust.method = "bonferroni")
pwc1

# Show adjusted p-values
pwc1 <- pwc1 %>% add_xy_position(x = "landscape")

#### FIGURE S8a ----
plot.list[[1]] <-
  ggboxplot(site.summary.lapi, x = "landscape", y = "mean", fill = "landscape", notch = TRUE) +
  stat_pvalue_manual(pwc1, label = "p.adj", tip.length = 0, step.increase = 0.1) +
  labs(subtitle = get_test_label(res.aov1, detailed = TRUE)) +
  ggtitle("B.lapidarius") + # for the main title
  xlab("landscape") +
  ylab("Mean") +
  scale_x_discrete(labels = c('rural', 'urban')) +
  theme_classic(base_size = 20) +
  theme(aspect.ratio = 1) +
  scale_fill_manual(values = c("#E69F00", "#56B4E9"))

### compare the mean of landscape for B.pascuorum ----
res.aov2 <- site.summary.pasc %>% anova_test(mean ~ landscape)
res.aov2
pwc2 <- site.summary.pasc %>%
  pairwise_t_test(mean ~ landscape, p.adjust.method = "bonferroni")
pwc2

# Show adjusted p-values
pwc2 <- pwc2 %>% add_xy_position(x = "landscape")

#### FIGURE S8b ----
plot.list[[2]] <-
  ggboxplot(site.summary.pasc, x = "landscape", y = "mean", fill = "landscape", notch = TRUE) +
  stat_pvalue_manual(pwc2, label = "p.adj", tip.length = 0, step.increase = 0.1) +
  labs(subtitle = get_test_label(res.aov2, detailed = TRUE)) +
  ggtitle("B.pascuorum") + # for the main title
  xlab("landscape") +
  ylab("Mean") +
  scale_x_discrete(labels = c('rural', 'urban')) +
  theme_classic(base_size = 20) +
  theme(aspect.ratio = 1) +
  scale_fill_manual(values = c("#E69F00", "#56B4E9"))

### compare the sd of landscape for B.lapidarius ----
res.aov5 <- site.summary.lapi %>% anova_test(sd ~ landscape)
res.aov5
pwc5 <- site.summary.lapi %>%
  pairwise_t_test(sd ~ landscape, p.adjust.method = "bonferroni")
pwc5

# Show adjusted p-values
pwc5 <- pwc5 %>% add_xy_position(x = "landscape")

#### FIGURE S8c ----
plot.list[[3]] <-
  ggboxplot(site.summary.lapi, x = "landscape", y = "sd", fill = "landscape", notch = TRUE) +
  stat_pvalue_manual(pwc5, label = "p.adj", tip.length = 0, step.increase = 0.1) +
  labs(subtitle = get_test_label(res.aov5, detailed = TRUE)) +
  ggtitle("B.lapidarius") + # for the main title
  xlab("landscape") +
  ylab("Standard deviation") +
  scale_x_discrete(labels = c('rural', 'urban')) +
  theme_classic(base_size = 20) +
  theme(aspect.ratio = 1) +
  scale_fill_manual(values = c("#E69F00", "#56B4E9"))
### compare the sd of landscape for B.pascuorum ----
res.aov6 <- site.summary.pasc %>% anova_test(sd ~ landscape)
res.aov6
pwc6 <- site.summary.pasc %>%
  pairwise_t_test(sd ~ landscape, p.adjust.method = "bonferroni")
pwc6

# Show adjusted p-values
pwc6 <- pwc6 %>% add_xy_position(x = "landscape")

#### FIGURE S8d ----
plot.list[[4]] <-
  ggboxplot(site.summary.pasc, x = "landscape", y = "sd", fill = "landscape", notch = TRUE) +
  stat_pvalue_manual(pwc6, label = "p.adj", tip.length = 0, step.increase = 0.1) +
  labs(subtitle = get_test_label(res.aov6, detailed = TRUE)) +
  ggtitle("B.pascuorum") + # for the main title
  xlab("landscape") +
  ylab("Standard deviation") +
  scale_x_discrete(labels = c('rural', 'urban')) +
  theme_classic(base_size = 20) +
  theme(aspect.ratio = 1) +
  scale_fill_manual(values = c("#E69F00", "#56B4E9")) 

### plotting ----
#### FIGURE S8 ----

# arrange them into one file to export
setwd(output)
plot <- ggarrange(plot.list[[1]], plot.list[[2]],
                  plot.list[[3]], plot.list[[4]],
                  ncol = 2, nrow = 2,
                  labels = c("A", "B", "C", "D"))
ggsave("./diet pattern/distances/summary_stats_mean_overview.png", 
       width = 20, height = 20)
setwd(input)

### compare the mean of species in urban areas ----
res.aov3 <- site.summary.urban %>% anova_test(mean ~ bbspecies)
res.aov3
pwc3 <- site.summary.urban %>%
  pairwise_t_test(mean ~ bbspecies, p.adjust.method = "bonferroni")
pwc3

# Show adjusted p-values
pwc3 <- pwc3 %>% add_xy_position(x = "bbspecies")

#### FIGURE S9a ----
plot.list[[5]] <-
  ggboxplot(site.summary.urban, x = "bbspecies", y = "mean", fill = "bbspecies", notch = TRUE) +
  stat_pvalue_manual(pwc3, label = "p.adj", tip.length = 0, step.increase = 0.1) +
  labs(subtitle = get_test_label(res.aov3, detailed = TRUE)) +
  ggtitle("Urban areas") + # for the main title
  xlab("Bumblebee species") +
  ylab("Mean") +
  scale_x_discrete(labels = c('B.lapidarius', 'B.pascuorum')) +
  theme_classic(base_size = 20) +
  theme(aspect.ratio = 1) +
  scale_fill_manual(values = c("#291600", "#e0b802"))

### compare the mean of species in rural areas ----
res.aov4 <- site.summary.rural %>% anova_test(mean ~ bbspecies)
res.aov4
pwc4 <- site.summary.rural %>%
  pairwise_t_test(mean ~ bbspecies, p.adjust.method = "bonferroni")
pwc4

# Show adjusted p-values
pwc4 <- pwc4 %>% add_xy_position(x = "bbspecies")

#### FIGURE S9b ----
plot.list[[6]] <-
  ggboxplot(site.summary.rural, x = "bbspecies", y = "mean", fill = "bbspecies", notch = TRUE) +
  stat_pvalue_manual(pwc4, label = "p.adj", tip.length = 0, step.increase = 0.1) +
  labs(subtitle = get_test_label(res.aov4, detailed = TRUE)) +
  ggtitle("Rural areas") + # for the main title
  xlab("Bumblebee species") +
  ylab("Mean") +
  scale_x_discrete(labels = c('B.lapidarius', 'B.pascuorum')) +
  theme_classic(base_size = 20) +
  theme(aspect.ratio =1 ) +
  scale_fill_manual(values = c("#291600","#e0b802")) 

### compare the sd of species in urban areas ----
res.aov7 <- site.summary.urban %>% anova_test(sd ~ bbspecies)
res.aov7
pwc7 <- site.summary.urban %>%
  pairwise_t_test(sd ~ bbspecies, p.adjust.method = "bonferroni")
pwc7

# Show adjusted p-values
pwc7 <- pwc7 %>% add_xy_position(x = "bbspecies")
#### FIGURE S9c ----
plot.list[[7]] <-
  ggboxplot(site.summary.urban, x = "bbspecies", y = "sd", fill = "bbspecies", notch = TRUE) +
  stat_pvalue_manual(pwc7, label = "p.adj", tip.length = 0, step.increase = 0.1) +
  labs(subtitle = get_test_label(res.aov7, detailed = TRUE)) +
  ggtitle("Urban areas") + # for the main title
  xlab("Bumblebee species") +
  ylab("Standard deviation") +
  scale_x_discrete(labels = c('B.lapidarius', 'B.pascuorum')) +
  theme_classic(base_size = 20) +
  theme(aspect.ratio = 1) +
  scale_fill_manual(values = c("#291600","#e0b802"))

### compare the sd of species in rural areas ----
res.aov8 <- site.summary.rural %>% anova_test(sd ~ bbspecies)
res.aov8
pwc8 <- site.summary.rural %>%
  pairwise_t_test(sd ~ bbspecies, p.adjust.method = "bonferroni")
pwc8

# Show adjusted p-values
pwc8 <- pwc8 %>% add_xy_position(x = "bbspecies")

#### FIGURE S9d ----
plot.list[[8]] <-
  ggboxplot(site.summary.rural, x = "bbspecies", y = "sd", fill = "bbspecies", notch = TRUE) +
  stat_pvalue_manual(pwc8, label = "p.adj", tip.length = 0, step.increase = 0.1) +
  labs(subtitle = get_test_label(res.aov8, detailed = TRUE)) +
  ggtitle("Rural areas") + # for the main title
  xlab("Bumblebee species") +
  ylab("Standard deviation") +
  scale_x_discrete(labels = c('B.lapidarius', 'B.pascuorum')) +
  theme_classic(base_size = 20) +
  theme(aspect.ratio = 1) +
  scale_fill_manual(values = c("#291600","#e0b802")) 


### plotting ----
#### FIGURE S9 ----

# arrange them into one file to export
setwd(output)
plot <- ggarrange(plot.list[[5]], plot.list[[6]],
                  plot.list[[7]], plot.list[[8]],
                  ncol = 2, nrow = 2,
                  labels = c("A", "B", "C", "D"))
ggsave("./diet pattern/distances/summary_stats_sd_overview.png", 
       width = 20, height = 20)
setwd(input)

# FUNCTIONAL DIVERSITY METRICS BETWEEN SPECIES ----

## preparation ----
# clear work environment
rm(list = ls()) 

#load libraries
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(glmmTMB)
library(DHARMa)

# set working directory to main repository
input <- "~/Library/CloudStorage/GoogleDrive-simo1996s@gmail.com/My Drive/ETH/Master Thesis/bumblebee_try/01_DATA"
output <- "~/Library/CloudStorage/GoogleDrive-simo1996s@gmail.com/My Drive/ETH/Master Thesis/bumblebee_try/03_OUTPUT"

# load data
setwd(input)
BB22.bb.traits <- read_csv("BB22_traits.csv")

## combined species ----

### data preparation ----

# load data
setwd(input)
BB22.bb.traits <- read_csv("BB22_traits.csv")

# replace the names for urban and rural
for (i in 1:nrow(BB22.bb.traits)) {
  if (BB22.bb.traits$landscape[i] == "urban") {
    BB22.bb.traits$landscape[i] = "U"
  } else {
    BB22.bb.traits$landscape[i] = "R"
  }
}

# rename column site to match other dataframes
BB22.bb.traits <- BB22.bb.traits %>%
  mutate(site = paste(location, landscape, replicate, sep = ""),
         region = paste(location, landscape, sep = "")) 

# import data with spatial information on sites
BB22.sites.meta <- read_csv("BB22_sites_2016.csv")

# import data with functional diversity of plants per site
BB22.fun.ID.pasc <- read_csv("./FD/FD_B.pascuorum_ID.short.csv") %>%
  rename_with(.cols = 1, ~"ID") %>%
  rename_with(.cols = 2, ~"sp_richn") %>%
  rename_with(.cols = 3, ~"fric") %>%
  rename_with(.cols = 5, ~"fdiv") %>%
  rename_with(.cols = 4, ~"feve")
BB22.fun.ID.pasc <- BB22.fun.ID.pasc %>%
  mutate(species = rep("B.pascuorum", length(BB22.fun.ID.pasc$ID)))

BB22.fun.ID.lapi <- read_csv("./FD/FD_B.lapidarius_ID.short.csv") %>%
  rename_with(.cols = 1, ~"ID") %>%
  rename_with(.cols = 2, ~"sp_richn") %>%
  rename_with(.cols = 3, ~"fric") %>%
  rename_with(.cols = 5, ~"fdiv") %>%
  rename_with(.cols = 4, ~"feve")
BB22.fun.ID.lapi <- BB22.fun.ID.lapi %>%
  mutate(species = rep("B.lapidarius", length(BB22.fun.ID.lapi$ID)))

### analysis ----
# combine both species
BB22.fun.ID <- rbind(BB22.fun.ID.pasc, BB22.fun.ID.lapi)

# filter functional metrics of plant traits to compare with traits per individual
BB22.fun.ID <- subset(BB22.fun.ID, ID %in% BB22.bb.traits$ID)

# add site coordinates to the trait data frame (in LV95)
BB22.ID <- merge(BB22.bb.traits, BB22.fun.ID, by  = "ID", all.x=TRUE)
BB22.ID <- merge(BB22.ID, BB22.sites.meta[, c(1,2,3)], by  = "site", all.x=TRUE) 

BB22.ID$bbspecies <- as.factor(BB22.ID$bbspecies)
levels(BB22.ID$bbspecies)

# plot the relationship of plants traits of one site and bumblebee traits of one site
traits <- colnames(BB22.ID[, 7:15]) # bumblebee traits to look at; prepare for loop
metrics <- colnames(BB22.ID[, c(17:20)]) # plant FD to look at (see file BB22_compute_FD); prepare for loop

# wilcox test
library(rstatix)
plot_list <- list()
w.test.list <- list()
mean.lapi <- list()
mean.pasc <- list()

for (j in metrics) {
  gg.data <- data.frame(species = BB22.ID$bbspecies, 
                        value = BB22.ID[[j]])
  w.test <- wilcox_test(gg.data,value~species)
  w.test.list[[j]] <- w.test
  mean.lapi[[j]] <- mean(gg.data$value[gg.data$species == "B.lapidarius"], na.rm=TRUE)
  mean.pasc[[j]] <- mean(gg.data$value[gg.data$species == "B.pascuorum"], na.rm=TRUE)
  
  p <- ggplot(gg.data, aes(x = species, y = value, fill = species)) + 
    geom_boxplot(notch = T) +
    xlab("") +
    scale_x_discrete(labels = c("B.lapidarius", "B.pascuroum")) +
    theme_classic(base_size = 20) +
    theme(aspect.ratio = 1) +
    guides(alpha = "none") +
    scale_fill_manual(values = c("#291600","#e0b802"), 
                      guide = "none") + 
    labs(subtitle = paste("p = ", w.test$p, sep = ""))
  if (j == "sp_richn") {
    p <- p + ylim(0, 18) + ylab("species richness")
  } else if (j == "fric") {
    p <- p + ylim(0, 1) + ylab("funtional richness")
  } else if (j == "fdiv") {
    p <- p + ylim(0.6, 0.85) + ylab("funtional divergence")
  } else {
    p <- p + ylim(0.5, 0.75) + ylab("funtional evenness")
  }
  plot_list[[j]] <- p
}

#### FIGURE 2 ----
# arrange them into one file to export
setwd(output)
plot <- ggarrange(plot_list[[1]],plot_list[[2]],
                  plot_list[[3]],plot_list[[4]],
                  ncol = 1, nrow = 4,
                  labels = c("A", "B", "C", "D"))
annotate_figure(plot, top = text_grob("species richness and funtional diversity between species", 
                                      face = "bold", size = 22))
ggsave("./diet pattern/FD/FD_box_BBtraits_species.png", 
       width = 6, height = 24)
setwd(input)

# export statistics of W tests
w.stats <- rbind(c(mean.lapi$sp_richn, mean.pasc$sp_richn, w.test.list$sp_richn),
                 c(mean.lapi$fric, mean.pasc$fric, w.test.list$fric),
                 c(mean.lapi$feve, mean.pasc$feve, w.test.list$feve),
                 c(mean.lapi$fdiv, mean.pasc$fdiv, w.test.list$fdiv))
setwd(output)
write.csv(w.stats,"./diet pattern/FD/FD_W_BBtraits_species.csv")
setwd(input)

# explore relationships
library(nlme)

#### FIGURES SUPPLEMENT SECTION "Bumblebee Morphology and Diet: Species" ----
# preparation for caption string
fmt <- "%s: adj.R^2 = %.3f, p = %.3f"

# perform loop to output plots per relationship summarized per FD
for (i in metrics) {
  x <- 1 # for naming the plots
  for (j in traits) {
    f <- formula(paste(i,"~", j))
    
    # fit lm for lapi and get R2 and p
    lm.lapi <- lm (f, data = BB22.ID[BB22.ID$bbspecies == "B.lapidarius",])
    sum.lapi <- summary(lm.lapi)
    lab.lapi <- sprintf(fmt, "B.lapidarius", sum.lapi$adj.r.squared, coef(sum.lapi)[2, 4])
    
    # fit lm for pasc and get R2 and p
    lm.pasc <- lm (f, data = BB22.ID[BB22.ID$bbspecies == "B.pascuorum",])
    sum.pasc <- summary(lm.pasc)
    lab.pasc <- sprintf(fmt, "B.pascuorum", sum.pasc$adj.r.squared, coef(sum.pasc)[2, 4])
    
    assign(paste("a", x, sep = ""), # assign the ggplot to plot name
           # define the ggplot
           ggplot(BB22.ID, aes_string(j, i, colour = "bbspecies")) + 
             geom_point() + 
             theme_classic(base_size = 20) + 
             theme(aspect.ratio = 1) + 
             labs(caption = paste(lab.lapi, "\n", lab.pasc)) +
             geom_smooth(method="lm", se = FALSE) +
             scale_color_manual(values = c("#291600","#e0b802"), labels = c("B.lapidarius", "B.pascuorum"))
           )
    x <- x+1
  } # end loop j
  
  setwd(output)
  plot4 <- ggarrange(a1, a2, a3, a4, a5, a6, a7, a8, a9, # arrange to plots nicely and export them
                     ncol = 3, nrow = 3,
                     labels = c(LETTERS[1:9]),
                     common.legend = TRUE)
  annotate_figure(plot4, top = text_grob(paste("comparison of ", i, " and traits between species", sep = ""),
                                         face = "bold", size = 22))
  ggsave(paste("./diet pattern/FD/FD_corr_", i, "_BBtraits_species.png", sep = ""),
         width = 15, height = 15)
  setwd(input)
} # end loop i


