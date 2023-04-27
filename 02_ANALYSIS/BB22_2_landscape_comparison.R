################################################
# Comparison between urban and rural bumblebee populations
# by Simonetta Selva
#
# Created: May 2nd, 2023
# Project: Bumblebee 2022
################################################

# HYPOTHESIS: The bumblebee's diet will change based on the availability of 
# plants in their landscape

# information: 
# 1) every subsection works in itself
# 2) all the plots not used in the main work are in the script but # are used to not print them

# FUNCTIONAL DIVERSITY METRICS ----
## preparation ----
# clear work environment
rm(list=ls()) 

#load libraries
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(glmmTMB)

# set working directory to main repository
input <- "YOUR PATH/Bumblebee_Simonetta_final/01_DATA" #change „YOUR PATH“ to where the repo is
output <- "YOUR PATH/Bumblebee_Simonetta_final/03_OUTPUT" #change „YOUR PATH“ to where the repo is
setwd(input)

# load data
BB22.bb.traits <- read_csv("BB22_traits.csv")

## B.PASCUORUM ---- 
# only B.pascuorum 
BB22.bb.traits.sp <- BB22.bb.traits[BB22.bb.traits$bbspecies == "B.pascuorum",]

# import data with spatial information on sites
BB22.sites.meta <- read_csv("BB22_sites_2016.csv")

# import data with functional diversity of plants per site
BB22.fun.ID <- read_csv("./FD/FD_B.pascuorum_ID.short.csv")%>% 
  rename_with(.cols = 1, ~"ID")%>%
  rename_with(.cols = 2, ~"sp_richn")%>%
  rename_with(.cols = 3, ~"fric")%>%
  rename_with(.cols = 5, ~"fdiv")%>%
  rename_with(.cols = 4, ~"feve")

# filter functional metrics of plant traits to compare with traits per individual
BB22.fun.ID <- subset(BB22.fun.ID, ID %in% BB22.bb.traits.sp$ID)

# add site coordinates to the trait data frame (in LV95)
BB22.ID <- merge(BB22.bb.traits.sp, BB22.fun.ID, 
                 by  = "ID", 
                 all.x=TRUE)
BB22.ID <- merge(BB22.ID, BB22.sites.meta[, c(1,2,3)], 
                 by  = "site", 
                 all.x=TRUE) 

### look at data ----

# Boxplots for all the variables we want to look at with Wilcoxon test
resp <- c("sp_richn", "fric", "fdiv", "feve")
library(psych)
library(rstatix)
plot_list  <- list()
w.test.list <- list()
mean.rural <- list()
mean.urban <- list()
palette.landscape <- c("#E69F00", "#56B4E9") #create color palette for landscape

for (j in resp) {
  gg.data <- data.frame(landscape=BB22.ID$landscape, value=BB22.ID[[j]])
  w.test <- wilcox_test(gg.data,value ~ landscape)
  w.test.list[[j]] <- w.test
  mean.rural[[j]] <- mean(gg.data$value[gg.data$landscape == "R"], 
                          na.rm=TRUE)
  mean.urban[[j]] <- mean(gg.data$value[gg.data$landscape == "U"], 
                          na.rm=TRUE)
  
  p <- ggplot(gg.data, aes(x=landscape, y = value, fill=landscape)) + 
    geom_boxplot(notch = T) + 
    xlab("") +
    scale_x_discrete(labels=c('rural', 'urban'))+
    theme_classic(base_size = 20) +     
    theme(aspect.ratio=1) + 
    guides(alpha = "none") +
    scale_fill_manual(values=palette.landscape, guide = "none") + 
    labs(subtitle = paste("p = ", w.test$p, sep=""))
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

#### FIGURE 4b ----
# arrange them into one file to export
setwd(output)
plot <- ggarrange(plot_list[[1]],plot_list[[2]],
                  plot_list[[3]],plot_list[[4]],
                  ncol = 1, nrow = 4,
                  labels = c("A", "B", "C", "D"))
annotate_figure(plot, top = text_grob("B.pascuorum: species richness and funtional diversity across landscapes", 
                                      face = "bold", size = 22))
ggsave("./landscape comparison/functional diversity/pasc_ID/FD_B.pascuorum_ID.png", width = 6, height = 24)
setwd(input)

# export statistics of W tests
w.stats.pasc <- rbind(c(mean.rural$sp_richn, mean.urban$sp_richn, w.test.list$sp_richn),
                      c(mean.rural$fric, mean.urban$fric, w.test.list$fric),
                      c(mean.rural$feve, mean.urban$feve, w.test.list$feve),
                      c(mean.rural$fdiv, mean.urban$fdiv, w.test.list$fdiv))
setwd(output)
write.csv(w.stats.pasc,"./landscape comparison/functional diversity/pasc_ID/FD_W_pasc_BBtraits_species.csv")
setwd(input)

# plot the relationship of plants traits of one site and bumblebee traits of one site
traits <- colnames(BB22.ID[, 7:15]) # bumblebee traits to look at; prepare for loop
metrics <- colnames(BB22.ID[, c(16:19)]) # plant FD to look at (see file BB22_compute_FD); prepare for loop

library(nlme)

#### FIGURES SUPPLEMENT SECTION "Bumblebee Morphology and Diet: Landscape" ----

# perform loop to output plots per relationship summarized per FD
# preparation for caption string
fmt <- "%s: adj.R^2 = %.3f, p = %.3f"

for (i in metrics) {
  x <- 1 # for naming the plots
  for (j in traits) {
    f <- formula(paste(i,"~", j))
    
    # fit lm for urban and get R2 and p
    lm.urban<- lm (f, data = BB22.ID[BB22.ID$landscape == "urban",])
    sum.urban <- summary(lm.urban)
    lab.urban <- sprintf(fmt, "urban", sum.urban$adj.r.squared, coef(sum.urban)[2, 4])
    
    # fit lm for rural and get R2 and p
    lm.rural <- lm (f, data = BB22.ID[BB22.ID$landscape == "rural",])
    sum.rural <- summary(lm.rural)
    lab.rural <- sprintf(fmt, "rural", sum.rural$adj.r.squared, coef(sum.rural)[2, 4])
    
    assign(paste("a", x, sep=""), # assign the ggplot to plot name
           # define the ggplot
           ggplot(BB22.ID, aes_string(j, i, colour = "landscape")) + 
             geom_point() + 
             theme_classic(base_size = 20) + 
             theme(aspect.ratio=1) + 
             labs(caption = paste(lab.urban, "\n", lab.rural)) +
             geom_smooth(method="lm", se = FALSE) +
             scale_color_manual(values=palette.landscape, labels=c("rural", "urban"))
    )
    x <- x+1
  } # end loop j
  setwd(output)
  plot4 <- ggarrange(a1,a2,a3,a4,a5,a6,a7,a8,a9, # arrange to plots nicely and export them 
                     ncol = 3, nrow = 3, 
                     labels = c(LETTERS[1:9]),   
                     common.legend = TRUE)
  annotate_figure(plot4, top = text_grob(paste("B.pascuorum: comparison of ", i, " and traits across landscapes", sep = ""),
                                         face = "bold", size = 22))
  ggsave(paste("./landscape comparison/functional diversity/pasc_ID/FD_pasc_corr_", i, "_BBtraits_landscapes.png", sep = ""), width = 15, height = 15)
  setwd(input)
} # end loop i

## B.LAPIDARIUS ---- 
# only B.lapidarius 
BB22.bb.traits.sp <- BB22.bb.traits[BB22.bb.traits$bbspecies == "B.lapidarius",]

# import data with spatial information on sites
BB22.sites.meta <- read_csv("BB22_sites_2016.csv")

# import data with functional diversity of plants per site
BB22.fun.ID <- read_csv("./FD/FD_B.lapidarius_ID.short.csv")%>% 
  rename_with(.cols = 1, ~"ID")%>%
  rename_with(.cols = 2, ~"sp_richn")%>%
  rename_with(.cols = 3, ~"fric")%>%
  rename_with(.cols = 5, ~"fdiv")%>%
  rename_with(.cols = 4, ~"feve")

# filter functional metrics of plant traits to compare with traits per individual
BB22.fun.ID <- subset(BB22.fun.ID, ID %in% BB22.bb.traits.sp$ID)

# add site coordinates to the trait data frame (in LV95)
BB22.ID <- merge(BB22.bb.traits.sp, BB22.fun.ID, by  = "ID", all.x=TRUE)
BB22.ID <- merge(BB22.ID, BB22.sites.meta[, c(1,2,3)], by  = "site", all.x=TRUE) 

### look at data ----

# Boxplots for all the variables we want to look at with Wilcoxon test
resp <- c("sp_richn", "fric", "fdiv", "feve")
library(psych)
library(rstatix)
plot_list  <- list()
w.test.list <-  list()
mean.rural <- list()
mean.urban <- list()
palette.landscape <- c("#E69F00", "#56B4E9") #create color palette for landscape

for (j in resp) {
  # j <- "sp_richn"
  gg.data <- data.frame(landscape=BB22.ID$landscape, value=BB22.ID[[j]])
  w.test <- wilcox_test(gg.data, value ~ landscape)
  w.test.list[[j]] <- w.test
  mean.rural[[j]] <- mean(gg.data$value[gg.data$landscape == "R"], na.rm=TRUE)
  mean.urban[[j]] <- mean(gg.data$value[gg.data$landscape == "U"], na.rm=TRUE)
  
  p <- ggplot(gg.data, aes(x=landscape, y = value, fill=landscape)) + 
    geom_boxplot(notch = T) + 
    xlab("") +
    scale_x_discrete(labels=c('rural', 'urban'))+
    theme_classic(base_size = 20) +     
    theme(aspect.ratio=1) + 
    guides(alpha = "none") +
    scale_fill_manual(values=palette.landscape, guide = "none") + 
    labs(subtitle = paste("p = ", w.test$p, sep=""))
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

#### FIGURE 4a ----
# arrange them into one file to export
setwd(output)
plot <- ggarrange(plot_list[[1]], plot_list[[2]],
                  plot_list[[3]], plot_list[[4]],
                  ncol = 1, nrow = 4,
                  labels = c("A", "B", "C", "D"))
annotate_figure(plot, top = text_grob("B.lapidarius: species richness and funtional diversity across landscapes", 
                                      face = "bold", size = 22))
ggsave("./landscape comparison/functional diversity/lapi_ID/FD_B.lapidarius_ID.png", width = 6, height = 24)
setwd(input)

# export statistics of W tests
w.stats.lapi <- rbind(c(mean.rural$sp_richn, mean.urban$sp_richn, w.test.list$sp_richn),
                      c(mean.rural$fric, mean.urban$fric, w.test.list$fric),
                      c(mean.rural$feve, mean.urban$feve, w.test.list$feve),
                      c(mean.rural$fdiv, mean.urban$fdiv, w.test.list$fdiv))
setwd(output)
write.csv(w.stats.lapi,"./landscape comparison/functional diversity/lapi_ID/FD_lapi_W_BBtraits_species.csv")
setwd(input)

#### FIGURES SUPPLEMENT SECTION "Bumblebee Morphology and Diet: Landscape" ----
# plot the relationship of plants traits of one site and bumblebee traits of one site
traits <- colnames(BB22.ID[, 7:15]) # bumblebee traits to look at; prepare for loop
metrics <- colnames(BB22.ID[, c(16:19)]) # plant FD to look at (see file BB22_compute_FD); prepare for loop

library(nlme)

# perform loop to output plots per relationship summarized per FD
# preparation for caption string
fmt <- "%s: adj.R^2 = %.3f, p = %.3f"

for (i in metrics) {
  x <- 1 # for naming the plots
  for (j in traits) {
    f <- formula(paste(i,"~", j))
    
    # fit lm for urban and get R2 and p
    lm.urban<- lm (f, data = BB22.ID[BB22.ID$landscape == "urban",])
    sum.urban <- summary(lm.urban)
    lab.urban <- sprintf(fmt, "urban", sum.urban$adj.r.squared, coef(sum.urban)[2, 4])
    
    # fit lm for rural and get R2 and p
    lm.rural <- lm (f, data = BB22.ID[BB22.ID$landscape == "rural",])
    sum.rural <- summary(lm.rural)
    lab.rural <- sprintf(fmt, "rural", sum.rural$adj.r.squared, coef(sum.rural)[2, 4])
    
    assign(paste("a", x, sep=""), # assign the ggplot to plot name
           # define the ggplot
           ggplot(BB22.ID, aes_string(j, i, colour = "landscape")) + 
             geom_point() + 
             theme_classic(base_size = 20) + 
             theme(aspect.ratio=1) + 
             labs(caption = paste(lab.urban, "\n", lab.rural)) +
             geom_smooth(method="lm", se = FALSE) +
             scale_color_manual(values=palette.landscape, labels=c("rural", "urban"))
    )
    x <- x+1
  } # end loop j
  setwd(output)
  plot4 <- ggarrange(a1,a2,a3,a4,a5,a6,a7,a8,a9, # arrange to plots nicely and export them 
                     ncol = 3, nrow = 3, 
                     labels = c(LETTERS[1:9]),   
                     common.legend = TRUE)
  annotate_figure(plot4, top = text_grob(paste("B.lapidarius: comparison of ", i, " and traits across landscapes", sep = ""),
                                         face = "bold", size = 22))
  ggsave(paste("./landscape comparison/functional diversity/lapi_ID/FD_lapi_corr_", i, "_BBtraits_landscapes.png", sep = ""), width = 15, height = 15)
  setwd(input)
} # end loop i


# PHYLOGENETIC SPECIES VARIABILITY ----
## preparation ----
# clear work environment
rm(list=ls()) 

#load libraries
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(glmmTMB)

# set working directory to main repository
input <- "YOUR PATH/Bumblebee_Simonetta_final/01_DATA" #change „YOUR PATH“ to where the repo is
output <- "YOUR PATH/Bumblebee_Simonetta_final/03_OUTPUT" #change „YOUR PATH“ to where the repo is
setwd(input)

# load data
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


## B.pascuorum ----
#only B.pascuroum 
BB22.bb.traits.sp <- BB22.bb.traits[BB22.bb.traits$bbspecies == "B.pascuorum",]

# import data with spatial information on sites
BB22.sites.meta <- read_csv("BB22_sites_2016.csv")

# import data with phylogenetic diversity of plants per site
BB22.fun.ID <- read_csv("./PD/PD_B.pascuorum_ID.short.csv")%>% 
  dplyr::select(-1)%>%
  rename_with(.cols = 1, ~"ID")%>%
  rename_with(.cols = 2, ~"sp_richn")%>%
  rename_with(.cols = 3, ~"pvar")%>%
  rename_with(.cols = 4, ~"pric")%>%
  rename_with(.cols = 5, ~"peve")%>%
  rename_with(.cols = 6, ~"pclu")

# filter functional metrics of plant traits to compare with traits per individual
BB22.fun.ID <- subset(BB22.fun.ID, ID %in% BB22.bb.traits.sp$ID)

# add site coordinates to the trait data frame (in LV95)
BB22.ID <- merge(BB22.bb.traits.sp, BB22.fun.ID, by  = "ID", all.x=TRUE)%>%
  mutate(site = paste(location, landscape, replicate, sep=""))
BB22.ID <- merge(BB22.ID, BB22.sites.meta[, c(1,2,3)], by  = "site", all.x=TRUE) 

### look at data ----------------------------------------------------------------------------------------

# Boxplots for all the variables we want to look at with Wilcoxon test
metrics <- colnames(BB22.ID[, c(18:21)]) # plant PD to look at; prepare for loop

library(rstatix)
plot_list  <- list()
w.test.list <- list()
mean.rural <- list()
mean.urban <- list()
palette.landscape <- c("#E69F00", "#56B4E9") #create color palette for landscape

for (j in metrics) {
  gg.data <- data.frame(landscape=BB22.ID$landscape,value=BB22.ID[[j]])
  w.test <- wilcox_test(gg.data,value~landscape)
  w.test.list[[j]] <- w.test
  mean.rural[[j]] <- mean(gg.data$value[gg.data$landscape == "R"], na.rm=TRUE)
  mean.urban[[j]] <- mean(gg.data$value[gg.data$landscape == "U"], na.rm=TRUE)
  
  p <- ggplot(gg.data, aes(x=landscape, y = value, fill=landscape)) + 
    geom_boxplot(notch = T) + 
    xlab("") +
    scale_x_discrete(labels=c('rural', 'urban'))+
    theme_classic(base_size = 20) +     
    theme(aspect.ratio=1) + 
    guides(alpha = "none") +
    scale_fill_manual(values=palette.landscape, guide = "none") + 
    labs(subtitle = paste("p = ", w.test$p, sep=""))
  if (j == "pric") {
    p <- p + ylim(0, 1) + ylab("phylogenetic species richness")
  } else if (j == "pclu") {
    p <- p + ylim(0, 1) + ylab("phylogenetic species clustering")
  } else if (j == "peve") {
    p <- p + ylim(0, 1) + ylab("phylogenetic species evenness")
  } else {
    p <- p + ylim(0, 1) + ylab("phylogenetic species variability")
  }
  plot_list[[j]] <- p
}

#### FIGURE 4b ----
# arrange them into one file to export
setwd(output)
plot <- ggarrange(plot_list[[1]],plot_list[[2]],
                  plot_list[[3]],plot_list[[4]],
                  ncol = 1, nrow = 4,
                  labels = c("A", "B", "C", "D"))
annotate_figure(plot, top = text_grob("B.pascuorum: phylogenetic diversity across landscapes", 
                                      face = "bold", size = 22))
ggsave("./landscape comparison/phylogenetic diversity/pasc_ID/PD_B.pascuorum_ID.png", width = 6, height = 24)
setwd(input)

# export statistics of W tests
w.stats.pasc <- rbind(c(mean.rural$pvar, mean.urban$pvar, w.test.list$pvar),
                      c(mean.rural$pric, mean.urban$pric, w.test.list$pric),
                      c(mean.rural$pclu, mean.urban$pclu, w.test.list$pclu),
                      c(mean.rural$peve, mean.urban$peve, w.test.list$peve))
setwd(output)
write.csv(w.stats.pasc,"./landscape comparison/phylogenetic diversity/pasc_ID/PD_W_pasc_BBtraits_species.csv")
setwd(input)

# plot the relationship of plants traits of one site and bumblebee traits of one site
traits <- colnames(BB22.ID[, 7:15]) # bumblebee traits to look at; prepare for loop

library(nlme)

# perform loop to output plots per relationship summarized per PD
# preparation for caption string
fmt <- "%s: adj.R^2 = %.3f, p = %.3f"

for (i in metrics) {
  x <- 1 # for naming the plots
  for (j in traits) {
    f <- formula(paste(i,"~", j))
  
    # fit lm for urban and get R2 and p
    lm.urban <- lm (f, data = BB22.ID[BB22.ID$landscape == "U",], na.action=na.exclude)
    sum.urban <- summary(lm.urban)
    lab.urban <- sprintf(fmt, "urban", sum.urban$adj.r.squared, coef(sum.urban)[2, 4])
    
    # fit lm for rural and get R2 and p
    lm.rural <- lm (f, data = BB22.ID[BB22.ID$landscape == "R",], na.action=na.exclude)
    sum.rural <- summary(lm.rural)
    lab.rural <- sprintf(fmt, "rural", sum.rural$adj.r.squared, coef(sum.rural)[2, 4])
    
    assign(paste("a", x, sep=""), # assign the ggplot to plot name
           # define the ggplot
           ggplot(BB22.ID, aes_string(j, i, colour = "landscape")) + 
             geom_point() + 
             theme_classic(base_size = 20) + 
             theme(aspect.ratio=1) + 
             labs(caption = paste(lab.urban, "\n", lab.rural)) +
             geom_smooth(method="lm", se = FALSE) +
             scale_color_manual(values=palette.landscape, labels=c("rural", "urban"))
    )
    x <- x+1
  } # end loop j
  setwd(output)
  plot4 <- ggarrange(a1,a2,a3,a4,a5,a6,a7,a8,a9, # arrange to plots nicely and export them 
                     ncol = 3, nrow = 3, 
                     labels = c(LETTERS[1:9]),   
                     common.legend = TRUE)
  annotate_figure(plot4, top = text_grob(paste("B.pascuorum: comparison of ", i, " and traits across landscapes", sep = ""),
                                         face = "bold", size = 22))
  ggsave(paste("./landscape comparison/phylogenetic diversity/pasc_ID/PD_pasc_corr_", i, "_BBtraits_landscapes.png", sep = ""), width = 15, height = 15)
  setwd(input)
} # end loop i

## B.LAPIDARIUS ----
#only B.lapidarius 
BB22.bb.traits.sp <- BB22.bb.traits[BB22.bb.traits$bbspecies == "B.lapidarius",]

# import data with spatial information on sites
BB22.sites.meta <- read_csv("BB22_sites_2016.csv")

# import data with phylogenetic diversity of plants per site
BB22.fun.ID <- read_csv("./PD/PD_B.lapidarius_ID.short.csv")%>% 
  dplyr::select(-1)%>%
  rename_with(.cols = 1, ~"ID")%>%
  rename_with(.cols = 2, ~"sp_richn")%>%
  rename_with(.cols = 3, ~"pvar")%>%
  rename_with(.cols = 4, ~"pric")%>%
  rename_with(.cols = 5, ~"peve")%>%
  rename_with(.cols = 6, ~"pclu")

# filter functional metrics of plant traits to compare with traits per individual
BB22.fun.ID <- subset(BB22.fun.ID, ID %in% BB22.bb.traits.sp$ID)

# add site coordinates to the trait data frame (in LV95)
BB22.ID <- merge(BB22.bb.traits.sp, BB22.fun.ID, by  = "ID", all.x=TRUE)%>%
  mutate(site = paste(location, landscape, replicate, sep=""))
BB22.ID <- merge(BB22.ID, BB22.sites.meta[, c(1,2,3)], by  = "site", all.x=TRUE) 

### look at data ----------------------------------------------------------------------------------------

# Boxplots for all the variables we want to look at with Wilcoxon test
metrics <- colnames(BB22.ID[, c(18:21)]) # plant PD to look at; prepare for loop

library(rstatix)
plot_list  <- list()
w.test.list <- list()
mean.rural <- list()
mean.urban <- list()
palette.landscape <- c("#E69F00", "#56B4E9") #create color palette for landscape

for (j in metrics) {
  gg.data <- data.frame(landscape=BB22.ID$landscape,value=BB22.ID[[j]])
  w.test <- wilcox_test(gg.data,value~landscape)
  w.test.list[[j]] <- w.test
  mean.rural[[j]] <- mean(gg.data$value[gg.data$landscape == "R"], na.rm=TRUE)
  mean.urban[[j]] <- mean(gg.data$value[gg.data$landscape == "U"], na.rm=TRUE)
  
  p <- ggplot(gg.data, aes(x=landscape, y = value, fill=landscape)) + 
    geom_boxplot(notch = T) + 
    xlab("") +
    scale_x_discrete(labels=c('rural', 'urban'))+
    theme_classic(base_size = 20) +     
    theme(aspect.ratio=1) + 
    guides(alpha = "none") +
    scale_fill_manual(values=palette.landscape, guide = "none") + 
    labs(subtitle = paste("p = ", w.test$p, sep=""))
  if (j == "pric") {
    p <- p + ylim(0, 1) + ylab("phylogenetic species richness")
  } else if (j == "pclu") {
    p <- p + ylim(0, 1) + ylab("phylogenetic species clustering")
  } else if (j == "peve") {
    p <- p + ylim(0, 1) + ylab("phylogenetic species evenness")
  } else {
    p <- p + ylim(0, 1) + ylab("phylogenetic species variability")
  }
  plot_list[[j]] <- p
}

#### FIGURE 4a ----
# arrange them into one file to export
setwd(output)
plot <- ggarrange(plot_list[[1]],plot_list[[2]],
                  plot_list[[3]],plot_list[[4]],
                  ncol = 1, nrow = 4,
                  labels = c("A", "B", "C", "D"))
annotate_figure(plot, top = text_grob("B.lapidarius: phylogenetic diversity across landscapes", 
                                      face = "bold", size = 22))
ggsave("./landscape comparison/phylogenetic diversity/lapi_ID/PD_B.lapidarius_ID.png", width = 6, height = 24)
setwd(input)

# export statistics of W tests
w.stats.lapi <- rbind(c(mean.rural$pvar, mean.urban$pvar, w.test.list$pvar),
                      c(mean.rural$pric, mean.urban$pric, w.test.list$pric),
                      c(mean.rural$pclu, mean.urban$pclu, w.test.list$pclu),
                      c(mean.rural$peve, mean.urban$peve, w.test.list$peve))
setwd(output)
write.csv(w.stats.lapi,"./landscape comparison/phylogenetic diversity/lapi_ID/PD_W_lapi_BBtraits_species.csv")
setwd(input)

# plot the relationship of plants traits of one site and bumblebee traits of one site
traits <- colnames(BB22.ID[, 7:15]) # bumblebee traits to look at; prepare for loop

library(nlme)

# perform loop to output plots per relationship summarized per FD
# preparation for caption string
fmt <- "%s: adj.R^2 = %.3f, p = %.3f"

for (i in metrics) {
  x <- 1 # for naming the plots
  for (j in traits) {
    f <- formula(paste(i,"~", j))
   
    # fit lm for urban and get R2 and p
    lm.urban <- lm (f, data = BB22.ID[BB22.ID$landscape == "U",], na.action=na.exclude)
    sum.urban <- summary(lm.urban)
    lab.urban <- sprintf(fmt, "urban", sum.urban$adj.r.squared, coef(sum.urban)[2, 4])
    
    # fit lm for rural and get R2 and p
    lm.rural <- lm (f, data = BB22.ID[BB22.ID$landscape == "R",], na.action=na.exclude)
    sum.rural <- summary(lm.rural)
    lab.rural <- sprintf(fmt, "rural", sum.rural$adj.r.squared, coef(sum.rural)[2, 4])
    
    assign(paste("a", x, sep=""), # assign the ggplot to plot name
           # define the ggplot
           ggplot(BB22.ID, aes_string(j, i, colour = "landscape")) + 
             geom_point() + 
             theme_classic(base_size = 20) + 
             theme(aspect.ratio=1) + 
             labs(caption = paste(lab.urban, "\n", lab.rural)) +
             geom_smooth(method="lm", se = FALSE) +
             scale_color_manual(values=palette.landscape, labels=c("rural", "urban"))
    )
    x <- x+1
  } # end loop j
  setwd(output)
  plot4 <- ggarrange(a1,a2,a3,a4,a5,a6,a7,a8,a9, # arrange to plots nicely and export them 
                     ncol = 3, nrow = 3, 
                     labels = c(LETTERS[1:9]),   
                     common.legend = TRUE)
  annotate_figure(plot4, top = text_grob(paste("B.lapidarius: comparison of ", i, " and traits across landscapes", sep = ""),
                                         face = "bold", size = 22))
  ggsave(paste("./landscape comparison/phylogenetic diversity/lapi_ID/PD_lapi_corr_", i, "_BBtraits_landscapes.png", sep = ""), 
         width = 15, height = 15)
  setwd(input)
} # end loop i


# DIET CONSISTENCY ----

## preparation ----
# clear environment
rm(list=ls())

# load required library
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(rlist)

# set working directory to main repository
input <- "YOUR PATH/Bumblebee_Simonetta_final/01_DATA" #change „YOUR PATH“ to where the repo is
output <- "YOUR PATH/Bumblebee_Simonetta_final/03_OUTPUT" #change „YOUR PATH“ to where the repo is
setwd(input)

# load data
BB22.full <- read_csv("BB22_full.csv") %>%
  mutate(site = as_factor(site),
         species = as_factor(species),
         region = substr(site, 1, 3)) 

## B.pascuorum ----
# only use B.pascuorum
BB22.full.species <- BB22.full[BB22.full$bbspecies == "B.pascuorum",]

### site level  ----
# prepare list element to store correlation matrices for species, genus and family
pasc.site <- list()

#### species ----
# convert into binary data frame for comaring and later calculate correlations between sites
BB22.full.table <- t(as.data.frame.matrix(+(table(BB22.full.species[, c("site", "plant.species")]) > 0)))
BB22.full.table <- as.data.frame(BB22.full.table)
BB22.full.table$species <- rownames(BB22.full.table)

# create species lists per site
sitenames <- c("ZHUA", "ZHUB", "ZHUC", "ZHRD", "ZHRE", "ZHRF", 
               "BEUA", "BEUB", "BEUC", "BERD", "BSUA", "BSUB", 
               "BSUC", "BSRD", "BSRE", "BSRF")

site.site.bb <- list()

for (i in sitenames) {
  site.site.bb[[i]]  <- unique(BB22.full.species$plant.species[BB22.full.species$site == i])
} # end loop i


# import species lists per site from pollen and GBIF/InfoFlora (file: data preparation)
# load("site_list_bb.RData")
# site.species.bb <- site.list
# 
# load("sp_list_gbif_infoflora.RData")
# site.species.occ <- site.list.1500

# import data from GBIF and InfoFlora and convert into species lists per site
# InfoFlora
#prepare list to store species lists
site.list.InfoFlora.species <- list()
for (i in sitenames) {
  # import data
  site.data.InfoFlora  <- read_csv(paste("./GBIF and InfoFlora/InfoFlora_", i, "_1500.csv", sep = ""))
  # summarize only needed variables (site and species)
  site.data.InfoFlora <- site.data.InfoFlora %>%
    summarise(Species = Taxon,
              site = rep(i, nrow(site.data.InfoFlora)),
              region = rep(substr(i, 1, 3), nrow(site.data.InfoFlora)))
  j <- 0
  site.data.InfoFlora$species <- c()
  for (j in 1:nrow(site.data.InfoFlora)){
    # species names are different than in GBIF -> need to be adapted
    site.data.InfoFlora$species[j] <- paste(unlist(strsplit(site.data.InfoFlora$Species[j], split=' ', fixed=TRUE))[1],
                                            unlist(strsplit(site.data.InfoFlora$Species[j], split=' ', fixed=TRUE))[2])
  } # end loop j

  # store species per site in a list
  site.list.InfoFlora.species[[i]] <- unique(site.data.InfoFlora$species)

} # end loop i

# GBIF
#prepare list to store species lists
site.list.gbif.species <- list()
for (i in sitenames) {
  # import data and filter for needed variables
  site.data.gbif  <- read_csv(paste("./GBIF and InfoFlora/GBIF_", i, "_1500.csv", sep = "")) %>%
    mutate(species = as_factor(species)) %>%
    filter(class == "Magnoliopsida" | class == "Liliopsida")
  site.list.gbif.species[[i]] <- unique(site.data.gbif$species)
} # end loop i

# combine species list GBIF and InfoFlora (site level) in new list
site.species.occ <- list()
for (i in sitenames) {
  x <- as.factor(site.list.InfoFlora.species[[i]])
  y <- as.factor(site.list.gbif.species[[i]])
  site.species.occ[[i]] <- unique(c(x,y))%>%
    droplevels()
} # end loop i

# create list with intersection of species lists for all sites
intersection.sites.species <- list()
table.sites.species <- list()
correlation.site.species <- data.frame(species = NA, site1 = NA, site2 = NA, Correlation = NA)

# compute and store correlations between each site
for (i in sitenames) {
  for (j in sitenames) {
    # find common species of two sites
    intersection.sites.species[[paste(i,"/",j, sep="")]] <- intersect(site.species.occ[[i]], site.species.occ[[j]])
    
    # filter binary dataframe for the two sites
    table.sites.species[[paste(i,"/",j, sep="")]] <- dplyr::select(BB22.full.table, i, j)
    
    # compare bumblebee pollen with pollen lists of landscape
    comparison.table <- BB22.full.table[BB22.full.table$species %in% intersection.sites.species[[paste(i,"/",j, sep="")]], c(i,j)]
    
    # compute correlation and store them in a dataframe
    correlation <- c("B.pascuorum", i, j, round(cor(comparison.table[,1], comparison.table[,2], method=c("pearson")),3))
    correlation.site.species <- rbind(correlation.site.species, correlation)
    
  } # end loop j
} # end loop i

correlation.site.species$site1 <- as.factor(correlation.site.species$site1)
matrix.site.species <- split(correlation.site.species, f = correlation.site.species$site1)
temp <- list.cbind(matrix.site.species)[, seq(4, 67, 4)]
colnames(temp) <- levels(correlation.site.species$site1)
rownames(temp) <- c("ZHUA", "ZHUB", "ZHUC", "ZHRD", "ZHRE", "ZHRF", "BEUA", "BEUB", "BEUC", "BERD", "BSUA", "BSUB", "BSUC", "BSRD", "BSRE", "BSRF")
matrix.site.species <- as.matrix(temp[,match(sitenames, colnames(temp))])
class(matrix.site.species) <- "numeric"

# store in list
pasc.site[["species"]] <- abs(matrix.site.species)


#### genus ----
# convert into binary data frame
BB22.full.table <- t(as.data.frame.matrix(+(table(BB22.full.species[, c("site", "genus")]) > 0)))
BB22.full.table <- as.data.frame(BB22.full.table)
BB22.full.table$genus <- rownames(BB22.full.table)

# create genus lists per site
sitenames <- c("ZHUA", "ZHUB", "ZHUC", "ZHRD", "ZHRE", "ZHRF", "BEUA", "BEUB", 
               "BEUC", "BERD", "BSUA", "BSUB", "BSUC", "BSRD", "BSRE", "BSRF")

site.genus.bb <- list()

for (i in sitenames) {
  site.genus.bb[[i]]  <- unique(BB22.full.species$genus[BB22.full.species$site == i])
} # end loop i

# import data from GBIF and InfoFlora and convert into species lists per site
# InfoFlora
library(stringr)

#prepare list to store species lists
site.list.InfoFlora.genus <- list()
for (i in sitenames) {
  # import data
  site.data.InfoFlora  <- read_csv(paste("./GBIF and InfoFlora/InfoFlora_", i, "_1500.csv", sep = ""))
  # summarize only needed variables (site and species)
  site.data.InfoFlora <- site.data.InfoFlora %>%
    summarise(genus = word(Taxon, 1),
              site = rep(i, nrow(site.data.InfoFlora)),
              region = rep(substr(i, 1, 3), nrow(site.data.InfoFlora)))

  site.list.InfoFlora.genus[[i]] <- unique(site.data.InfoFlora$genus)
} # end loop i

# GBIF
#prepare list to store species lists
site.list.gbif.genus <- list()
for (i in sitenames) {
  # import data and filter for needed variables
  site.data.gbif  <- read_csv(paste("./GBIF and InfoFlora/GBIF_", i, "_1500.csv", sep = "")) %>%
    mutate(genus = as_factor(genus)) %>%
    filter(class == "Magnoliopsida" | class == "Liliopsida")
  site.list.gbif.genus[[i]] <- unique(site.data.gbif$genus)
} # end loop i

# combine species list GBIF and InfoFlora (site level) in new list
site.genus.occ <- list()
for (i in sitenames) {
  x <- as.factor(site.list.InfoFlora.genus[[i]])
  y <- as.factor(site.list.gbif.genus[[i]])
  site.genus.occ[[i]] <- unique(c(x,y)) %>%
    droplevels()
} # end loop i

# create list with intersection of genus lists for all sites
intersection.sites.genus <- list()
table.sites.genus <- list()
correlation.site.genus <- data.frame(genus = NA, site1 = NA, site2 = NA, Correlation = NA)

# compute and store correlations between each site
for (i in sitenames) {
  for (j in sitenames) {
    # find common genus of two sites
    intersection.sites.genus[[paste(i,"/",j, sep="")]] <- intersect(site.genus.occ[[i]], site.genus.occ[[j]])
    
    # filter binary dataframe for the two sites
    table.sites.genus[[paste(i,"/",j, sep="")]] <- dplyr::select(BB22.full.table, i, j)
    
    # compare bumblebee pollen with pollen lists of landscape
    comparison.table <- BB22.full.table[BB22.full.table$genus %in% intersection.sites.genus[[paste(i,"/",j, sep="")]], c(i,j)]
    
    # compute correlation and store them in a dataframe
    correlation <- c("B.pascuorum", i, j, round(cor(comparison.table[,1], comparison.table[,2], method=c("pearson")),3))
    correlation.site.genus <- rbind(correlation.site.genus, correlation)
    
  } # end loop j
} # end loop i

correlation.site.genus$site1 <- as.factor(correlation.site.genus$site1)
matrix.site.genus <- split(correlation.site.genus, f = correlation.site.genus$site1)
temp <- list.cbind(matrix.site.genus)[, seq(4, 67, 4)]
colnames(temp) <- levels(correlation.site.genus$site1)
rownames(temp) <- c("ZHUA", "ZHUB", "ZHUC", "ZHRD", "ZHRE", "ZHRF", "BEUA", "BEUB", "BEUC", "BERD", "BSUA", "BSUB", "BSUC", "BSRD", "BSRE", "BSRF")
matrix.site.genus <- as.matrix(temp[,match(sitenames, colnames(temp))])
class(matrix.site.genus) <- "numeric"

# store in list
pasc.site[["genus"]] <- abs(matrix.site.genus)


#### family ----

# convert into binary data frame
BB22.full.table <- t(as.data.frame.matrix(+(table(BB22.full.species[, c("site", "family")]) > 0)))
BB22.full.table <- as.data.frame(BB22.full.table)
BB22.full.table$family <- rownames(BB22.full.table)

# create genus lists per site
sitenames <- c("ZHUA", "ZHUB", "ZHUC", "ZHRD", "ZHRE", "ZHRF", "BEUA", "BEUB", "BEUC", "BERD", "BSUA", "BSUB", "BSUC", "BSRD", "BSRE", "BSRF")
site.family.bb <- list()

for (i in sitenames) {
  site.family.bb[[i]]  <- unique(BB22.full.species$family[BB22.full.species$site == i])
} # end loop i

# import data from GBIF and InfoFlora and convert into species lists per site
# InfoFlora
#prepare list to store species lists
site.list.InfoFlora.family <- list()
for (i in sitenames) {
  # import data
  site.data.InfoFlora  <- read_csv(paste("./GBIF and InfoFlora/InfoFlora_", i, "_1500.csv", sep = ""))
  # summarize only needed variables (site and species)
  site.data.InfoFlora <- site.data.InfoFlora %>%
    summarise(family = Familie,
              site = rep(i, nrow(site.data.InfoFlora)),
              region = rep(substr(i, 1, 3), nrow(site.data.InfoFlora)))
  
  site.list.InfoFlora.family[[i]] <- unique(site.data.InfoFlora$family)
} # end loop i

# GBIF
#prepare list to store species lists
site.list.gbif.family <- list()
for (i in sitenames) {
  # import data and filter for needed variables
  site.data.gbif  <- read_csv(paste("./GBIF and InfoFlora/GBIF_", i, "_1500.csv", sep = "")) %>%
    mutate(family = as_factor(family)) %>%
    filter(class == "Magnoliopsida" | class == "Liliopsida")
  site.list.gbif.family[[i]] <- unique(site.data.gbif$family)
} # end loop i

# combine species list GBIF and InfoFlora (site level) in new list
site.family.occ <- list()
for (i in sitenames) {
  x <- as.factor(site.list.InfoFlora.family[[i]])
  y <- as.factor(site.list.gbif.family[[i]])
  site.family.occ[[i]] <- unique(c(x,y)) %>% 
    droplevels()
} # end loop i

# create list with intersection of family lists for all sites
intersection.sites.family <- list()
table.sites.family <- list()
correlation.site.family <- data.frame(family = NA, site1 = NA, site2 = NA, Correlation = NA)

# compute and store correlations between each site
for (i in sitenames) {
  for (j in sitenames) {
    # find common family of two sites
    intersection.sites.family[[paste(i,"/",j, sep="")]] <- intersect(site.family.occ[[i]], site.family.occ[[j]])
    
    # filter binary dataframe for the two sites
    table.sites.family[[paste(i,"/",j, sep="")]] <- dplyr::select(BB22.full.table, i, j)
    
    # compare bumblebee pollen with pollen lists of landscape
    comparison.table <- BB22.full.table[BB22.full.table$family %in% intersection.sites.family[[paste(i,"/",j, sep="")]], c(i,j)]
    
    # compute correlation and store them in a dataframe
    correlation <- c("B.pascuorum", i, j, round(cor(comparison.table[,1], comparison.table[,2], method=c("pearson")),3))
    correlation.site.family <- rbind(correlation.site.family, correlation)
    
  } # end loop j
} # end loop i

correlation.site.family$site1 <- as.factor(correlation.site.family$site1)
matrix.site.family <- split(correlation.site.family, f = correlation.site.family$site1)
temp <- list.cbind(matrix.site.family)[, seq(4, 67, 4)]
colnames(temp) <- levels(correlation.site.family$site1)
rownames(temp) <- c("ZHUA", "ZHUB", "ZHUC", "ZHRD", "ZHRE", "ZHRF", "BEUA", "BEUB", "BEUC", "BERD", "BSUA", "BSUB", "BSUC", "BSRD", "BSRE", "BSRF")
matrix.site.family <- as.matrix(temp[,match(sitenames, colnames(temp))])
class(matrix.site.family) <- "numeric"

# store in list
pasc.site[["family"]] <- abs(matrix.site.family)

#### plotting ----
library(ggcorrplot)
axis.order.site <- c("BERD", 
                     "BSRD", "BSRE", "BSRF",
                     "ZHRD", "ZHRE", "ZHRF",
                     "BEUA", "BEUB", "BEUC",
                     "BSUA", "BSUB", "BSUC", 
                     "ZHUA", "ZHUB", "ZHUC") # to consistently order the axis

a1 <- ggcorrplot(pasc.site[["species"]][axis.order.site, axis.order.site], 
                 hc.order = F, type = "lower", 
                 outline.col = "white")+ 
  labs(title ='species', fill = "", x = "", y = "") +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_gradient2(limit = c(0,1), low = "white", high =  "#07575B")

a2 <- ggcorrplot(pasc.site[["genus"]][axis.order.site, axis.order.site], 
                 hc.order = F, type = "lower", 
                 outline.col = "white")+ 
  labs(title ='genus', fill = "", x = "", y = "") +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_gradient2(limit = c(0,1), low = "white", high =  "#07575B")

a3 <- ggcorrplot(pasc.site[["family"]][axis.order.site, axis.order.site], 
                 hc.order = F, type = "lower", 
                 outline.col = "white")+ 
  labs(title ='family', fill = "", x = "", y = "") +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_gradient2(limit = c(0,1), low = "white", high =  "#07575B")

##### FIGURE S10 ----
# arrange them into one file to export
setwd(output)
plot1 <- ggarrange(a1, a2, a3, ncol = 1, nrow = 3,
                   labels = c("A", "B", "C"), common.legend = TRUE, legend = "right")
annotate_figure(plot1, top = text_grob("B.pascuroum", face = "bold", size = 22))
ggsave("./landscape comparison/diet consistency/corr_pasc_site.png", width = 10, height = 30)
setwd(input)

### region level ----
# prepare list element to store correlation matrices for species, genus and family
pasc.region <- list()

#### species ----

# convert into binary data frame for comaring and later calculate correlations between regions
BB22.full.table <- t(as.data.frame.matrix(+(table(BB22.full.species[, c("region", "plant.species")]) > 0)))
BB22.full.table <- as.data.frame(BB22.full.table)
BB22.full.table$species <- rownames(BB22.full.table)

# create species lists per region
regionnames <- c("ZHU", "ZHR", "BSU", "BSR", "BEU", "BER")
region.species.bb <- list()

for (i in regionnames) {
  region.species.bb[[i]]  <- unique(BB22.full.species$plant.species[BB22.full.species$region == i])
} # end loop i

# import data from GBIF and InfoFlora and convert into species lists per site
# InfoFlora
# combine species list per site into regional level
region.list.InfoFlora <- list()
region.list.InfoFlora[["ZHU"]] <- unique(c(site.list.InfoFlora.species$ZHUA, site.list.InfoFlora.species$ZHUB, site.list.InfoFlora.species$ZHUC))
region.list.InfoFlora[["ZHR"]] <- unique(c(site.list.InfoFlora.species$ZHRD, site.list.InfoFlora.species$ZHRE, site.list.InfoFlora.species$ZHRF))
region.list.InfoFlora[["BSU"]] <- unique(c(site.list.InfoFlora.species$BSUA, site.list.InfoFlora.species$BSUB, site.list.InfoFlora.species$BSUC))
region.list.InfoFlora[["BSR"]] <- unique(c(site.list.InfoFlora.species$BSRD, site.list.InfoFlora.species$BSRE, site.list.InfoFlora.species$BSRF))
region.list.InfoFlora[["BEU"]] <- unique(c(site.list.InfoFlora.species$BEUA, site.list.InfoFlora.species$BEUB, site.list.InfoFlora.species$BEUC))
region.list.InfoFlora[["BER"]] <- unique(c(site.list.InfoFlora.species$BERD))

# GBIF
# combine species list per site into regional level
region.list.gbif <- list()
region.list.gbif[["ZHU"]] <- unique(c(site.list.gbif.species$ZHUA, site.list.gbif.species$ZHUB, site.list.gbif.species$ZHUC))
region.list.gbif[["ZHR"]] <- unique(c(site.list.gbif.species$ZHRD, site.list.gbif.species$ZHRE, site.list.gbif.species$ZHRF))
region.list.gbif[["BSU"]] <- unique(c(site.list.gbif.species$BSUA, site.list.gbif.species$BSUB, site.list.gbif.species$BSUC))
region.list.gbif[["BSR"]] <- unique(c(site.list.gbif.species$BSRD, site.list.gbif.species$BSRE, site.list.gbif.species$BSRF))
region.list.gbif[["BEU"]] <- unique(c(site.list.gbif.species$BEUA, site.list.gbif.species$BEUB, site.list.gbif.species$BEUC))
region.list.gbif[["BER"]] <- unique(c(site.list.gbif.species$BERD))


# combine species list GBIF and InfoFlora (region level) in new list
region.species.occ <- list()
for (i in regionnames) {
  x <- as.factor(region.list.InfoFlora[[i]])
  y <- as.factor(region.list.gbif[[i]])
  region.species.occ[[i]] <- unique(c(x,y))%>% 
    droplevels()
} # end loop i

# create list with intersection of species lists for all regions
intersection.regions.species <- list()
table.regions.species <- list()
correlation.region.species <- data.frame(species = NA, region1 = NA, region2 = NA, Correlation = NA)

# compute and store correlations between each region
for (i in regionnames) {
  for (j in regionnames) {
    # find common species of two regions
    intersection.regions.species[[paste(i,"/",j, sep="")]] <- intersect(region.species.occ[[i]], region.species.occ[[j]])
    
    # filter binary dataframe for the two regions
    table.regions.species[[paste(i,"/",j, sep="")]] <- dplyr::select(BB22.full.table, i, j)
    
    # compare bumblebee pollen with pollen lists of landscape
    comparison.table <- BB22.full.table[BB22.full.table$species %in% intersection.regions.species[[paste(i,"/",j, sep="")]], c(i,j)]
    
    # compute correlation and store them in a dataframe
    correlation <- c("B.pascuorum", i, j, round(cor(comparison.table[,1], comparison.table[,2], method=c("pearson")),3))
    correlation.region.species <- rbind(correlation.region.species, correlation)
    
  } # end loop j
} # end loop i

correlation.region.species$region1 <- as.factor(correlation.region.species$region1)
matrix.region.species <- split(correlation.region.species, f = correlation.region.species$region1)
temp <- list.cbind(matrix.region.species)[, seq(4, 24, 4)]
colnames(temp) <- levels(correlation.region.species$region1)
rownames(temp) <- regionnames
matrix.region.species <- as.matrix(temp[,match(regionnames, colnames(temp))])
class(matrix.region.species) <- "numeric"

# store in list
pasc.region[["species"]] <- abs(matrix.region.species)


#### genus ----
# convert into binary data frame
BB22.full.table <- t(as.data.frame.matrix(+(table(BB22.full.species[, c("region", "genus")]) > 0)))
BB22.full.table <- as.data.frame(BB22.full.table)
BB22.full.table$genus <- rownames(BB22.full.table)

# create genus lists per region
region.genus.bb <- list()

for (i in regionnames) {
  region.genus.bb[[i]]  <- unique(BB22.full.species$genus[BB22.full.species$region == i])
} # end loop i

# import data from GBIF and InfoFlora and convert into species lists per site
# InfoFlora
# combine species list per site into regional level
region.list.InfoFlora <- list()
region.list.InfoFlora[["ZHU"]] <- unique(c(site.list.InfoFlora.genus$ZHUA, site.list.InfoFlora.genus$ZHUB, site.list.InfoFlora.genus$ZHUC))
region.list.InfoFlora[["ZHR"]] <- unique(c(site.list.InfoFlora.genus$ZHRD, site.list.InfoFlora.genus$ZHRE, site.list.InfoFlora.genus$ZHRF))
region.list.InfoFlora[["BSU"]] <- unique(c(site.list.InfoFlora.genus$BSUA, site.list.InfoFlora.genus$BSUB, site.list.InfoFlora.genus$BSUC))
region.list.InfoFlora[["BSR"]] <- unique(c(site.list.InfoFlora.genus$BSRD, site.list.InfoFlora.genus$BSRE, site.list.InfoFlora.genus$BSRF))
region.list.InfoFlora[["BEU"]] <- unique(c(site.list.InfoFlora.genus$BEUA, site.list.InfoFlora.genus$BEUB, site.list.InfoFlora.genus$BEUC))
region.list.InfoFlora[["BER"]] <- unique(c(site.list.InfoFlora.genus$BERD))

# GBIF
# combine species list per site into regional level
region.list.gbif <- list()
region.list.gbif[["ZHU"]] <- unique(c(site.list.gbif.genus$ZHUA, site.list.gbif.genus$ZHUB, site.list.gbif.genus$ZHUC))
region.list.gbif[["ZHR"]] <- unique(c(site.list.gbif.genus$ZHRD, site.list.gbif.genus$ZHRE, site.list.gbif.genus$ZHRF))
region.list.gbif[["BSU"]] <- unique(c(site.list.gbif.genus$BSUA, site.list.gbif.genus$BSUB, site.list.gbif.genus$BSUC))
region.list.gbif[["BSR"]] <- unique(c(site.list.gbif.genus$BSRD, site.list.gbif.genus$BSRE, site.list.gbif.genus$BSRF))
region.list.gbif[["BEU"]] <- unique(c(site.list.gbif.genus$BEUA, site.list.gbif.genus$BEUB, site.list.gbif.genus$BEUC))
region.list.gbif[["BER"]] <- unique(c(site.list.gbif.genus$BERD))

# combine species list GBIF and InfoFlora (region level) in new list
region.genus.occ <- list()
for (i in regionnames) {
  x <- as.factor(region.list.InfoFlora[[i]])
  y <- as.factor(region.list.gbif[[i]])
  region.genus.occ[[i]] <- unique(c(x,y)) %>% 
    droplevels()
} # end loop i

# create list with intersection of genus lists for all regions
intersection.regions.genus <- list()
table.regions.genus <- list()
correlation.region.genus <- data.frame(genus = NA, region1 = NA, region2 = NA, Correlation = NA)

# compute and store correlations between each region
for (i in regionnames) {
  for (j in regionnames) {
    # find common genus of two regions
    intersection.regions.genus[[paste(i,"/",j, sep="")]] <- intersect(region.genus.occ[[i]], region.genus.occ[[j]])
    
    # filter binary dataframe for the two regions
    table.regions.genus[[paste(i,"/",j, sep="")]] <- dplyr::select(BB22.full.table, i, j)
    
    # compare bumblebee pollen with pollen lists of landscape
    comparison.table <- BB22.full.table[BB22.full.table$genus %in% intersection.regions.genus[[paste(i,"/",j, sep="")]], c(i,j)]
    
    # compute correlation and store them in a dataframe
    correlation <- c("B.pascuorum", i, j, round(cor(comparison.table[,1], comparison.table[,2], method=c("pearson")),3))
    correlation.region.genus <- rbind(correlation.region.genus, correlation)
    
  } # end loop j
} # end loop i

correlation.region.genus$region1 <- as.factor(correlation.region.genus$region1)
matrix.region.genus <- split(correlation.region.genus, f = correlation.region.genus$region1)
temp <- list.cbind(matrix.region.genus)[, seq(4, 24, 4)]
colnames(temp) <- levels(correlation.region.genus$region1)
rownames(temp) <- regionnames
matrix.region.genus <- as.matrix(temp[,match(regionnames, colnames(temp))])
class(matrix.region.genus) <- "numeric"

# store in list
pasc.region[["genus"]] <- abs(matrix.region.genus)



#### family ----
# convert into binary data frame
BB22.full.table <- t(as.data.frame.matrix(+(table(BB22.full.species[, c("region", "family")]) > 0)))
BB22.full.table <- as.data.frame(BB22.full.table)
BB22.full.table$family <- rownames(BB22.full.table)

# create genus lists per site
region.family.bb <- list()

for (i in regionnames) {
  region.family.bb[[i]]  <- unique(BB22.full.species$family[BB22.full.species$region == i])
} # end loop i

# InfoFlora
# combine species list per site into regional level
region.list.InfoFlora <- list()
region.list.InfoFlora[["ZHU"]] <- unique(c(site.list.InfoFlora.family$ZHUA, site.list.InfoFlora.family$ZHUB, site.list.InfoFlora.family$ZHUC))
region.list.InfoFlora[["ZHR"]] <- unique(c(site.list.InfoFlora.family$ZHRD, site.list.InfoFlora.family$ZHRE, site.list.InfoFlora.family$ZHRF))
region.list.InfoFlora[["BSU"]] <- unique(c(site.list.InfoFlora.family$BSUA, site.list.InfoFlora.family$BSUB, site.list.InfoFlora.family$BSUC))
region.list.InfoFlora[["BSR"]] <- unique(c(site.list.InfoFlora.family$BSRD, site.list.InfoFlora.family$BSRE, site.list.InfoFlora.family$BSRF))
region.list.InfoFlora[["BEU"]] <- unique(c(site.list.InfoFlora.family$BEUA, site.list.InfoFlora.family$BEUB, site.list.InfoFlora.family$BEUC))
region.list.InfoFlora[["BER"]] <- unique(c(site.list.InfoFlora.family$BERD))

# GBIF
# combine species list per site into regional level
region.list.gbif <- list()
region.list.gbif[["ZHU"]] <- unique(c(site.list.gbif.family$ZHUA, site.list.gbif.family$ZHUB, site.list.gbif.family$ZHUC))
region.list.gbif[["ZHR"]] <- unique(c(site.list.gbif.family$ZHRD, site.list.gbif.family$ZHRE, site.list.gbif.family$ZHRF))
region.list.gbif[["BSU"]] <- unique(c(site.list.gbif.family$BSUA, site.list.gbif.family$BSUB, site.list.gbif.family$BSUC))
region.list.gbif[["BSR"]] <- unique(c(site.list.gbif.family$BSRD, site.list.gbif.family$BSRE, site.list.gbif.family$BSRF))
region.list.gbif[["BEU"]] <- unique(c(site.list.gbif.family$BEUA, site.list.gbif.family$BEUB, site.list.gbif.family$BEUC))
region.list.gbif[["BER"]] <- unique(c(site.list.gbif.family$BERD))

# combine species list GBIF and InfoFlora (region level) in new list
region.family.occ <- list()
for (i in regionnames) {
  x <- as.factor(region.list.InfoFlora[[i]])
  y <- as.factor(region.list.gbif[[i]])
  region.family.occ[[i]] <- unique(c(x,y)) %>% 
    droplevels()
} # end loop i

# create list with intersection of family lists for all regions
intersection.regions.family <- list()
table.regions.family <- list()
correlation.region.family <- data.frame(family = NA, region1 = NA, region2 = NA, Correlation = NA)

# compute and store correlations between each region
for (i in regionnames) {
  for (j in regionnames) {
    # find common family of two regions
    intersection.regions.family[[paste(i,"/",j, sep="")]] <- intersect(region.family.occ[[i]], region.family.occ[[j]])
    
    # filter binary dataframe for the two regions
    table.regions.family[[paste(i,"/",j, sep="")]] <- dplyr::select(BB22.full.table, i, j)
    
    # compare bumblebee pollen with pollen lists of landscape
    comparison.table <- BB22.full.table[BB22.full.table$family %in% intersection.regions.family[[paste(i,"/",j, sep="")]], c(i,j)]
    
    # compute correlation and store them in a dataframe
    correlation <- c("B.pascuorum", i, j, round(cor(comparison.table[,1], comparison.table[,2], method=c("pearson")),3))
    correlation.region.family <- rbind(correlation.region.family, correlation)
    
  } # end loop j
} # end loop i

correlation.region.family$region1 <- as.factor(correlation.region.family$region1)
matrix.region.family <- split(correlation.region.family, f = correlation.region.family$region1)
temp <- list.cbind(matrix.region.family)[, seq(4, 24, 4)]
colnames(temp) <- levels(correlation.region.family$region1)
rownames(temp) <- regionnames
matrix.region.family <- as.matrix(temp[,match(regionnames, colnames(temp))])
class(matrix.region.family) <- "numeric"

# store in list
pasc.region[["family"]] <- abs(matrix.region.family)

#### plotting ----
library(ggcorrplot)
axis.order.reg <- c("BER", "BSR", "ZHR", "BEU", "BSU", "ZHU") # to consistenly order the axis

b1 <- ggcorrplot(pasc.region[["species"]][axis.order.reg, axis.order.reg], 
                 hc.order = F, type = "lower",
                 outline.col = "white")+ 
  labs(title ='species', fill = "", x = "", y = "") +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_gradient2(limit = c(0,1), low = "white", high =  "#07575B")


b2 <- ggcorrplot(pasc.region[["genus"]][axis.order.reg, axis.order.reg], 
                 hc.order = F, type = "lower", 
                 outline.col = "white")+ 
  labs(title ='genus', fill = "", x = "", y = "") +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_gradient2(limit = c(0,1), low = "white", high =  "#07575B")

b3 <- ggcorrplot(pasc.region[["family"]][axis.order.reg, axis.order.reg], 
                 hc.order = F, type = "lower", 
                 outline.col = "white")+ 
  labs(title ='family', fill = "", x = "", y = "") +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_gradient2(limit = c(0,1), low = "white", high =  "#07575B")

##### FIGURE 5 ----
# arrange them into one file to export
setwd(output)
plot2 <- ggarrange(b1,b2,b3, ncol = 1, nrow = 3,
                   labels = c("A", "B", "C"), common.legend = TRUE, legend = "right")
annotate_figure(plot2, top = text_grob("B.pascuroum", face = "bold", size = 22))
ggsave("./landscape comparison/diet consistency/corr_pasc_region.png", width = 10, height = 30)
setwd(input)

## B.lapidarius ----
# only use B.lapidarius
BB22.full.species <- BB22.full[BB22.full$bbspecies == "B.lapidarius",]

### site level  ----
# prepare list element to store correlation matrices for species, genus and family
lapi.site <- list()

#### species ----
# convert into binary data frame for comaring and later calculate correlations between sites
BB22.full.table <- t(as.data.frame.matrix(+(table(BB22.full.species[, c("site", "plant.species")]) > 0)))
BB22.full.table <- as.data.frame(BB22.full.table)
BB22.full.table$species <- rownames(BB22.full.table)

# create species lists per site
sitenames <- c("ZHUA", "ZHUB", "ZHUC", "ZHRD", "ZHRE", "ZHRF", "BEUA", "BEUB", "BEUC", "BERD", "BSUA", "BSUB", "BSUC", "BSRD", "BSRE", "BSRF")
site.species.bb <- list()

for (i in sitenames) {
  site.species.bb[[i]]  <- unique(BB22.full.species$plant.species[BB22.full.species$site == i])
} # end loop i

# import data from GBIF and InfoFlora and convert into species lists per site
# InfoFlora
#prepare list to store species lists
site.list.InfoFlora.species <- list()
for (i in sitenames) {
  # import data
  site.data.InfoFlora  <- read_csv(paste("./GBIF and InfoFlora/InfoFlora_", i, "_1500.csv", sep = ""))
  # summarize only needed variables (site and species)
  site.data.InfoFlora <- site.data.InfoFlora %>%
    summarise(Species = Taxon,
              site = rep(i, nrow(site.data.InfoFlora)),
              region = rep(substr(i, 1, 3), nrow(site.data.InfoFlora)))
  j <- 0
  site.data.InfoFlora$species <- c()
  for (j in 1:nrow(site.data.InfoFlora)){
    # species names are different than in GBIF -> need to be adapted
    site.data.InfoFlora$species[j] <- paste(unlist(strsplit(site.data.InfoFlora$Species[j], split=' ', fixed=TRUE))[1],
                                            unlist(strsplit(site.data.InfoFlora$Species[j], split=' ', fixed=TRUE))[2])
  } # end loop j
  
  # store species per site in a list
  site.list.InfoFlora.species[[i]] <- unique(site.data.InfoFlora$species)
  
} # end loop i

# GBIF
#prepare list to store species lists
site.list.gbif.species <- list()
for (i in sitenames) {
  # import data and filter for needed variables
  site.data.gbif  <- read_csv(paste("./GBIF and InfoFlora/GBIF_", i, "_1500.csv", sep = "")) %>%
    mutate(species = as_factor(species)) %>%
    filter(class == "Magnoliopsida" | class == "Liliopsida")
  site.list.gbif.species[[i]] <- unique(site.data.gbif$species)
} # end loop i

# combine species list GBIF and InfoFlora (site level) in new list
site.species.occ <- list()
for (i in sitenames) {
  x <- as.factor(site.list.InfoFlora.species[[i]])
  y <- as.factor(site.list.gbif.species[[i]])
  site.species.occ[[i]] <- unique(c(x,y))%>% 
    droplevels()
} # end loop i

# create list with intersection of species lists for all sites
intersection.sites.species <- list()
table.sites.species <- list()
correlation.site.species <- data.frame(species = NA, site1 = NA, site2 = NA, Correlation = NA)

# compute and store correlations between each site
for (i in sitenames) {
  for (j in sitenames) {
    # find common species of two sites
    intersection.sites.species[[paste(i,"/",j, sep="")]] <- intersect(site.species.occ[[i]], site.species.occ[[j]])
    
    # filter binary dataframe for the two sites
    table.sites.species[[paste(i,"/",j, sep="")]] <- dplyr::select(BB22.full.table, i, j)
    
    # compare bumblebee pollen with pollen lists of landscape
    comparison.table <- BB22.full.table[BB22.full.table$species %in% intersection.sites.species[[paste(i,"/",j, sep="")]], c(i,j)]
    
    # compute correlation and store them in a dataframe
    correlation <- c("B.lapidarius", i, j, round(cor(comparison.table[,1], comparison.table[,2], method=c("pearson")),3))
    correlation.site.species <- rbind(correlation.site.species, correlation)
    
  } # end loop j
} # end loop i

correlation.site.species$site1 <- as.factor(correlation.site.species$site1)
matrix.site.species <- split(correlation.site.species, f = correlation.site.species$site1)
temp <- list.cbind(matrix.site.species)[, seq(4, 67, 4)]
colnames(temp) <- levels(correlation.site.species$site1)
rownames(temp) <- c("ZHUA", "ZHUB", "ZHUC", "ZHRD", "ZHRE", "ZHRF", "BEUA", "BEUB", "BEUC", "BERD", "BSUA", "BSUB", "BSUC", "BSRD", "BSRE", "BSRF")
matrix.site.species <- as.matrix(temp[,match(sitenames, colnames(temp))])
class(matrix.site.species) <- "numeric"

# store in list
lapi.site[["species"]] <- abs(matrix.site.species)

#### genus ----
# convert into binary data frame
BB22.full.table <- t(as.data.frame.matrix(+(table(BB22.full.species[, c("site", "genus")]) > 0)))
BB22.full.table <- as.data.frame(BB22.full.table)
BB22.full.table$genus <- rownames(BB22.full.table)

# create genus lists per site
sitenames <- c("ZHUA", "ZHUB", "ZHUC", "ZHRD", "ZHRE", "ZHRF", "BEUA", "BEUB", "BEUC", "BERD", "BSUA", "BSUB", "BSUC", "BSRD", "BSRE", "BSRF")
site.genus.bb <- list()

for (i in sitenames) {
  site.genus.bb[[i]]  <- unique(BB22.full.species$genus[BB22.full.species$site == i])
} # end loop i

# import data from GBIF and InfoFlora and convert into species lists per site
# InfoFlora
library(stringr)
#prepare list to store species lists
site.list.InfoFlora.genus <- list()
for (i in sitenames) {
  # import data
  site.data.InfoFlora  <- read_csv(paste("./GBIF and InfoFlora/InfoFlora_", i, "_1500.csv", sep = ""))
  # summarize only needed variables (site and species)
  site.data.InfoFlora <- site.data.InfoFlora %>%
    summarise(genus = word(Taxon, 1),
              site = rep(i, nrow(site.data.InfoFlora)),
              region = rep(substr(i, 1, 3), nrow(site.data.InfoFlora)))
  
  site.list.InfoFlora.genus[[i]] <- unique(site.data.InfoFlora$genus)
} # end loop i

# GBIF
#prepare list to store species lists
site.list.gbif.genus <- list()
for (i in sitenames) {
  # import data and filter for needed variables
  site.data.gbif  <- read_csv(paste("./GBIF and InfoFlora/GBIF_", i, "_1500.csv", sep = "")) %>%
    mutate(genus = as_factor(genus)) %>%
    filter(class == "Magnoliopsida" | class == "Liliopsida")
  site.list.gbif.genus[[i]] <- unique(site.data.gbif$genus)
} # end loop i

# combine species list GBIF and InfoFlora (site level) in new list
site.genus.occ <- list()
for (i in sitenames) {
  x <- as.factor(site.list.InfoFlora.genus[[i]])
  y <- as.factor(site.list.gbif.genus[[i]])
  site.genus.occ[[i]] <- unique(c(x,y)) %>% 
    droplevels()
} # end loop i

# create list with intersection of genus lists for all sites
intersection.sites.genus <- list()
table.sites.genus <- list()
correlation.site.genus <- data.frame(genus = NA, site1 = NA, site2 = NA, Correlation = NA)

# compute and store correlations between each site
for (i in sitenames) {
  for (j in sitenames) {
    # find common genus of two sites
    intersection.sites.genus[[paste(i,"/",j, sep="")]] <- intersect(site.genus.occ[[i]], site.genus.occ[[j]])
    
    # filter binary dataframe for the two sites
    table.sites.genus[[paste(i,"/",j, sep="")]] <- dplyr::select(BB22.full.table, i, j)
    
    # compare bumblebee pollen with pollen lists of landscape
    comparison.table <- BB22.full.table[BB22.full.table$genus %in% intersection.sites.genus[[paste(i,"/",j, sep="")]], c(i,j)]
    
    # compute correlation and store them in a dataframe
    correlation <- c("B.lapidarius", i, j, round(cor(comparison.table[,1], comparison.table[,2], method=c("pearson")),3))
    correlation.site.genus <- rbind(correlation.site.genus, correlation)
    
  } # end loop j
} # end loop i

correlation.site.genus$site1 <- as.factor(correlation.site.genus$site1)
matrix.site.genus <- split(correlation.site.genus, f = correlation.site.genus$site1)
temp <- list.cbind(matrix.site.genus)[, seq(4, 67, 4)]
colnames(temp) <- levels(correlation.site.genus$site1)
rownames(temp) <- c("ZHUA", "ZHUB", "ZHUC", "ZHRD", "ZHRE", "ZHRF", "BEUA", "BEUB", "BEUC", "BERD", "BSUA", "BSUB", "BSUC", "BSRD", "BSRE", "BSRF")
matrix.site.genus <- as.matrix(temp[,match(sitenames, colnames(temp))])
class(matrix.site.genus) <- "numeric"

# store in list
lapi.site[["genus"]] <- abs(matrix.site.genus)


#### family ----

# convert into binary data frame
BB22.full.table <- t(as.data.frame.matrix(+(table(BB22.full.species[, c("site", "family")]) > 0)))
BB22.full.table <- as.data.frame(BB22.full.table)
BB22.full.table$family <- rownames(BB22.full.table)

# create genus lists per site
sitenames <- c("ZHUA", "ZHUB", "ZHUC", "ZHRD", "ZHRE", "ZHRF", "BEUA", "BEUB", "BEUC", "BERD", "BSUA", "BSUB", "BSUC", "BSRD", "BSRE", "BSRF")
site.family.bb <- list()

for (i in sitenames) {
  site.family.bb[[i]]  <- unique(BB22.full.species$family[BB22.full.species$site == i])
} # end loop i

# import data from GBIF and InfoFlora and convert into species lists per site
# InfoFlora
#prepare list to store species lists
site.list.InfoFlora.family <- list()
for (i in sitenames) {
  # import data
  site.data.InfoFlora  <- read_csv(paste("./GBIF and InfoFlora/InfoFlora_", i, "_1500.csv", sep = ""))
  # summarize only needed variables (site and species)
  site.data.InfoFlora <- site.data.InfoFlora %>%
    summarise(family = Familie,
              site = rep(i, nrow(site.data.InfoFlora)),
              region = rep(substr(i, 1, 3), nrow(site.data.InfoFlora)))
  
  site.list.InfoFlora.family[[i]] <- unique(site.data.InfoFlora$family)
} # end loop i

# GBIF
#prepare list to store species lists
site.list.gbif.family <- list()
for (i in sitenames) {
  # import data and filter for needed variables
  site.data.gbif  <- read_csv(paste("./GBIF and InfoFlora/GBIF_", i, "_1500.csv", sep = "")) %>%
    mutate(family = as_factor(family)) %>%
    filter(class == "Magnoliopsida" | class == "Liliopsida")
  site.list.gbif.family[[i]] <- unique(site.data.gbif$family)
} # end loop i

# combine species list GBIF and InfoFlora (site level) in new list
site.family.occ <- list()
for (i in sitenames) {
  x <- as.factor(site.list.InfoFlora.family[[i]])
  y <- as.factor(site.list.gbif.family[[i]])
  site.family.occ[[i]] <- unique(c(x,y)) %>% 
    droplevels()
} # end loop i

# create list with intersection of family lists for all sites
intersection.sites.family <- list()
table.sites.family <- list()
correlation.site.family <- data.frame(family = NA, site1 = NA, site2 = NA, Correlation = NA)

# compute and store correlations between each site
for (i in sitenames) {
  for (j in sitenames) {
    # find common family of two sites
    intersection.sites.family[[paste(i,"/",j, sep="")]] <- intersect(site.family.occ[[i]], site.family.occ[[j]])
    
    # filter binary dataframe for the two sites
    table.sites.family[[paste(i,"/",j, sep="")]] <- dplyr::select(BB22.full.table, i, j)
    
    # compare bumblebee pollen with pollen lists of landscape
    comparison.table <- BB22.full.table[BB22.full.table$family %in% intersection.sites.family[[paste(i,"/",j, sep="")]], c(i,j)]
    
    # compute correlation and store them in a dataframe
    correlation <- c("B.lapidarius", i, j, round(cor(comparison.table[,1], comparison.table[,2], method=c("pearson")),3))
    correlation.site.family <- rbind(correlation.site.family, correlation)
    
  } # end loop j
} # end loop i

correlation.site.family$site1 <- as.factor(correlation.site.family$site1)
matrix.site.family <- split(correlation.site.family, f = correlation.site.family$site1)
temp <- list.cbind(matrix.site.family)[, seq(4, 67, 4)]
colnames(temp) <- levels(correlation.site.family$site1)
rownames(temp) <- c("ZHUA", "ZHUB", "ZHUC", "ZHRD", "ZHRE", "ZHRF", "BEUA", "BEUB", "BEUC", "BERD", "BSUA", "BSUB", "BSUC", "BSRD", "BSRE", "BSRF")
matrix.site.family <- as.matrix(temp[,match(sitenames, colnames(temp))])
class(matrix.site.family) <- "numeric"

# store in list
lapi.site[["family"]] <- abs(matrix.site.family)

#### plotting ----
library(ggcorrplot)
axis.order.site <- c("BERD", 
                     "BSRD", "BSRE", "BSRF",
                     "ZHRD", "ZHRE", "ZHRF",
                     "BEUA", "BEUB", "BEUC",
                     "BSUA", "BSUB", "BSUC", 
                     "ZHUA", "ZHUB", "ZHUC") # to consistenly order the axis

c1 <- ggcorrplot(lapi.site[["species"]][axis.order.site, axis.order.site], 
                 hc.order = F, type = "lower", 
                 outline.col = "white")+ 
  labs(title ='species', fill = "", x = "", y = "") +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_gradient2(limit = c(0,1), low = "white", high =  "#07575B")

c2 <- ggcorrplot(lapi.site[["genus"]][axis.order.site, axis.order.site], 
                 hc.order = F, type = "lower", 
                 outline.col = "white")+ 
  labs(title ='genus', fill = "", x = "", y = "") +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_gradient2(limit = c(0,1), low = "white", high =  "#07575B")

c3 <- ggcorrplot(lapi.site[["family"]][axis.order.site, axis.order.site], 
                 hc.order = F, type = "lower", 
                 outline.col = "white")+ 
  labs(title ='family', fill = "", x = "", y = "") +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_gradient2(limit = c(0,1), low = "white", high =  "#07575B")

##### FIGURE S10 ----
# arrange them into one file to export
setwd(output)
plot3 <- ggarrange(c1,c2,c3, ncol = 1, nrow = 3,
                   labels = c("A", "B", "C"), common.legend = TRUE, legend = "right")
annotate_figure(plot3, top = text_grob("B.lapidarius", face = "bold", size = 22))
ggsave("./landscape comparison/diet consistency/corr_lapi_site.png", width = 10, height = 30)
setwd(input)

### region level ----
# prepare list element to store correlation matrices for species, genus and family
lapi.region <- list()

#### species ----

# convert into binary data frame for comaring and later calculate correlations between regions
BB22.full.table <- t(as.data.frame.matrix(+(table(BB22.full.species[, c("region", "plant.species")]) > 0)))
BB22.full.table <- as.data.frame(BB22.full.table)
BB22.full.table$species <- rownames(BB22.full.table)

# create species lists per region
regionnames <- c("ZHU", "ZHR", "BSU", "BSR", "BEU", "BER")
region.species.bb <- list()

for (i in regionnames) {
  region.species.bb[[i]]  <- unique(BB22.full.species$plant.species[BB22.full.species$region == i])
} # end loop i

# import data from GBIF and InfoFlora and convert into species lists per site
# InfoFlora
# combine species list per site into regional level
region.list.InfoFlora <- list()
region.list.InfoFlora[["ZHU"]] <- unique(c(site.list.InfoFlora.species$ZHUA, site.list.InfoFlora.species$ZHUB, site.list.InfoFlora.species$ZHUC))
region.list.InfoFlora[["ZHR"]] <- unique(c(site.list.InfoFlora.species$ZHRD, site.list.InfoFlora.species$ZHRE, site.list.InfoFlora.species$ZHRF))
region.list.InfoFlora[["BSU"]] <- unique(c(site.list.InfoFlora.species$BSUA, site.list.InfoFlora.species$BSUB, site.list.InfoFlora.species$BSUC))
region.list.InfoFlora[["BSR"]] <- unique(c(site.list.InfoFlora.species$BSRD, site.list.InfoFlora.species$BSRE, site.list.InfoFlora.species$BSRF))
region.list.InfoFlora[["BEU"]] <- unique(c(site.list.InfoFlora.species$BEUA, site.list.InfoFlora.species$BEUB, site.list.InfoFlora.species$BEUC))
region.list.InfoFlora[["BER"]] <- unique(c(site.list.InfoFlora.species$BERD))

# GBIF
# combine species list per site into regional level
region.list.gbif <- list()
region.list.gbif[["ZHU"]] <- unique(c(site.list.gbif.species$ZHUA, site.list.gbif.species$ZHUB, site.list.gbif.species$ZHUC))
region.list.gbif[["ZHR"]] <- unique(c(site.list.gbif.species$ZHRD, site.list.gbif.species$ZHRE, site.list.gbif.species$ZHRF))
region.list.gbif[["BSU"]] <- unique(c(site.list.gbif.species$BSUA, site.list.gbif.species$BSUB, site.list.gbif.species$BSUC))
region.list.gbif[["BSR"]] <- unique(c(site.list.gbif.species$BSRD, site.list.gbif.species$BSRE, site.list.gbif.species$BSRF))
region.list.gbif[["BEU"]] <- unique(c(site.list.gbif.species$BEUA, site.list.gbif.species$BEUB, site.list.gbif.species$BEUC))
region.list.gbif[["BER"]] <- unique(c(site.list.gbif.species$BERD))


# combine species list GBIF and InfoFlora (region level) in new list
region.species.occ <- list()
for (i in regionnames) {
  x <- as.factor(region.list.InfoFlora[[i]])
  y <- as.factor(region.list.gbif[[i]])
  region.species.occ[[i]] <- unique(c(x,y))%>% 
    droplevels()
} # end loop i

# create list with intersection of species lists for all regions
intersection.regions.species <- list()
table.regions.species <- list()
correlation.region.species <- data.frame(species = NA, region1 = NA, region2 = NA, Correlation = NA)

# compute and store correlations between each region
for (i in regionnames) {
  for (j in regionnames) {
    # find common species of two regions
    intersection.regions.species[[paste(i,"/",j, sep="")]] <- intersect(region.species.occ[[i]], region.species.occ[[j]])
    
    # filter binary dataframe for the two regions
    table.regions.species[[paste(i,"/",j, sep="")]] <- dplyr::select(BB22.full.table, i, j)
    
    # compare bumblebee pollen with pollen lists of landscape
    comparison.table <- BB22.full.table[BB22.full.table$species %in% intersection.regions.species[[paste(i,"/",j, sep="")]], c(i,j)]
    
    # compute correlation and store them in a dataframe
    correlation <- c("B.lapidarius", i, j, round(cor(comparison.table[,1], comparison.table[,2], method=c("pearson")),3))
    correlation.region.species <- rbind(correlation.region.species, correlation)
    
  } # end loop j
} # end loop i

correlation.region.species$region1 <- as.factor(correlation.region.species$region1)
matrix.region.species <- split(correlation.region.species, f = correlation.region.species$region1)
temp <- list.cbind(matrix.region.species)[, seq(4, 24, 4)]
colnames(temp) <- levels(correlation.region.species$region1)
rownames(temp) <- regionnames
matrix.region.species <- as.matrix(temp[,match(regionnames, colnames(temp))])
class(matrix.region.species) <- "numeric"

# store in list
lapi.region[["species"]] <- abs(matrix.region.species)


#### genus ----
# convert into binary data frame
BB22.full.table <- t(as.data.frame.matrix(+(table(BB22.full.species[, c("region", "genus")]) > 0)))
BB22.full.table <- as.data.frame(BB22.full.table)
BB22.full.table$genus <- rownames(BB22.full.table)

# create genus lists per region
region.genus.bb <- list()

for (i in regionnames) {
  region.genus.bb[[i]]  <- unique(BB22.full.species$genus[BB22.full.species$region == i])
} # end loop i

# import data from GBIF and InfoFlora and convert into species lists per site
# InfoFlora
# combine species list per site into regional level
region.list.InfoFlora <- list()
region.list.InfoFlora[["ZHU"]] <- unique(c(site.list.InfoFlora.genus$ZHUA, site.list.InfoFlora.genus$ZHUB, site.list.InfoFlora.genus$ZHUC))
region.list.InfoFlora[["ZHR"]] <- unique(c(site.list.InfoFlora.genus$ZHRD, site.list.InfoFlora.genus$ZHRE, site.list.InfoFlora.genus$ZHRF))
region.list.InfoFlora[["BSU"]] <- unique(c(site.list.InfoFlora.genus$BSUA, site.list.InfoFlora.genus$BSUB, site.list.InfoFlora.genus$BSUC))
region.list.InfoFlora[["BSR"]] <- unique(c(site.list.InfoFlora.genus$BSRD, site.list.InfoFlora.genus$BSRE, site.list.InfoFlora.genus$BSRF))
region.list.InfoFlora[["BEU"]] <- unique(c(site.list.InfoFlora.genus$BEUA, site.list.InfoFlora.genus$BEUB, site.list.InfoFlora.genus$BEUC))
region.list.InfoFlora[["BER"]] <- unique(c(site.list.InfoFlora.genus$BERD))

# GBIF
# combine species list per site into regional level
region.list.gbif <- list()
region.list.gbif[["ZHU"]] <- unique(c(site.list.gbif.genus$ZHUA, site.list.gbif.genus$ZHUB, site.list.gbif.genus$ZHUC))
region.list.gbif[["ZHR"]] <- unique(c(site.list.gbif.genus$ZHRD, site.list.gbif.genus$ZHRE, site.list.gbif.genus$ZHRF))
region.list.gbif[["BSU"]] <- unique(c(site.list.gbif.genus$BSUA, site.list.gbif.genus$BSUB, site.list.gbif.genus$BSUC))
region.list.gbif[["BSR"]] <- unique(c(site.list.gbif.genus$BSRD, site.list.gbif.genus$BSRE, site.list.gbif.genus$BSRF))
region.list.gbif[["BEU"]] <- unique(c(site.list.gbif.genus$BEUA, site.list.gbif.genus$BEUB, site.list.gbif.genus$BEUC))
region.list.gbif[["BER"]] <- unique(c(site.list.gbif.genus$BERD))

# combine species list GBIF and InfoFlora (region level) in new list
region.genus.occ <- list()
for (i in regionnames) {
  x <- as.factor(region.list.InfoFlora[[i]])
  y <- as.factor(region.list.gbif[[i]])
  region.genus.occ[[i]] <- unique(c(x,y)) %>% 
    droplevels()
} # end loop i

# create list with intersection of genus lists for all regions
intersection.regions.genus <- list()
table.regions.genus <- list()
correlation.region.genus <- data.frame(genus = NA, region1 = NA, region2 = NA, Correlation = NA)

# compute and store correlations between each region
for (i in regionnames) {
  for (j in regionnames) {
    # find common genus of two regions
    intersection.regions.genus[[paste(i,"/",j, sep="")]] <- intersect(region.genus.occ[[i]], region.genus.occ[[j]])
    
    # filter binary dataframe for the two regions
    table.regions.genus[[paste(i,"/",j, sep="")]] <- dplyr::select(BB22.full.table, i, j)
    
    # compare bumblebee pollen with pollen lists of landscape
    comparison.table <- BB22.full.table[BB22.full.table$genus %in% intersection.regions.genus[[paste(i,"/",j, sep="")]], c(i,j)]
    
    # compute correlation and store them in a dataframe
    correlation <- c("B.lapidarius", i, j, round(cor(comparison.table[,1], comparison.table[,2], method=c("pearson")),3))
    correlation.region.genus <- rbind(correlation.region.genus, correlation)
    
  } # end loop j
} # end loop i

correlation.region.genus$region1 <- as.factor(correlation.region.genus$region1)
matrix.region.genus <- split(correlation.region.genus, f = correlation.region.genus$region1)
temp <- list.cbind(matrix.region.genus)[, seq(4, 24, 4)]
colnames(temp) <- levels(correlation.region.genus$region1)
rownames(temp) <- regionnames
matrix.region.genus <- as.matrix(temp[,match(regionnames, colnames(temp))])
class(matrix.region.genus) <- "numeric"

# store in list
lapi.region[["genus"]] <- abs(matrix.region.genus)

#### family ----
# convert into binary data frame
BB22.full.table <- t(as.data.frame.matrix(+(table(BB22.full.species[, c("region", "family")]) > 0)))
BB22.full.table <- as.data.frame(BB22.full.table)
BB22.full.table$family <- rownames(BB22.full.table)

# create genus lists per site
region.family.bb <- list()

for (i in regionnames) {
  region.family.bb[[i]]  <- unique(BB22.full.species$family[BB22.full.species$region == i])
} # end loop i

# InfoFlora
# combine species list per site into regional level
region.list.InfoFlora <- list()
region.list.InfoFlora[["ZHU"]] <- unique(c(site.list.InfoFlora.family$ZHUA, site.list.InfoFlora.family$ZHUB, site.list.InfoFlora.family$ZHUC))
region.list.InfoFlora[["ZHR"]] <- unique(c(site.list.InfoFlora.family$ZHRD, site.list.InfoFlora.family$ZHRE, site.list.InfoFlora.family$ZHRF))
region.list.InfoFlora[["BSU"]] <- unique(c(site.list.InfoFlora.family$BSUA, site.list.InfoFlora.family$BSUB, site.list.InfoFlora.family$BSUC))
region.list.InfoFlora[["BSR"]] <- unique(c(site.list.InfoFlora.family$BSRD, site.list.InfoFlora.family$BSRE, site.list.InfoFlora.family$BSRF))
region.list.InfoFlora[["BEU"]] <- unique(c(site.list.InfoFlora.family$BEUA, site.list.InfoFlora.family$BEUB, site.list.InfoFlora.family$BEUC))
region.list.InfoFlora[["BER"]] <- unique(c(site.list.InfoFlora.family$BERD))

# GBIF
# combine species list per site into regional level
region.list.gbif <- list()
region.list.gbif[["ZHU"]] <- unique(c(site.list.gbif.family$ZHUA, site.list.gbif.family$ZHUB, site.list.gbif.family$ZHUC))
region.list.gbif[["ZHR"]] <- unique(c(site.list.gbif.family$ZHRD, site.list.gbif.family$ZHRE, site.list.gbif.family$ZHRF))
region.list.gbif[["BSU"]] <- unique(c(site.list.gbif.family$BSUA, site.list.gbif.family$BSUB, site.list.gbif.family$BSUC))
region.list.gbif[["BSR"]] <- unique(c(site.list.gbif.family$BSRD, site.list.gbif.family$BSRE, site.list.gbif.family$BSRF))
region.list.gbif[["BEU"]] <- unique(c(site.list.gbif.family$BEUA, site.list.gbif.family$BEUB, site.list.gbif.family$BEUC))
region.list.gbif[["BER"]] <- unique(c(site.list.gbif.family$BERD))

# combine species list GBIF and InfoFlora (region level) in new list
region.family.occ <- list()
for (i in regionnames) {
  x <- as.factor(region.list.InfoFlora[[i]])
  y <- as.factor(region.list.gbif[[i]])
  region.family.occ[[i]] <- unique(c(x,y)) %>% 
    droplevels()
} # end loop i

# create list with intersection of family lists for all regions
intersection.regions.family <- list()
table.regions.family <- list()
correlation.region.family <- data.frame(family = NA, region1 = NA, region2 = NA, Correlation = NA)

# compute and store correlations between each region
for (i in regionnames) {
  for (j in regionnames) {
    # find common family of two regions
    intersection.regions.family[[paste(i,"/",j, sep="")]] <- intersect(region.family.occ[[i]], region.family.occ[[j]])
    
    # filter binary dataframe for the two regions
    table.regions.family[[paste(i,"/",j, sep="")]] <- dplyr::select(BB22.full.table, i, j)
    
    # compare bumblebee pollen with pollen lists of landscape
    comparison.table <- BB22.full.table[BB22.full.table$family %in% intersection.regions.family[[paste(i,"/",j, sep="")]], c(i,j)]
    
    # compute correlation and store them in a dataframe
    correlation <- c("B.lapidarius", i, j, round(cor(comparison.table[,1], comparison.table[,2], method=c("pearson")),3))
    correlation.region.family <- rbind(correlation.region.family, correlation)
    
  } # end loop j
} # end loop i

correlation.region.family$region1 <- as.factor(correlation.region.family$region1)
matrix.region.family <- split(correlation.region.family, f = correlation.region.family$region1)
temp <- list.cbind(matrix.region.family)[, seq(4, 24, 4)]
colnames(temp) <- levels(correlation.region.family$region1)
rownames(temp) <- regionnames
matrix.region.family <- as.matrix(temp[,match(regionnames, colnames(temp))])
class(matrix.region.family) <- "numeric"

# store in list
lapi.region[["family"]] <- abs(matrix.region.family)

#### plotting ----
library(ggcorrplot)
axis.order.reg <- c("BER", "BSR", "ZHR", "BEU", "BSU", "ZHU") # to consistenly order the axis

d1 <- ggcorrplot(lapi.region[["species"]][axis.order.reg, axis.order.reg], 
                 hc.order = F, type = "lower", 
                 outline.col = "white")+ 
  labs(title ='species', fill = "", x = "", y = "") +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_gradient2(limit = c(0,1), low = "white", high =  "#07575B")


d2 <- ggcorrplot(lapi.region[["genus"]][axis.order.reg, axis.order.reg], 
                 hc.order = F, type = "lower", 
                 outline.col = "white")+ 
  labs(title ='genus', fill = "", x = "", y = "") +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_gradient2(limit = c(0,1), low = "white", high =  "#07575B")

d3 <- ggcorrplot(lapi.region[["family"]][axis.order.reg, axis.order.reg], 
                 hc.order = F, type = "lower", 
                 outline.col = "white")+ 
  labs(title ='family', fill = "", x = "", y = "") +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_gradient2(limit = c(0,1), low = "white", high =  "#07575B")

##### FIGURE 5 ----
# arrange them into one file to export
setwd(output)
plot4 <- ggarrange(d1,d2,d3, ncol = 1, nrow = 3,
                   labels = c("A", "B", "C"), common.legend = TRUE, legend = "right")
annotate_figure(plot4, top = text_grob("B.lapidarius", face = "bold", size = 22))
ggsave("./landscape comparison/diet consistency/corr_lapi_region.png", width = 10, height = 30)
setwd(input)

