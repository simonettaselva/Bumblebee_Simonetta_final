################################################
# Statistical Models script chemistry analysis
# by Simonetta Selva
#
# Created: May 2nd, 2023
# Project: Bumblebee 2022
################################################

# preparation ----
# clear work environment
rm(list=ls()) 

#load libraries
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(glmmTMB)
library(DHARMa)
library(vegan)

# load function
#function to produce model-checking plots for the fixed effects of an lmer model
fix.check <- function(mod){
  par(mfrow = c(1,3))
  plot(fitted(mod),resid(mod),main="Scale-location plot")	#should have no pattern
  abline(h = 0, col="red", lty=2)
  print(anova(lm(fitted(mod)~resid(mod))))	#should be non-significant
  qqnorm(resid(mod), ylab="Residuals")		#should be approximately straight line
  qqline(resid(mod), col="red")
  plot(density(resid(mod)))					#should be roughly normally distributed
  rug(resid(mod))
  par(mfrow = c(1,1))
}


# set working directory to main repository
input <- "YOUR PATH/Bumblebee_Simonetta_final/01_DATA" #change „YOUR PATH“ to where the repo is
output <- "YOUR PATH/Bumblebee_Simonetta_final/03_OUTPUT" #change „YOUR PATH“ to where the repo is
setwd(input)

# load data
BB22.bb.traits.input <- read_csv("BB22_traits.csv")
BB22.chemical.input <- read.csv("BB22_chemical.csv")

# summarize chemical data on site level
BB22.chemical <- BB22.chemical.input %>% 
  group_by(site, bbspecies) %>%
  summarise(location = as.factor(location),
            landscape = as.factor(landscape),
            replicate = as.factor(replicate),
            bbspecies = as.factor(bbspecies),
            site = as.factor(paste(location, landscape, replicate, sep="")),
            AS.mean = mean(AS)) %>%
  mutate(replicate = fct_relevel(replicate,"A", "B", "C", "D", "E", "F")) %>%
  distinct()

# summarize bb traits data on site level
BB22.bb.traits <- BB22.bb.traits.input %>% 
  mutate(landscape = as.factor(landscape))
BB22.bb.traits <- BB22.bb.traits %>%
  mutate(landscape = recode_factor(landscape, rural = "R", urban = "U"),
        site = as.factor(paste(location, landscape, replicate, sep="")))
BB22.bb.traits.mean <- BB22.bb.traits %>%
  group_by(site, bbspecies) %>%
  summarise(location = as.factor(location),
            landscape = as.factor(landscape),
            replicate = as.factor(replicate),
            bbspecies = as.factor(bbspecies),
            intertegular_distance.mean = mean(intertegular_distance),
            glossa.mean = mean(glossa),
            prementum.mean = mean(prementum),
            proboscis_length.mean = mean(proboscis_length),
            proboscis_ratio.mean = mean(proboscis_ratio),
            fore_wing_length.mean = mean(fore_wing_length),
            fore_wing_ratio.mean = mean(fore_wing_ratio),
            corbicula_length.mean = mean(corbicula_length),
            corbicula_ratio.mean = mean(corbicula_ratio))%>%
  distinct()

# # combine traits and cwm in new dataframes
# BB22.chem.site <- merge(BB22.chemical, BB22.bb.traits.mean[, c(1,6:14)], by  = "site", all.x=TRUE)

# import data with spatial information on sites
BB22.sites.meta <- read_csv("BB22_sites_2016.csv")

# B.PASCUORUM ----------------------------------------------------------------------------------------
#only B.pascuroum 
BB22.bb.sp <- BB22.bb.traits.mean[BB22.bb.traits.mean$bbspecies == "B.pascuorum",]
BB22.chemical.sp <- BB22.chemical[BB22.chemical$bbspecies == "B.pascuorum",]

# combine traits and cwm in new dataframes
BB22.chem.site <- merge(BB22.chemical.sp, BB22.bb.sp[, c(1,6:14)], by  = "site", all.x=TRUE) %>%
  mutate(region = paste(location, landscape, sep=""))

## look at data ----------------------------------------------------------------------------------------

# Boxplots for all the variables we want to look at with Wilcoxon test (landscape level)
library(psych)
library(rstatix)
palette.landscape <- c("#E69F00", "#56B4E9") #create color palette for landscape

w.test <- wilcox_test(BB22.chem.site,AS.mean~landscape)

setwd(output)
ggplot(BB22.chem.site, aes(x=landscape, y = AS.mean, fill=landscape)) + 
  geom_boxplot(notch = T) + 
  xlab("") +
  scale_x_discrete(labels=c('rural', 'urban'))+
  theme_classic(base_size = 20) +     
  theme(aspect.ratio=1) + 
  guides(alpha = "none") +
  scale_fill_manual(values=palette.landscape, guide = "none") + 
  labs(subtitle = paste("p = ", w.test$p, sep=""))
ggsave(paste("./chemistry/pasc_site/chemical_pasc_AS_landscapes.png", sep = ""), width = 8, height = 8)
setwd(input)


BB22.chem.site.comp <- BB22.chem.site %>% 
  filter(complete.cases(.)) %>% 
  droplevels()

# plot the relationship of plants traits of one site and bumblebee traits of one site
traits.bb <- colnames(BB22.chem.site.comp[, 7:15]) # bumblebee traits to look at; prepare for loop
plot_list <- list()

for (i in traits.bb) {
  plot_list[[i]] <- 
    ggplot(BB22.chem.site.comp, aes_string(i, "AS.mean", colour = "landscape")) + 
    geom_point() + 
    theme_classic(base_size = 20) + 
    theme(aspect.ratio=1) + 
    geom_smooth(method="lm", se = FALSE) +
    scale_color_manual(values=palette.landscape, labels=c("rural", "urban")) + 
    stat_cor(aes(color = landscape), size = 5)
  
} # end loop i

setwd(output)
plot <- ggarrange(plot_list[[1]],plot_list[[2]],plot_list[[3]],
                  plot_list[[4]],plot_list[[5]],plot_list[[6]],
                  plot_list[[7]],plot_list[[8]],plot_list[[9]],
                  ncol = 3, nrow = 3, 
                  labels = c(LETTERS[1:8]),   
                  common.legend = TRUE)

annotate_figure(plot, top = text_grob(paste("B.pascuorum: comparison of AS and traits across landscapes", sep = ""),
                                       face = "bold", size = 22))
ggsave(paste("./chemistry/pasc_site/chemical_pasc_AS_BBtraits_landscapes.png", sep = ""), width = 10, height = 10)
setwd(input)


## modelling ----------------------------------------------------------------------------------------
# center and scale all the variables
BB22.chem.site.comp[, 6:14] <- scale(BB22.chem.site.comp[, 6:14],center=TRUE,scale=TRUE)

# correlation analysis
# look at the correlation between the explanatory variables
library(corrplot)
M<-cor(BB22.chem.site.comp[, 7:14], use = "complete.obs") # subset data; only explanatory variables

# matrix of the p-value of the correlation
p.mat <- cor.mtest(BB22.chem.site.comp[, 7:14])
head(p.mat)

#correlation plot
corrplot::corrplot(M, type="upper", order="hclust", 
                   p.mat = p.mat$p, sig.level = 0.01, tl.col = "black",
                   col = COL2('RdBu', 10)) # plot correlation with p-values

# !!! a lot are multicollinear. in the model only use: proboscis_ratio, fore_wing_ratio

#  fit GLMM
#load the libraries
library(lme4)
library(car)
library(MuMIn)
library(arm)

# built an initial full model based on collinearity 
M1.full <- lmer(AS.mean ~ proboscis_ratio.mean + fore_wing_ratio.mean + (1|landscape),
                data = BB22.chem.site.comp)
fix.check(M1.full) # looks ok
vif(M1.full) # looks good

# subset meta information on the sites based on sites, we have AS data on
BB22.sites.meta.sub <- subset(BB22.sites.meta, site %in% BB22.chem.site.comp$site)

# formal test for spatial correlation
sims <- simulateResiduals(M1.full)
BB22.chem.site.comp$site <- as.factor(BB22.chem.site.comp$site)
simulationOutput <- recalculateResiduals(sims, group = BB22.chem.site.comp$site)
testSpatialAutocorrelation(simulationOutput,
                           x = BB22.sites.meta.sub$LV95_x, 
                           y = BB22.sites.meta.sub$LV95_y, 
                           plot = FALSE)
# there is no spatial autocorrelation

# update model: add intertegular distance to the model
M1.full.1 <- lmer(AS.mean ~ intertegular_distance.mean + proboscis_ratio.mean + fore_wing_ratio.mean + 
                    (1|landscape),
                  data = BB22.chem.site.comp)
fix.check(M1.full.1)
summary(M1.full.1)
vif(M1.full.1) # looks good
anova(M1.full.1, M1.full) # better fit without intertegular

# dredging
# Things to take into account when dredging:
# 1) epends on the models included in the candidate set. You can’t identify a model as being the 
# “best” fit to the data if you didn’t include the model to begin with!
# 2) The parameter estimates and predictions arising from the “best” model or set of best models 
# should be biologically meaningful.

options(na.action = "na.fail") # Required for dredge to run
std.model <- MuMIn::dredge(M1.full.1)
options(na.action = "na.omit") # set back to default
# Get the top best models
top.mod <- get.models(std.model, subset = delta < 6) ## Delta-AICc < 6 (Burnham et al., 2011)
# Model averaging
avg.model <- model.avg(top.mod,revised.var = TRUE)
summary(avg.model)
avg.model$sw

#  no variable is significant


## relationship FDs & species richness vs. AS ----------------------------------------------------------------------------------------

# import FDs on site level for B.pascuorum
setwd(input)
BB22.fd.sites <- read.csv("./FD/FD_package_B.pascuorum_site.csv") %>% 
  rename(site = X) # ajust variable name of site

# merge AS and FD into one dataframe
BB22.fd.chem.sites <- merge(BB22.chemical.sp[, 1:6], BB22.fd.sites, by  = "site", all.x=TRUE) %>%
  mutate(region = paste(location, landscape, sep="")) # add variable region

 a <- ggplot(BB22.fd.chem.sites, aes(nbsp.w, AS.mean)) + 
  geom_point() + 
  ylim(25, 150) +
  theme_classic(base_size = 20) + 
  theme(aspect.ratio=1) + 
  xlab("species richness per site") +
  ylab("Aminoacid content (µg/mg)") +
  geom_smooth(method="lm", se = FALSE) +
  stat_cor(size = 5)

b <- ggplot(BB22.fd.chem.sites, aes(FRic.w, AS.mean)) + 
  geom_point() + 
  ylim(25, 150) +
  theme_classic(base_size = 20) + 
  theme(aspect.ratio=1) + 
  xlab("functional richness per site") +
  ylab("Aminoacid content (µg/mg)") +
  geom_smooth(method="lm", se = FALSE) +
  stat_cor(size = 5)

c <- ggplot(BB22.fd.chem.sites, aes(FEve.w, AS.mean)) + 
  geom_point() + 
  ylim(25, 150) +
  theme_classic(base_size = 20) + 
  theme(aspect.ratio=1) + 
  xlab("functional evenness per site") +
  ylab("Aminoacid content (µg/mg)") +
  geom_smooth(method="lm", se = FALSE) +
  stat_cor(size = 5)

d <- ggplot(BB22.fd.chem.sites, aes(FDiv.w, AS.mean)) + 
  geom_point() + 
  ylim(25, 150) +
  theme_classic(base_size = 20) + 
  theme(aspect.ratio=1) + 
  xlab("functional divergence per site") +
  ylab("Aminoacid content (µg/mg)") +
  geom_smooth(method="lm", se = FALSE) +
  stat_cor(size = 5)

setwd(output)
plot <- ggarrange(a, b, c, d, 
                  ncol = 2, nrow = 2, 
                  labels = c(LETTERS[1:4]),   
                  common.legend = TRUE)

annotate_figure(plot, top = text_grob(paste("B.pascuorum: comparison of AS and FDs", sep = ""),
                                      face = "bold", size = 22))
ggsave(paste("./chemistry/pasc_site/chemical_pasc_AS_BBtraits.png", sep = ""), width = 10, height = 10)
setwd(input)

# B.LAPIDARIUS ----------------------------------------------------------------------------------------
#only B.lapidarius
BB22.bb.sp <- BB22.bb.traits.mean[BB22.bb.traits.mean$bbspecies == "B.lapidarius",]
BB22.chemical.sp <- BB22.chemical[BB22.chemical$bbspecies == "B.lapidarius",]

# combine traits and cwm in new dataframes
BB22.chem.site <- merge(BB22.chemical.sp, BB22.bb.sp[, c(1,6:14)], by  = "site", all.x=TRUE) %>%
  mutate(region = paste(location, landscape, sep=""))

## look at data ----------------------------------------------------------------------------------------

# Boxplots for all the variables we want to look at with Wilcoxon test (landscape level)
library(psych)
library(rstatix)
palette.landscape <- c("#E69F00", "#56B4E9") #create color palette for landscape

w.test <- wilcox_test(BB22.chem.site,AS.mean~landscape)

setwd(output)
ggplot(BB22.chem.site, aes(x=landscape, y = AS.mean, fill=landscape)) + 
  geom_boxplot(notch = T) + 
  xlab("") +
  scale_x_discrete(labels=c('rural', 'urban'))+
  theme_classic(base_size = 20) +     
  theme(aspect.ratio=1) + 
  guides(alpha = "none") +
  scale_fill_manual(values=palette.landscape, guide = "none") + 
  labs(subtitle = paste("p = ", w.test$p, sep=""))
ggsave(paste("./chemistry/lapi_site/chemical_lapi_AS_landscapes.png", sep = ""), width = 8, height = 8)
setwd(input)


BB22.chem.site.comp <- BB22.chem.site %>% 
  filter(complete.cases(.)) %>% 
  droplevels()

# plot the relationship of plants traits of one site and bumblebee traits of one site
traits.bb <- colnames(BB22.chem.site.comp[, 7:15]) # bumblebee traits to look at; prepare for loop
plot_list <- list()

for (i in traits.bb) {
  plot_list[[i]] <- 
    ggplot(BB22.chem.site.comp, aes_string(i, "AS.mean", colour = "landscape")) + 
    geom_point() + 
    theme_classic(base_size = 20) + 
    theme(aspect.ratio=1) + 
    geom_smooth(method="lm", se = FALSE) +
    scale_color_manual(values=palette.landscape, labels=c("rural", "urban")) + 
    stat_cor(aes(color = landscape), size = 5)
  
} # end loop i

setwd(output)
plot <- ggarrange(plot_list[[1]],plot_list[[2]],plot_list[[3]],
                  plot_list[[4]],plot_list[[5]],plot_list[[6]],
                  plot_list[[7]],plot_list[[8]],plot_list[[9]],
                  ncol = 3, nrow = 3, 
                  labels = c(LETTERS[1:8]),   
                  common.legend = TRUE)

annotate_figure(plot, top = text_grob(paste("B.lapidarius: comparison of AS and traits across landscapes", sep = ""),
                                      face = "bold", size = 22))
ggsave(paste("./chemistry/lapi_site/chemical_lapi_AS_BBtraits_landscapes.png", sep = ""), width = 10, height = 10)
setwd(input)


## modelling ----------------------------------------------------------------------------------------
# center and scale all the variables
BB22.chem.site.comp[, 6:14] <- scale(BB22.chem.site.comp[, 6:14],center=TRUE,scale=TRUE)

# correlation analysis
# look at the correlation between the explanatory variables
library(corrplot)
M<-cor(BB22.chem.site.comp[, 7:14], use = "complete.obs") # subset data; only explanatory variables

# matrix of the p-value of the correlation
p.mat <- cor.mtest(BB22.chem.site.comp[, 7:14])
head(p.mat)

#correlation plot
corrplot::corrplot(M, type="upper", order="hclust",
                   p.mat = p.mat$p, sig.level = 0.01, tl.col = "black",
                   col = COL2('RdBu', 10)) # plot correlation with p-values

# !!! a lot are multicollinear. in the model only use: proboscis_ratio, fore_wing_ratio

#  fit GLMM
#load the libraries
library(lme4)
library(car)
library(MuMIn)
library(arm)

# built an initial full model based on collinearity
M1.full <- lmer(AS.mean ~ proboscis_ratio.mean + fore_wing_ratio.mean + (1|landscape),
                data = BB22.chem.site.comp)
fix.check(M1.full) # looks ok
vif(M1.full) # looks good

# subset meta information on the sites based on sites, we have AS data on
BB22.sites.meta.sub <- subset(BB22.sites.meta, site %in% BB22.chem.site.comp$site)

# formal test for spatial correlation
sims <- simulateResiduals(M1.full)
BB22.chem.site.comp$site <- as.factor(BB22.chem.site.comp$site)
simulationOutput <- recalculateResiduals(sims, group = BB22.chem.site.comp$site)
testSpatialAutocorrelation(simulationOutput,
                           x = BB22.sites.meta.sub$LV95_x,
                           y = BB22.sites.meta.sub$LV95_y,
                           plot = FALSE)
# there is no spatial autocorrelation

# update model: add intertegular distance to the model
M1.full.1 <- lmer(AS.mean ~ intertegular_distance.mean + proboscis_ratio.mean + fore_wing_ratio.mean +
                    (1|landscape),
                  data = BB22.chem.site.comp)
fix.check(M1.full.1)
summary(M1.full.1)
vif(M1.full.1) # looks good
anova(M1.full.1, M1.full) # better fit without intertegular

# dredging
# Things to take into account when dredging:
# 1) epends on the models included in the candidate set. You can’t identify a model as being the
# “best” fit to the data if you didn’t include the model to begin with!
# 2) The parameter estimates and predictions arising from the “best” model or set of best models
# should be biologically meaningful.

options(na.action = "na.fail") # Required for dredge to run
std.model <- MuMIn::dredge(M1.full.1)
options(na.action = "na.omit") # set back to default
# Get the top best models
top.mod <- get.models(std.model, subset = delta < 6) ## Delta-AICc < 6 (Burnham et al., 2011)
# Model averaging
avg.model <- model.avg(top.mod,revised.var = TRUE)
summary(avg.model)
avg.model$sw

#  no variable is significant


## relationship FDs & species richness vs. AS ----------------------------------------------------------------------------------------

# import FDs on site level for B.lapidarius
setwd(input)
BB22.fd.sites <- read.csv("./FD/FD_package_B.lapidarius_site.csv") %>% 
  rename(site = X) # ajust variable name of site

# merge AS and FD into one dataframe
BB22.fd.chem.sites <- merge(BB22.chemical.sp[, 1:6], BB22.fd.sites, by  = "site", all.x=TRUE) %>%
  mutate(region = paste(location, landscape, sep="")) # add variable region

a <- ggplot(BB22.fd.chem.sites, aes(nbsp.w, AS.mean)) + 
  geom_point() + 
  ylim(25, 150) +
  theme_classic(base_size = 20) + 
  theme(aspect.ratio=1) + 
  xlab("species richness per site") +
  ylab("Aminoacid content (µg/mg)") +
  geom_smooth(method="lm", se = FALSE) +
  stat_cor(size = 5)

b <- ggplot(BB22.fd.chem.sites, aes(FRic.w, AS.mean)) + 
  geom_point() + 
  ylim(25, 150) +
  theme_classic(base_size = 20) + 
  theme(aspect.ratio=1) + 
  xlab("functional richness per site") +
  ylab("Aminoacid content (µg/mg)") +
  geom_smooth(method="lm", se = FALSE) +
  stat_cor(size = 5)

c <- ggplot(BB22.fd.chem.sites, aes(FEve.w, AS.mean)) + 
  geom_point() + 
  ylim(25, 150) +
  theme_classic(base_size = 20) + 
  theme(aspect.ratio=1) + 
  xlab("functional evenness per site") +
  ylab("Aminoacid content (µg/mg)") +
  geom_smooth(method="lm", se = FALSE) +
  stat_cor(size = 5)

d <- ggplot(BB22.fd.chem.sites, aes(FDiv.w, AS.mean)) + 
  geom_point() + 
  ylim(25, 150) +
  theme_classic(base_size = 20) + 
  theme(aspect.ratio=1) + 
  xlab("functional divergence per site") +
  ylab("Aminoacid content (µg/mg)") +
  geom_smooth(method="lm", se = FALSE) +
  stat_cor(size = 5)

setwd(output)
plot <- ggarrange(a, b, c, d, 
                  ncol = 2, nrow = 2, 
                  labels = c(LETTERS[1:4]),   
                  common.legend = TRUE)

annotate_figure(plot, top = text_grob(paste("B.lapidarius: comparison of AS and FDs", sep = ""),
                                      face = "bold", size = 22))
ggsave(paste("./chemistry/lapi_site/chemical_lapi_AS_BBtraits.png", sep = ""), width = 10, height = 10)
setwd(input)
