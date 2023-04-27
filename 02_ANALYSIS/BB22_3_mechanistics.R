################################################
# Mechansims Shaping the Diet
# by Simonetta Selva
#
# Created: May 2nd, 2023
# Project: Bumblebee 2022
################################################

# HYPOTHESIS: Bumblebees with more generalistic feeding habits, characterised 
# by shorter tongue length, will visit a more diverse range of plant species

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
library(DHARMa)

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
  title(main=print(paste0(summary(mod)$call[2])),outer=T, line=-1)
  par(mfrow = c(1,1))
}

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

### modelling ----------------------------------------------------------------------------------------
# center and scale all the variables
BB22.ID[, c(7:15, 17:20)] <- scale(BB22.ID[, c(7:15, 17:20)],center=TRUE,scale=TRUE)

# correlation analysis
# look at the correlation between the explanatory variables
library(corrplot)
M <- cor(BB22.ID[, 7:15], use = "complete.obs") # subset data; only explanatory variables

# matrix of the p-value of the correlation
p.mat <- cor.mtest(BB22.ID[, 7:15])
head(p.mat)

#### FIGURE S4 ----
# correlation plot
corrplot::corrplot(M, type="upper", order="hclust", 
                   p.mat = p.mat$p, sig.level = 0.01, tl.col = "black",
                   col = COL2('RdBu', 10)) # plot correlation with p-values

# !!! a lot are multicollinear. in the model only use: proboscis_ratio, fore_wing_ratio and corbicula_ratio


#  fit GLMM
#load the libraries
library(lme4)
library(car)
library(MuMIn)
library(arm)

#### Species Richness ----------------------------------------------------------------------------------------

# built an initial full model based on collinearity 
M1.full <- lmer(sp_richn ~ proboscis_ratio + fore_wing_ratio + corbicula_ratio + 
                  (1|landscape),
                data=BB22.ID)
summary(M1.full)

# check model assumptions
par(mar = c(8, 2, 5, 2))
fix.check(M1.full) # looks ok
vif <- vif(M1.full)
vif # looks good
mtext(side = 1, line = 6, adj = 0,
      paste("VIF:",
            paste(paste(names(vif)[1], round(vif[1], 3), sep=": "),
            paste(names(vif)[2], round(vif[2], 3), sep=": "),
            paste(names(vif)[3], round(vif[3], 3), sep=": "),
            sep = ", "), sep = " ")
      )
      
# formal test for spatial correlation
sims <- simulateResiduals(M1.full)
BB22.ID$site <- as.factor(BB22.ID$site)
simulationOutput <- recalculateResiduals(sims, group = BB22.ID$site)
testSpatialAutocorrelation(simulationOutput, x = BB22.sites.meta$LV95_x, y = BB22.sites.meta$LV95_y, plot = FALSE)
# there is no spatial autocorrelation

# update model: add intertegular distance to the model
M1.full.1 <- lmer(sp_richn ~ intertegular_distance + proboscis_ratio + fore_wing_ratio + 
                    corbicula_ratio + (1|landscape),
                  data=BB22.ID)
summary(M1.full.1)

# check model assumptions
par(mar = c(8, 2, 5, 2))
fix.check(M1.full.1) # looks ok
vif <- vif(M1.full.1)
vif # looks good
mtext(side = 1, line = 6, adj = 0,
      paste("VIF:",
            paste(paste(names(vif)[1], round(vif[1], 3), sep=": "),
                  paste(names(vif)[2], round(vif[2], 3), sep=": "),
                  paste(names(vif)[3], round(vif[3], 3), sep=": "),
                  paste(names(vif)[4], round(vif[4], 3), sep=": "),
                  sep = ", "), sep = " ")
)

# check which models fits better
anova(M1.full.1, M1.full) # with intertegular_distance better fit than M1.full

# dredging
# Things to take into account when dredging:
# 1) depends on the models included in the candidate set. You can’t identify a model as being the 
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

# the only variable having a significant effect in any model is intertegular distance

#### Functional Richness ----------------------------------------------------------------------------------------
# built an initial full model based on colinearity 
M2.full <- lmer(fric ~ proboscis_ratio + fore_wing_ratio + corbicula_ratio + 
                  (1|landscape),
                data=BB22.ID)
summary(M2.full)

# check model assumptions
par(mar = c(8, 2, 5, 2))
fix.check(M2.full) # looks ok
vif <- vif(M2.full)
vif # looks good
mtext(side = 1, line = 6, adj = 0,
      paste("VIF:",
            paste(paste(names(vif)[1], round(vif[1], 3), sep=": "),
                  paste(names(vif)[2], round(vif[2], 3), sep=": "),
                  paste(names(vif)[3], round(vif[3], 3), sep=": "),
                  sep = ", "), sep = " ")
)

# formal test for spatial correlation
sims <- simulateResiduals(M2.full)
BB22.ID$site <- as.factor(BB22.ID$site)
simulationOutput <- recalculateResiduals(sims, group = BB22.ID$site)
testSpatialAutocorrelation(simulationOutput, x = BB22.sites.meta$LV95_x, y = BB22.sites.meta$LV95_y, plot = FALSE)
# there is no spatial autocorrelation

# add intertegular distance to the model
M2.full.1 <- lmer(fric ~ intertegular_distance + proboscis_ratio + fore_wing_ratio + 
                    corbicula_ratio + (1|landscape),
                  data=BB22.ID) # boundary (singular) fit: see help('isSingular')
summary(M2.full.1)

# check model assumptions
par(mar = c(8, 2, 5, 2))
fix.check(M2.full.1) # looks ok
vif <- vif(M2.full.1)
vif # looks good
mtext(side = 1, line = 6, adj = 0,
      paste("VIF:",
            paste(paste(names(vif)[1], round(vif[1], 3), sep=": "),
                  paste(names(vif)[2], round(vif[2], 3), sep=": "),
                  paste(names(vif)[3], round(vif[3], 3), sep=": "),
                  paste(names(vif)[4], round(vif[4], 3), sep=": "),
                  sep = ", "), sep = " ")
)

# check which models fits better
anova(M2.full.1, M2.full) # with intertegular_distance better fit than M1.full

# dredging
options(na.action = "na.fail") # Required for dredge to run
std.model <- MuMIn::dredge(M2.full.1)
options(na.action = "na.omit") # set back to default
# Get the top best models
top.mod <- get.models(std.model, subset = delta < 6) ## Delta-AICc < 6 (Burnham et al., 2011)
# Model averaging
avg.model <- model.avg(top.mod,revised.var = TRUE)
summary(avg.model)
avg.model$sw

#### Functional Divergence ----------------------------------------------------------------------------------------
# built an initial full model based on collinearity 
M3.full <- lmer(fdiv ~ proboscis_ratio + fore_wing_ratio + corbicula_ratio + (1|landscape),
                data=BB22.ID)
summary(M3.full)

# check model assumptions
par(mar = c(8, 2, 5, 2))
fix.check(M3.full) # looks ok
vif <- vif(M3.full)
vif # looks good
mtext(side = 1, line = 6, adj = 0,
      paste("VIF:",
            paste(paste(names(vif)[1], round(vif[1], 3), sep=": "),
                  paste(names(vif)[2], round(vif[2], 3), sep=": "),
                  paste(names(vif)[3], round(vif[3], 3), sep=": "),
                  sep = ", "), sep = " ")
)

# formal test for spatial correlation
sims <- simulateResiduals(M3.full)
BB22.ID$site <- as.factor(BB22.ID$site)
simulationOutput <- recalculateResiduals(sims, group = BB22.ID$site)
testSpatialAutocorrelation(simulationOutput, x = BB22.sites.meta$LV95_x, y = BB22.sites.meta$LV95_y, plot = FALSE)
# there is no spatial autocorrelation

# add intertegular distance to the model
M3.full.1 <- lmer(fdiv ~ intertegular_distance + proboscis_ratio + fore_wing_ratio + 
                    corbicula_ratio + (1|landscape),
                  data=BB22.ID)
summary(M3.full.1)

# check model assumptions
par(mar = c(8, 2, 5, 2))
fix.check(M3.full.1) # looks ok
vif <- vif(M3.full.1)
vif # looks good
mtext(side = 1, line = 6, adj = 0,
      paste("VIF:",
            paste(paste(names(vif)[1], round(vif[1], 3), sep=": "),
                  paste(names(vif)[2], round(vif[2], 3), sep=": "),
                  paste(names(vif)[3], round(vif[3], 3), sep=": "),
                  paste(names(vif)[4], round(vif[4], 3), sep=": "),
                  sep = ", "), sep = " ")
)

# check which models fits better
anova(M3.full.1, M3.full) # with intertegular_distance not a better fit than M1.full

# dredging
options(na.action = "na.fail") # Required for dredge to run
std.model <- MuMIn::dredge(M3.full.1)
options(na.action = "na.omit") # set back to default
# Get the top best models
top.mod <- get.models(std.model, subset = delta < 6) ## Delta-AICc < 6 (Burnham et al., 2011)
# Model averaging
avg.model <- model.avg(top.mod,revised.var = TRUE)
summary(avg.model)
avg.model$sw

#### Functional Evenness ----------------------------------------------------------------------------------------
# built an initial full model based on collinearity 
M4.full <- lmer(feve ~ proboscis_ratio + fore_wing_ratio + corbicula_ratio + (1|landscape),
                data=BB22.ID) #boundary (singular) fit: see help('isSingular')
summary(M4.full)

# check model assumptions
par(mar = c(8, 2, 5, 2))
fix.check(M4.full) # looks ok
vif <- vif(M4.full)
vif # looks good
mtext(side = 1, line = 6, adj = 0,
      paste("VIF:",
            paste(paste(names(vif)[1], round(vif[1], 3), sep=": "),
                  paste(names(vif)[2], round(vif[2], 3), sep=": "),
                  paste(names(vif)[3], round(vif[3], 3), sep=": "),
                  sep = ", "), sep = " ")
)

# formal test for spatial correlation
sims <- simulateResiduals(M4.full)
BB22.ID$site <- as.factor(BB22.ID$site)
simulationOutput = recalculateResiduals(sims, group = BB22.ID$site)
testSpatialAutocorrelation(sims, x = BB22.sites.meta$LV95_x, y = BB22.sites.meta$LV95_y, plot = FALSE)
# there is no spatial autocorrelation

# add intertegular distance to the model
M4.full.1 <- lmer(feve ~ intertegular_distance + proboscis_ratio + fore_wing_ratio + 
                    corbicula_ratio + (1|landscape),
                  data=BB22.ID) # boundary (singular) fit: see help('isSingular')
summary(M4.full.1)

# check model assumptions
par(mar = c(8, 2, 5, 2))
fix.check(M4.full.1) # looks ok
vif <- vif(M4.full.1)
vif # looks good
mtext(side = 1, line = 6, adj = 0,
      paste("VIF:",
            paste(paste(names(vif)[1], round(vif[1], 3), sep=": "),
                  paste(names(vif)[2], round(vif[2], 3), sep=": "),
                  paste(names(vif)[3], round(vif[3], 3), sep=": "),
                  paste(names(vif)[4], round(vif[4], 3), sep=": "),
                  sep = ", "), sep = " ")
)

# check which models fits better
anova(M4.full.1, M4.full) # with intertegular_distance not a better fit than M1.full

# dredging
options(na.action = "na.fail") # Required for dredge to run
std.model <- MuMIn::dredge(M4.full.1)
options(na.action = "na.omit") # set back to default
# Get the top best models
top.mod <- get.models(std.model, subset = delta < 6) ## Delta-AICc < 6 (Burnham et al., 2011)
# Model averaging
avg.model <- model.avg(top.mod,revised.var = TRUE)
summary(avg.model)
avg.model$sw

# save the last eight produced plots in a selected path
plots.dir.path <- list.files(tempdir(), pattern="rs-graphics", full.names = TRUE);
plots.png.paths <- list.files(plots.dir.path, pattern=".png", full.names = TRUE)
file.copy(from = plots.png.paths,
          to= paste(output, "/mechanistics/functional diversity/pasc_ID/models", sep = ""))

#clear all plots
dev.off(dev.list()["RStudioGD"])

### Model testing ----------------------------------------------------------------------------------------

#### Data splitting  ----------------------------------------------------------------------------------------
library(rsample)
set.seed(123) # for reproducibility
BB22.ID.complete <- BB22.ID %>%
  filter(complete.cases(.))
split <- initial_split(BB22.ID.complete, prop = 0.8)
ddf_train <- training(split)
ddf_test <- testing(split) 

#### Model training and predicting ----------------------------------------------------------------------------------------
library(caret) 
library(yardstick)

##### Linear Model ----------------------------------------------------------------------------------------
metrics <- colnames(BB22.ID[, c(17:20)]) # plant FD to look at (see file BB22_compute_FD); prepare for loop

for (i in metrics) {
  f <- as.formula(paste(i,  "~ intertegular_distance + proboscis_ratio + fore_wing_ratio + corbicula_ratio + landscape", sep = ""))
  train <- train(form = f,
                 data = ddf_train, 
                 method = "lm",
                 trControl = trainControl(method = "cv", number = 10))
  
  ## Get variable importance, and turn into a data frame
  var_imp <- varImp(train, scale=FALSE)$importance
  var_imp <- data.frame(variables=row.names(var_imp), importance=var_imp$Overall)
  
  ## Create a plot of variable importance
  var_imp %>%
    
    ## Sort the data by importance
    arrange(importance) %>%
    
    ## Create a ggplot object for aesthetic
    ggplot(aes(x=reorder(variables, importance), y=importance)) + 
    
    ## Plot the bar graph
    geom_bar(stat='identity') + 
    
    ## Flip the graph to make a horizontal bar plot
    coord_flip() + 
    
    ## Add x-axis label
    xlab('Variables') +
    
    ## Add a title
    labs(title='Linear model variable importance') + 
    
    ## Some layout for the plot
    theme_minimal() + 
    theme(axis.text = element_text(size = 10), 
          axis.title = element_text(size = 15), 
          plot.title = element_text(size = 20), 
    )
  setwd(output)
  ggsave(paste("./mechanistics/functional diversity/pasc_ID/training and predicting/linear model/", i, "_lm_var_imp_pasc.png", sep = ""), width = 8, height = 5)
  setwd(input)
  
  predict_train <- predict(
    ## lm object
    object=train, 
    ## Data to use for predictions
    newdata=ddf_train)
  
  predict_test <- predict(
    ## lm object
    object=train, 
    ## Data to use for predictions
    newdata=ddf_test)
  
  data_metrics_train <- data.frame(truth=ddf_train[,i], pred=predict_train) #mistake
  metrics_train <- metrics(data_metrics_train, truth, pred)
  
  data_metrics_test <- data.frame(truth=ddf_test[,i], pred=predict_test)
  metrics_test <- metrics(data_metrics_test, truth, pred)
  
  # plotting prediction
  library(patchwork)
  
  gg_test <- ggplot(data_metrics_test, aes(x=truth, y=pred))+
    geom_point()+
    labs(title = "Density of Data Point for Testing Data",
         subtitle = paste("Linear Model: Rsq = ", round(metrics_test$.estimate[2], 3), sep = ""), 
         x = paste("predicted values for ", i, sep = ""),
         y = paste("observed values for ", i, sep = "")) +
    theme(aspect.ratio=1)+
    geom_abline(intercept = 0, slope = 1, color = "#fc5e03")
  
  gg_train <- ggplot(data_metrics_train, aes(x=truth, y=pred))+
    geom_point()+
    labs(title = "Density of Data Point for Training Data",
         subtitle = paste("Linear Model: Rsq = ", round(metrics_train$.estimate[2], 3), sep = ""),
         x = paste("predicted values for ", i, sep = ""),
         y = paste("observed values for ", i, sep = "")) +
    theme(aspect.ratio=1)+
    geom_abline(intercept = 0, slope = 1, color = "#fc5e03")
  
  # arrange them into one file to export
  setwd(output)
  ggarrange(gg_test, gg_train, ncol = 2, nrow = 1,
            labels = c("A", "B"), common.legend = TRUE, legend = "right")
  ggsave(paste("./mechanistics/functional diversity/pasc_ID/training and predicting/linear model/", i, "_lm_prediction_pasc.png", sep = ""), width = 8, height = 5)
  setwd(input)
} # end i loop

##### Random Forest ----------------------------------------------------------------------------------------
for (i in metrics) {
  f <- as.formula(paste(i,  "~ intertegular_distance + proboscis_ratio + fore_wing_ratio + corbicula_ratio + landscape", sep = ""))
  train <- train(form = f,
                 data = ddf_train, 
                 method = "rf",
                 trControl = trainControl(method = "cv", number = 10))
  
  ## Get variable importance, and turn into a data frame
  var_imp <- varImp(train, scale=FALSE)$importance
  var_imp <- data.frame(variables=row.names(var_imp), importance=var_imp$Overall)
  
  ## Create a plot of variable importance
  var_imp %>%
    
    ## Sort the data by importance
    arrange(importance) %>%
    
    ## Create a ggplot object for aesthetic
    ggplot(aes(x=reorder(variables, importance), y=importance)) + 
    
    ## Plot the bar graph
    geom_bar(stat='identity') + 
    
    ## Flip the graph to make a horizontal bar plot
    coord_flip() + 
    
    ## Add x-axis label
    xlab('Variables') +
    
    ## Add a title
    labs(title='Random forest variable importance') + 
    
    ## Some layout for the plot
    theme_minimal() + 
    theme(axis.text = element_text(size = 10), 
          axis.title = element_text(size = 15), 
          plot.title = element_text(size = 20), 
    )
  setwd(output)
  ggsave(paste("./mechanistics/functional diversity/pasc_ID/training and predicting/random forest/", i, "_rf_var_imp_pasc.png", sep = ""), width = 8, height = 8)
  setwd(input)
  
  predict_train <- predict(
    ## lm object
    object=train, 
    ## Data to use for predictions
    newdata=ddf_train)
  
  predict_test <- predict(
    ## lm object
    object=train, 
    ## Data to use for predictions
    newdata=ddf_test)
  
  data_metrics_train <- data.frame(truth=ddf_train[, i], pred=predict_train)
  metrics_train <- metrics(data_metrics_train, truth, pred)
  
  data_metrics_test <- data.frame(truth=ddf_test[, i], pred=predict_test)
  metrics_test <- metrics(data_metrics_test, truth, pred)
  
  # plotting prediction
  library(patchwork)
  
  gg_test <- ggplot(data_metrics_test, aes(x=truth, y=pred))+
    geom_point()+
    labs(title = "Density of Data Point for Testing Data",
         subtitle = paste("Random Forest: Rsq = ", round(metrics_test$.estimate[2], 3), sep = ""),
         x = paste("predicted values for ", i, sep = ""),
         y = paste("observed values for ", i, sep = "")) +
    theme(aspect.ratio=1)+
    geom_abline(intercept = 0, slope = 1, color = "#fc5e03")
  
  gg_train <- ggplot(data_metrics_train, aes(x=truth, y=pred))+
    geom_point()+ 
    labs(title = "Density of Data Point for Training Data",
         subtitle = paste("Random Forest: Rsq = ", round(metrics_train$.estimate[2], 3), sep = ""),
         x = paste("predicted values for ", i, sep = ""),
         y = paste("observed values for ", i, sep = "")) +
    theme(aspect.ratio=1)+
    geom_abline(intercept = 0, slope = 1, color = "#fc5e03")
  
  # arrange them into one file to export
  setwd(output)
  ggarrange(gg_test, gg_train, ncol = 2, nrow = 1,
            labels = c("A", "B"), common.legend = TRUE, legend = "right")
  ggsave(paste("./mechanistics/functional diversity/pasc_ID/training and predicting/random forest/", i, "_rf_prediction_pasc.png", sep = ""), width = 8, height = 5)
  setwd(input)
  
} # end i loop

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
BB22.ID <- merge(BB22.bb.traits.sp, BB22.fun.ID, 
                 by  = "ID", 
                 all.x=TRUE)
BB22.ID <- merge(BB22.ID, BB22.sites.meta[, c(1,2,3)], 
                 by  = "site", 
                 all.x=TRUE) 

### modelling ----------------------------------------------------------------------------------------
# center and scale all the variables
BB22.ID[, c(7:15, 17:20)] <- scale(BB22.ID[, c(7:15, 17:20)],center=TRUE,scale=TRUE)

# correlation analysis
# look at the correlation between the explanatory variables
library(corrplot)
M <- cor(BB22.ID[, 7:15], use = "complete.obs") # subset data; only explanatory variables

# matrix of the p-value of the correlation
p.mat <- cor.mtest(BB22.ID[, 7:15])
head(p.mat)

#### FIGURE S4 ----
#correlation plot
corrplot::corrplot(M, type="upper", order="hclust", 
                   p.mat = p.mat$p, sig.level = 0.01, tl.col = "black",
                   col = COL2('RdBu', 10)) # plot correlation with p-values

# !!! a lot are multicollinear. in the model only use: proboscis_ratio, fore_wing_ratio and corbicula_ratio


#  fit GLMM
#load the libraries
library(lme4)
library(car)
library(MuMIn)
library(arm)

#### Species Richness ----------------------------------------------------------------------------------------

# built an initial full model based on collinearity 
M1.full <- lmer(sp_richn ~ proboscis_ratio + fore_wing_ratio + corbicula_ratio + 
                  (1|landscape),
                data=BB22.ID)
summary(M1.full)

# check model assumptions
par(mar = c(8, 2, 5, 2))
fix.check(M1.full) # looks ok
vif <- vif(M1.full)
vif # looks good
mtext(side = 1, line = 6, adj = 0,
      paste("VIF:",
            paste(paste(names(vif)[1], round(vif[1], 3), sep=": "),
                  paste(names(vif)[2], round(vif[2], 3), sep=": "),
                  paste(names(vif)[3], round(vif[3], 3), sep=": "),
                  sep = ", "), sep = " ")
)

# formal test for spatial correlation
sims <- simulateResiduals(M1.full)
BB22.ID$site <- as.factor(BB22.ID$site)
simulationOutput <- recalculateResiduals(sims, group = BB22.ID$site)
testSpatialAutocorrelation(simulationOutput, x = BB22.sites.meta$LV95_x, y = BB22.sites.meta$LV95_y, plot = FALSE)
# there is no spatial autocorrelation

# update model: add intertegular distance to the model
M1.full.1 <- lmer(sp_richn ~ intertegular_distance + proboscis_ratio + fore_wing_ratio + 
                    corbicula_ratio + (1|landscape),
                  data=BB22.ID)
summary(M1.full.1)

# check model assumptions
par(mar = c(8, 2, 5, 2))
fix.check(M1.full.1) # looks ok
vif <- vif(M1.full.1)
vif # looks good
mtext(side = 1, line = 6, adj = 0,
      paste("VIF:",
            paste(paste(names(vif)[1], round(vif[1], 3), sep=": "),
                  paste(names(vif)[2], round(vif[2], 3), sep=": "),
                  paste(names(vif)[3], round(vif[3], 3), sep=": "),
                  paste(names(vif)[4], round(vif[4], 3), sep=": "),
                  sep = ", "), sep = " ")
)

# check which models fits better
anova(M1.full.1, M1.full) # with intertegular_distance better fit than M1.full

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

# the only variable having a significant effect in any model is intertegular distance

#### Functional Richness ----------------------------------------------------------------------------------------
# built an initial full model based on colinearity 
M2.full <- lmer(fric ~ proboscis_ratio + fore_wing_ratio + corbicula_ratio + 
                  (1|landscape),
                data=BB22.ID)
summary(M2.full)

# check model assumptions
par(mar = c(8, 2, 5, 2))
fix.check(M2.full) # looks ok
vif <- vif(M2.full)
vif # looks good
mtext(side = 1, line = 6, adj = 0,
      paste("VIF:",
            paste(paste(names(vif)[1], round(vif[1], 3), sep=": "),
                  paste(names(vif)[2], round(vif[2], 3), sep=": "),
                  paste(names(vif)[3], round(vif[3], 3), sep=": "),
                  sep = ", "), sep = " ")
)

# formal test for spatial correlation
sims <- simulateResiduals(M2.full)
BB22.ID$site <- as.factor(BB22.ID$site)
simulationOutput <- recalculateResiduals(sims, group = BB22.ID$site)
testSpatialAutocorrelation(simulationOutput, x = BB22.sites.meta$LV95_x, y = BB22.sites.meta$LV95_y, plot = FALSE)
# there is no spatial autocorrelation

# add intertegular distance to the model
M2.full.1 <- lmer(fric ~ intertegular_distance + proboscis_ratio + fore_wing_ratio + 
                    corbicula_ratio + (1|landscape),
                  data=BB22.ID) # boundary (singular) fit: see help('isSingular')
summary(M2.full.1)

# check model assumptions
par(mar = c(8, 2, 5, 2))
fix.check(M2.full.1) # looks ok
vif <- vif(M2.full.1)
vif # looks good
mtext(side = 1, line = 6, adj = 0,
      paste("VIF:",
            paste(paste(names(vif)[1], round(vif[1], 3), sep=": "),
                  paste(names(vif)[2], round(vif[2], 3), sep=": "),
                  paste(names(vif)[3], round(vif[3], 3), sep=": "),
                  paste(names(vif)[4], round(vif[4], 3), sep=": "),
                  sep = ", "), sep = " ")
)

# check which models fits better
anova(M2.full.1, M2.full) # with intertegular_distance better fit than M1.full

# dredging
options(na.action = "na.fail") # Required for dredge to run
std.model <- MuMIn::dredge(M2.full.1)
options(na.action = "na.omit") # set back to default
# Get the top best models
top.mod <- get.models(std.model, subset = delta < 6) ## Delta-AICc < 6 (Burnham et al., 2011)
# Model averaging
avg.model <- model.avg(top.mod,revised.var = TRUE)
summary(avg.model)
avg.model$sw

#### Functional Divergence ----------------------------------------------------------------------------------------
# built an initial full model based on collinearity 
M3.full <- lmer(fdiv ~ proboscis_ratio + fore_wing_ratio + corbicula_ratio + (1|landscape),
                data=BB22.ID)
summary(M3.full)

# check model assumptions
par(mar = c(8, 2, 5, 2))
fix.check(M3.full) # looks ok
vif <- vif(M3.full)
vif # looks good
mtext(side = 1, line = 6, adj = 0,
      paste("VIF:",
            paste(paste(names(vif)[1], round(vif[1], 3), sep=": "),
                  paste(names(vif)[2], round(vif[2], 3), sep=": "),
                  paste(names(vif)[3], round(vif[3], 3), sep=": "),
                  sep = ", "), sep = " ")
)

# formal test for spatial correlation
sims <- simulateResiduals(M3.full)
BB22.ID$site <- as.factor(BB22.ID$site)
simulationOutput <- recalculateResiduals(sims, group = BB22.ID$site)
testSpatialAutocorrelation(simulationOutput, x = BB22.sites.meta$LV95_x, y = BB22.sites.meta$LV95_y, plot = FALSE)
# there is no spatial autocorrelation

# add intertegular distance to the model
M3.full.1 <- lmer(fdiv ~ intertegular_distance + proboscis_ratio + fore_wing_ratio + 
                    corbicula_ratio + (1|landscape),
                  data=BB22.ID)
summary(M3.full.1)

# check model assumptions
par(mar = c(8, 2, 5, 2))
fix.check(M3.full.1) # looks ok
vif <- vif(M3.full.1)
vif # looks good
mtext(side = 1, line = 6, adj = 0,
      paste("VIF:",
            paste(paste(names(vif)[1], round(vif[1], 3), sep=": "),
                  paste(names(vif)[2], round(vif[2], 3), sep=": "),
                  paste(names(vif)[3], round(vif[3], 3), sep=": "),
                  paste(names(vif)[4], round(vif[4], 3), sep=": "),
                  sep = ", "), sep = " ")
)

# check which models fits better
anova(M3.full.1, M3.full) # with intertegular_distance not a better fit than M1.full

# dredging
options(na.action = "na.fail") # Required for dredge to run
std.model <- MuMIn::dredge(M3.full.1)
options(na.action = "na.omit") # set back to default
# Get the top best models
top.mod <- get.models(std.model, subset = delta < 6) ## Delta-AICc < 6 (Burnham et al., 2011)
# Model averaging
avg.model <- model.avg(top.mod,revised.var = TRUE)
summary(avg.model)
avg.model$sw

#### Functional Evenness ----------------------------------------------------------------------------------------
# built an initial full model based on collinearity 
M4.full <- lmer(feve ~ proboscis_ratio + fore_wing_ratio + corbicula_ratio + (1|landscape),
                data=BB22.ID) #boundary (singular) fit: see help('isSingular')
summary(M4.full)

# check model assumptions
par(mar = c(8, 2, 5, 2))
fix.check(M4.full) # looks ok
vif <- vif(M4.full)
vif # looks good
mtext(side = 1, line = 6, adj = 0,
      paste("VIF:",
            paste(paste(names(vif)[1], round(vif[1], 3), sep=": "),
                  paste(names(vif)[2], round(vif[2], 3), sep=": "),
                  paste(names(vif)[3], round(vif[3], 3), sep=": "),
                  sep = ", "), sep = " ")
)

# formal test for spatial correlation
sims <- simulateResiduals(M4.full)
BB22.ID$site <- as.factor(BB22.ID$site)
simulationOutput = recalculateResiduals(sims, group = BB22.ID$site)
testSpatialAutocorrelation(sims, x = BB22.sites.meta$LV95_x, y = BB22.sites.meta$LV95_y, plot = FALSE)
# there is no spatial autocorrelation

# add intertegular distance to the model
M4.full.1 <- lmer(feve ~ intertegular_distance + proboscis_ratio + fore_wing_ratio + 
                    corbicula_ratio + (1|landscape),
                  data=BB22.ID) # boundary (singular) fit: see help('isSingular')
summary(M4.full.1)

# check model assumptions
par(mar = c(8, 2, 5, 2))
fix.check(M4.full.1) # looks ok
vif <- vif(M4.full.1)
vif # looks good
mtext(side = 1, line = 6, adj = 0,
      paste("VIF:",
            paste(paste(names(vif)[1], round(vif[1], 3), sep=": "),
                  paste(names(vif)[2], round(vif[2], 3), sep=": "),
                  paste(names(vif)[3], round(vif[3], 3), sep=": "),
                  paste(names(vif)[4], round(vif[4], 3), sep=": "),
                  sep = ", "), sep = " ")
)

# check which models fits better
anova(M4.full.1, M4.full) # with intertegular_distance not a better fit than M1.full

# dredging
options(na.action = "na.fail") # Required for dredge to run
std.model <- MuMIn::dredge(M4.full.1)
options(na.action = "na.omit") # set back to default
# Get the top best models
top.mod <- get.models(std.model, subset = delta < 6) ## Delta-AICc < 6 (Burnham et al., 2011)
# Model averaging
avg.model <- model.avg(top.mod,revised.var = TRUE)
summary(avg.model)
avg.model$sw


# save the last eight produced plots in a selected path
plots.dir.path <- list.files(tempdir(), pattern="rs-graphics", full.names = TRUE);
plots.png.paths <- list.files(plots.dir.path, pattern=".png", full.names = TRUE)
file.copy(from = plots.png.paths,
          to= paste(output, "/mechanistics/functional diversity/lapi_ID/models", sep = ""))

### Model testing ----------------------------------------------------------------------------------------

#### Data splitting  ----------------------------------------------------------------------------------------
library(rsample)
set.seed(123) # for reproducibility
BB22.ID.complete <- BB22.ID %>%
  filter(complete.cases(.))
split <- initial_split(BB22.ID.complete, prop = 0.8)
ddf_train <- training(split)
ddf_test <- testing(split) 

#### Model training and predicting ----------------------------------------------------------------------------------------
library(caret) 
library(yardstick)

##### Linear Model ----------------------------------------------------------------------------------------
metrics <- colnames(BB22.ID[, c(17:20)]) # plant FD to look at (see file BB22_compute_FD); prepare for loop

for (i in metrics) {
  f <- as.formula(paste(i,  "~ intertegular_distance + proboscis_ratio + fore_wing_ratio + corbicula_ratio + landscape", sep = ""))
  train <- train(form = f,
                 data = ddf_train, 
                 method = "lm",
                 trControl = trainControl(method = "cv", number = 10))
  
  ## Get variable importance, and turn into a data frame
  var_imp <- varImp(train, scale=FALSE)$importance
  var_imp <- data.frame(variables=row.names(var_imp), importance=var_imp$Overall)
  
  ## Create a plot of variable importance
  var_imp %>%
    
    ## Sort the data by importance
    arrange(importance) %>%
    
    ## Create a ggplot object for aesthetic
    ggplot(aes(x=reorder(variables, importance), y=importance)) + 
    
    ## Plot the bar graph
    geom_bar(stat='identity') + 
    
    ## Flip the graph to make a horizontal bar plot
    coord_flip() + 
    
    ## Add x-axis label
    xlab('Variables') +
    
    ## Add a title
    labs(title='Linear model variable importance') + 
    
    ## Some layout for the plot
    theme_minimal() + 
    theme(axis.text = element_text(size = 10), 
          axis.title = element_text(size = 15), 
          plot.title = element_text(size = 20), 
    )
  setwd(output)
  ggsave(paste("./mechanistics/functional diversity/lapi_ID/training and predicting/linear model/", i, "_lm_var_imp_lapi.png", sep = ""), width = 8, height = 5)
  setwd(input)
  
  predict_train <- predict(
    ## lm object
    object=train, 
    ## Data to use for predictions
    newdata=ddf_train)
  
  predict_test <- predict(
    ## lm object
    object=train, 
    ## Data to use for predictions
    newdata=ddf_test)
  
  data_metrics_train <- data.frame(truth=ddf_train[,i], pred=predict_train) #mistake
  metrics_train <- metrics(data_metrics_train, truth, pred)
  
  data_metrics_test <- data.frame(truth=ddf_test[,i], pred=predict_test)
  metrics_test <- metrics(data_metrics_test, truth, pred)
  
  # plotting prediction
  library(patchwork)
  
  gg_test <- ggplot(data_metrics_test, aes(x=truth, y=pred))+
    geom_point()+
    labs(title = "Density of Data Point for Testing Data",
         subtitle = paste("Linear Model: Rsq = ", round(metrics_test$.estimate[2], 3), sep = ""), 
         x = paste("predicted values for ", i, sep = ""),
         y = paste("observed values for ", i, sep = "")) +
    theme(aspect.ratio=1)+
    geom_abline(intercept = 0, slope = 1, color = "#fc5e03")
  
  gg_train <- ggplot(data_metrics_train, aes(x=truth, y=pred))+
    geom_point()+
    labs(title = "Density of Data Point for Training Data",
         subtitle = paste("Linear Model: Rsq = ", round(metrics_train$.estimate[2], 3), sep = ""),
         x = paste("predicted values for ", i, sep = ""),
         y = paste("observed values for ", i, sep = "")) +
    theme(aspect.ratio=1)+
    geom_abline(intercept = 0, slope = 1, color = "#fc5e03")
  
  # arrange them into one file to export
  setwd(output)
  ggarrange(gg_test, gg_train, ncol = 2, nrow = 1,
            labels = c("A", "B"), common.legend = TRUE, legend = "right")
  ggsave(paste("./mechanistics/functional diversity/lapi_ID/training and predicting/linear model/", i, "_lm_prediction_lapi.png", sep = ""), width = 8, height = 5)
  setwd(input)
} # end i loop

##### Random Forest ----------------------------------------------------------------------------------------
for (i in metrics) {
  f <- as.formula(paste(i,  "~ intertegular_distance + proboscis_ratio + fore_wing_ratio + corbicula_ratio + landscape", sep = ""))
  train <- train(form = f,
                 data = ddf_train, 
                 method = "rf",
                 trControl = trainControl(method = "cv", number = 10))
  
  ## Get variable importance, and turn into a data frame
  var_imp <- varImp(train, scale=FALSE)$importance
  var_imp <- data.frame(variables=row.names(var_imp), importance=var_imp$Overall)
  
  ## Create a plot of variable importance
  var_imp %>%
    
    ## Sort the data by importance
    arrange(importance) %>%
    
    ## Create a ggplot object for aesthetic
    ggplot(aes(x=reorder(variables, importance), y=importance)) + 
    
    ## Plot the bar graph
    geom_bar(stat='identity') + 
    
    ## Flip the graph to make a horizontal bar plot
    coord_flip() + 
    
    ## Add x-axis label
    xlab('Variables') +
    
    ## Add a title
    labs(title='Random forest variable importance') + 
    
    ## Some layout for the plot
    theme_minimal() + 
    theme(axis.text = element_text(size = 10), 
          axis.title = element_text(size = 15), 
          plot.title = element_text(size = 20), 
    )
  setwd(output)
  ggsave(paste("./mechanistics/functional diversity/lapi_ID/training and predicting/random forest/", i, "_rf_var_imp_lapi.png", sep = ""), width = 8, height = 8)
  setwd(input)
  
  predict_train <- predict(
    ## lm object
    object=train, 
    ## Data to use for predictions
    newdata=ddf_train)
  
  predict_test <- predict(
    ## lm object
    object=train, 
    ## Data to use for predictions
    newdata=ddf_test)
  
  data_metrics_train <- data.frame(truth=ddf_train[, i], pred=predict_train)
  metrics_train <- metrics(data_metrics_train, truth, pred)
  
  data_metrics_test <- data.frame(truth=ddf_test[, i], pred=predict_test)
  metrics_test <- metrics(data_metrics_test, truth, pred)
  
  # plotting prediction
  library(patchwork)
  
  gg_test <- ggplot(data_metrics_test, aes(x=truth, y=pred))+
    geom_point()+
    labs(title = "Density of Data Point for Testing Data",
         subtitle = paste("Random Forest: Rsq = ", round(metrics_test$.estimate[2], 3), sep = ""),
         x = paste("predicted values for ", i, sep = ""),
         y = paste("observed values for ", i, sep = "")) +
    theme(aspect.ratio=1)+
    geom_abline(intercept = 0, slope = 1, color = "#fc5e03")
  
  gg_train <- ggplot(data_metrics_train, aes(x=truth, y=pred))+
    geom_point()+ 
    labs(title = "Density of Data Point for Training Data",
         subtitle = paste("Random Forest: Rsq = ", round(metrics_train$.estimate[2], 3), sep = ""),
         x = paste("predicted values for ", i, sep = ""),
         y = paste("observed values for ", i, sep = "")) +
    theme(aspect.ratio=1)+
    geom_abline(intercept = 0, slope = 1, color = "#fc5e03")
  
  # arrange them into one file to export
  setwd(output)
  ggarrange(gg_test, gg_train, ncol = 2, nrow = 1,
            labels = c("A", "B"), common.legend = TRUE, legend = "right")
  ggsave(paste("./mechanistics/functional diversity/lapi_ID/training and predicting/random forest/", i, "_rf_prediction_lapi.png", sep = ""), width = 8, height = 5)
  setwd(input)
  
} # end i loop

# PHYLOGENETICS ----
## preparation ----
# clear work environment
rm(list=ls()) 

#load libraries
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(glmmTMB)
library(DHARMa)

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
  title(main=print(paste0(summary(mod)$call[2])),outer=T, line=-1)
  par(mfrow = c(1,1))
}

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

### modelling ----------------------------------------------------------------------------------------
# center and scale all the variables
BB22.ID[, c(7:15, 18:21)] <- scale(BB22.ID[, c(7:15, 17:20)],center=TRUE,scale=TRUE)
BB22.ID <- BB22.ID %>% filter(complete.cases(.)) # remove all entries with NAs

# correlation analysis
# look at the correlation between the explanatory variables
library(corrplot)
M <- cor(BB22.ID[, 7:15], use = "complete.obs") # subset data; only explanatory variables

# matrix of the p-value of the correlation
p.mat <- cor.mtest(BB22.ID[, 7:15])
head(p.mat)

#correlation plot
corrplot::corrplot(M, type="upper", order="hclust", 
                   p.mat = p.mat$p, sig.level = 0.01, tl.col = "black",
                   col = COL2('RdBu', 10)) # plot correlation with p-values

# !!! a lot are multicollinear. in the model only use: proboscis_ratio, fore_wing_ratio and corbicula_ratio


#  fit GLMM
#load the libraries
library(lme4)
library(car)
library(MuMIn)
library(arm)

#### Phylogenetic Species Variability  ----------------------------------------------------------------------------------------

# built an initial full model based on collinearity 
M1.full <- lmer(pvar ~ proboscis_ratio + fore_wing_ratio + corbicula_ratio + 
                  (1|landscape),
                data=BB22.ID)
summary(M1.full)

# check model assumptions
par(mar = c(8, 2, 5, 2))
fix.check(M1.full) # looks ok
vif <- vif(M1.full)
vif # looks good
mtext(side = 1, line = 6, adj = 0,
      paste("VIF:",
            paste(paste(names(vif)[1], round(vif[1], 3), sep=": "),
                  paste(names(vif)[2], round(vif[2], 3), sep=": "),
                  paste(names(vif)[3], round(vif[3], 3), sep=": "),
                  sep = ", "), sep = " ")
)

# formal test for spatial correlation
sims <- simulateResiduals(M1.full)
BB22.ID$site <- as.factor(BB22.ID$site)
simulationOutput <- recalculateResiduals(sims, group = BB22.ID$site)
testSpatialAutocorrelation(simulationOutput, x = BB22.sites.meta$LV95_x, y = BB22.sites.meta$LV95_y, plot = FALSE)
# there is no spatial autocorrelation

# update model: add intertegular distance to the model
M1.full.1 <- lmer(pvar ~ intertegular_distance + proboscis_ratio + fore_wing_ratio + 
                    corbicula_ratio + (1|landscape),
                  data=BB22.ID)
summary(M1.full.1)

# check model assumptions
par(mar = c(8, 2, 5, 2))
fix.check(M1.full.1) # looks ok
vif <- vif(M1.full.1)
vif # looks good
mtext(side = 1, line = 6, adj = 0,
      paste("VIF:",
            paste(paste(names(vif)[1], round(vif[1], 3), sep=": "),
                  paste(names(vif)[2], round(vif[2], 3), sep=": "),
                  paste(names(vif)[3], round(vif[3], 3), sep=": "),
                  paste(names(vif)[4], round(vif[4], 3), sep=": "),
                  sep = ", "), sep = " ")
)

# check which models fits better
anova(M1.full.1, M1.full) # without intertegular_distance better fit than M1.full

# dredging
# Things to take into account when dredging:
# 1) depends on the models included in the candidate set. You can’t identify a model as being the 
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

# no variable is having a significant effect in any model

# save the last eight produced plots in a selected path
# plots.dir.path <- list.files(tempdir(), pattern="rs-graphics", full.names = TRUE);
# plots.png.paths <- list.files(plots.dir.path, pattern=".png", full.names = TRUE)
# file.copy(from = plots.png.paths,
#           to= paste(output, "/mechanistics/phylogenetic diversity/pasc_ID/models", sep = ""))

## Model testing ----------------------------------------------------------------------------------------

### Data splitting  ----------------------------------------------------------------------------------------
library(rsample)
set.seed(123) # for reproducibility
BB22.ID.complete <- BB22.ID %>%
  filter(complete.cases(.))
split <- initial_split(BB22.ID.complete, prop = 0.8)
ddf_train <- training(split)
ddf_test <- testing(split) 

### Model training and predicting ----------------------------------------------------------------------------------------
library(caret) 
library(yardstick)

#### Linear Model ----------------------------------------------------------------------------------------
metrics <- "pvar"
for (i in metrics) {
  f <- as.formula(paste(i,  "~ intertegular_distance + proboscis_ratio + fore_wing_ratio + corbicula_ratio + landscape", sep = ""))
  train <- train(form = f,
                 data = ddf_train, 
                 method = "lm",
                 trControl = trainControl(method = "cv", number = 10))
  
  ## Get variable importance, and turn into a data frame
  var_imp <- varImp(train, scale=FALSE)$importance
  var_imp <- data.frame(variables=row.names(var_imp), importance=var_imp$Overall)
  
  ## Create a plot of variable importance
  var_imp %>%
    
    ## Sort the data by importance
    arrange(importance) %>%
    
    ## Create a ggplot object for aesthetic
    ggplot(aes(x=reorder(variables, importance), y=importance)) + 
    
    ## Plot the bar graph
    geom_bar(stat='identity') + 
    
    ## Flip the graph to make a horizontal bar plot
    coord_flip() + 
    
    ## Add x-axis label
    xlab('Variables') +
    
    ## Add a title
    labs(title='Linear model variable importance') + 
    
    ## Some layout for the plot
    theme_minimal() + 
    theme(axis.text = element_text(size = 10), 
          axis.title = element_text(size = 15), 
          plot.title = element_text(size = 20), 
    )
  setwd(output)
  ggsave(paste("./mechanistics/phylogenetic diversity/pasc_ID/training and predicting/linear model/", i, "_lm_var_imp_pasc.png", sep = ""), width = 8, height = 5)
  setwd(input)
  
  predict_train <- predict(
    ## lm object
    object=train, 
    ## Data to use for predictions
    newdata=ddf_train)
  
  predict_test <- predict(
    ## lm object
    object=train, 
    ## Data to use for predictions
    newdata=ddf_test)
  
  data_metrics_train <- data.frame(truth=ddf_train[,i], pred=predict_train) #mistake
  metrics_train <- metrics(data_metrics_train, truth, pred)
  
  data_metrics_test <- data.frame(truth=ddf_test[,i], pred=predict_test)
  metrics_test <- metrics(data_metrics_test, truth, pred)
  
  # plotting prediction
  library(patchwork)
  
  gg_test <- ggplot(data_metrics_test, aes(x=truth, y=pred))+
    geom_point()+
    labs(title = "Density of Data Point for Testing Data",
         subtitle = paste("Linear Model: Rsq = ", round(metrics_test$.estimate[2], 3), sep = ""), 
         x = paste("predicted values for ", i, sep = ""),
         y = paste("observed values for ", i, sep = "")) +
    theme(aspect.ratio=1)+
    geom_abline(intercept = 0, slope = 1, color = "#fc5e03")
  
  gg_train <- ggplot(data_metrics_train, aes(x=truth, y=pred))+
    geom_point()+
    labs(title = "Density of Data Point for Training Data",
         subtitle = paste("Linear Model: Rsq = ", round(metrics_train$.estimate[2], 3), sep = ""),
         x = paste("predicted values for ", i, sep = ""),
         y = paste("observed values for ", i, sep = "")) +
    theme(aspect.ratio=1)+
    geom_abline(intercept = 0, slope = 1, color = "#fc5e03")
  
  # arrange them into one file to export
  setwd(output)
  ggarrange(gg_test, gg_train, ncol = 2, nrow = 1,
            labels = c("A", "B"), common.legend = TRUE, legend = "right")
  ggsave(paste("./mechanistics/phylogenetic diversity/pasc_ID/training and predicting/linear model/", i, "_lm_prediction_pasc.png", sep = ""), width = 8, height = 5)
  setwd(input)
} # end i loop

#### Random Forest ----------------------------------------------------------------------------------------

for (i in metrics) {
  f <- as.formula(paste(i,  "~ intertegular_distance + proboscis_ratio + fore_wing_ratio + corbicula_ratio + landscape", sep = ""))
  train <- train(form = f,
                 data = ddf_train, 
                 method = "rf",
                 trControl = trainControl(method = "cv", number = 10))
  
  ## Get variable importance, and turn into a data frame
  var_imp <- varImp(train, scale=FALSE)$importance
  var_imp <- data.frame(variables=row.names(var_imp), importance=var_imp$Overall)
  
  ## Create a plot of variable importance
  var_imp %>%
    
    ## Sort the data by importance
    arrange(importance) %>%
    
    ## Create a ggplot object for aesthetic
    ggplot(aes(x=reorder(variables, importance), y=importance)) + 
    
    ## Plot the bar graph
    geom_bar(stat='identity') + 
    
    ## Flip the graph to make a horizontal bar plot
    coord_flip() + 
    
    ## Add x-axis label
    xlab('Variables') +
    
    ## Add a title
    labs(title='Random forest variable importance') + 
    
    ## Some layout for the plot
    theme_minimal() + 
    theme(axis.text = element_text(size = 10), 
          axis.title = element_text(size = 15), 
          plot.title = element_text(size = 20), 
    )
  setwd(output)
  ggsave(paste("./mechanistics/phylogenetic diversity/pasc_ID/training and predicting/random forest/", i, "_rf_var_imp_pasc.png", sep = ""), width = 8, height = 8)
  setwd(input)
  
  predict_train <- predict(
    ## lm object
    object=train, 
    ## Data to use for predictions
    newdata=ddf_train)
  
  predict_test <- predict(
    ## lm object
    object=train, 
    ## Data to use for predictions
    newdata=ddf_test)
  
  data_metrics_train <- data.frame(truth=ddf_train[, i], pred=predict_train)
  metrics_train <- metrics(data_metrics_train, truth, pred)
  
  data_metrics_test <- data.frame(truth=ddf_test[, i], pred=predict_test)
  metrics_test <- metrics(data_metrics_test, truth, pred)
  
  # plotting prediction
  library(patchwork)
  
  gg_test <- ggplot(data_metrics_test, aes(x=truth, y=pred))+
    geom_point()+
    labs(title = "Density of Data Point for Testing Data",
         subtitle = paste("Random Forest: Rsq = ", round(metrics_test$.estimate[2], 3), sep = ""),
         x = paste("predicted values for ", i, sep = ""),
         y = paste("observed values for ", i, sep = "")) +
    theme(aspect.ratio=1)+
    geom_abline(intercept = 0, slope = 1, color = "#fc5e03")
  
  gg_train <- ggplot(data_metrics_train, aes(x=truth, y=pred))+
    geom_point()+ 
    labs(title = "Density of Data Point for Training Data",
         subtitle = paste("Random Forest: Rsq = ", round(metrics_train$.estimate[2], 3), sep = ""),
         x = paste("predicted values for ", i, sep = ""),
         y = paste("observed values for ", i, sep = "")) +
    theme(aspect.ratio=1)+
    geom_abline(intercept = 0, slope = 1, color = "#fc5e03")
  
  # arrange them into one file to export
  setwd(output)
  ggarrange(gg_test, gg_train, ncol = 2, nrow = 1,
            labels = c("A", "B"), common.legend = TRUE, legend = "right")
  ggsave(paste("./mechanistics/phylogenetic diversity/pasc_ID/training and predicting/random forest/", i, "_rf_prediction_pasc.png", sep = ""), width = 8, height = 5)
  setwd(input)
  
} # end i loop

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

### modelling ----------------------------------------------------------------------------------------
# center and scale all the variables
BB22.ID[, c(7:15, 18:21)] <- scale(BB22.ID[, c(7:15, 17:20)],center=TRUE,scale=TRUE)
BB22.ID <- BB22.ID %>% filter(complete.cases(.)) # remove all entries with NAs

# correlation analysis
# look at the correlation between the explanatory variables
library(corrplot)
M <- cor(BB22.ID[, 7:15], use = "complete.obs") # subset data; only explanatory variables

# matrix of the p-value of the correlation
p.mat <- cor.mtest(BB22.ID[, 7:15])
head(p.mat)

#correlation plot
corrplot::corrplot(M, type="upper", order="hclust", 
                   p.mat = p.mat$p, sig.level = 0.01, tl.col = "black",
                   col = COL2('RdBu', 10)) # plot correlation with p-values

# !!! a lot are multicollinear. in the model only use: proboscis_ratio, fore_wing_ratio and corbicula_ratio


#  fit GLMM
#load the libraries
library(lme4)
library(car)
library(MuMIn)
library(arm)

#### Phylogenetic Species Variability  ----------------------------------------------------------------------------------------

# built an initial full model based on collinearity 
M1.full <- lmer(pvar ~ proboscis_ratio + fore_wing_ratio + corbicula_ratio + 
                  (1|landscape),
                data=BB22.ID)
summary(M1.full)

# check model assumptions
par(mar = c(8, 2, 5, 2))
fix.check(M1.full) # looks ok
vif <- vif(M1.full)
vif # looks good
mtext(side = 1, line = 6, adj = 0,
      paste("VIF:",
            paste(paste(names(vif)[1], round(vif[1], 3), sep=": "),
                  paste(names(vif)[2], round(vif[2], 3), sep=": "),
                  paste(names(vif)[3], round(vif[3], 3), sep=": "),
                  sep = ", "), sep = " ")
)

# formal test for spatial correlation
sims <- simulateResiduals(M1.full)
BB22.ID$site <- as.factor(BB22.ID$site)
simulationOutput <- recalculateResiduals(sims, group = BB22.ID$site)
testSpatialAutocorrelation(simulationOutput, x = BB22.sites.meta$LV95_x, y = BB22.sites.meta$LV95_y, plot = FALSE)
# there is no spatial autocorrelation

# update model: add intertegular distance to the model
M1.full.1 <- lmer(pvar ~ intertegular_distance + proboscis_ratio + fore_wing_ratio + 
                    corbicula_ratio + (1|landscape),
                  data=BB22.ID)
summary(M1.full.1)

# check model assumptions
par(mar = c(8, 2, 5, 2))
fix.check(M1.full.1) # looks ok
vif <- vif(M1.full.1)
vif # looks good
mtext(side = 1, line = 6, adj = 0,
      paste("VIF:",
            paste(paste(names(vif)[1], round(vif[1], 3), sep=": "),
                  paste(names(vif)[2], round(vif[2], 3), sep=": "),
                  paste(names(vif)[3], round(vif[3], 3), sep=": "),
                  paste(names(vif)[4], round(vif[4], 3), sep=": "),
                  sep = ", "), sep = " ")
)

# check which models fits better
anova(M1.full.1, M1.full) # without intertegular_distance better fit than M1.full

# dredging
# Things to take into account when dredging:
# 1) depends on the models included in the candidate set. You can’t identify a model as being the 
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

# no variable is having a significant effect in any model

# save the last eight produced plots in a selected path
plots.dir.path <- list.files(tempdir(), pattern="rs-graphics", full.names = TRUE);
plots.png.paths <- list.files(plots.dir.path, pattern=".png", full.names = TRUE)
file.copy(from = plots.png.paths,
          to= paste(output, "./mechanistics/phylogenetic diversity/lapi_ID/models", sep = ""))


## Model testing ----------------------------------------------------------------------------------------

### Data splitting  ----------------------------------------------------------------------------------------
library(rsample)
set.seed(123) # for reproducibility
BB22.ID.complete <- BB22.ID %>%
  filter(complete.cases(.))
split <- initial_split(BB22.ID.complete, prop = 0.8)
ddf_train <- training(split)
ddf_test <- testing(split) 

### Model training and predicting ----------------------------------------------------------------------------------------
library(caret) 
library(yardstick)

#### Linear Model ----------------------------------------------------------------------------------------
metrics <- "pvar"
for (i in metrics) {
  f <- as.formula(paste(i,  "~ intertegular_distance + proboscis_ratio + fore_wing_ratio + corbicula_ratio + landscape", sep = ""))
  train <- train(form = f,
                 data = ddf_train, 
                 method = "lm",
                 trControl = trainControl(method = "cv", number = 10))
  
  ## Get variable importance, and turn into a data frame
  var_imp <- varImp(train, scale=FALSE)$importance
  var_imp <- data.frame(variables=row.names(var_imp), importance=var_imp$Overall)
  
  ## Create a plot of variable importance
  var_imp %>%
    
    ## Sort the data by importance
    arrange(importance) %>%
    
    ## Create a ggplot object for aesthetic
    ggplot(aes(x=reorder(variables, importance), y=importance)) + 
    
    ## Plot the bar graph
    geom_bar(stat='identity') + 
    
    ## Flip the graph to make a horizontal bar plot
    coord_flip() + 
    
    ## Add x-axis label
    xlab('Variables') +
    
    ## Add a title
    labs(title='Linear model variable importance') + 
    
    ## Some layout for the plot
    theme_minimal() + 
    theme(axis.text = element_text(size = 10), 
          axis.title = element_text(size = 15), 
          plot.title = element_text(size = 20), 
    )
  setwd(output)
  ggsave(paste("./mechanistics/phylogenetic diversity/lapi_ID/training and predicting/linear model/", i, "_lm_var_imp_pasc.png", sep = ""), width = 8, height = 5)
  setwd(input)
  
  predict_train <- predict(
    ## lm object
    object=train, 
    ## Data to use for predictions
    newdata=ddf_train)
  
  predict_test <- predict(
    ## lm object
    object=train, 
    ## Data to use for predictions
    newdata=ddf_test)
  
  data_metrics_train <- data.frame(truth=ddf_train[,i], pred=predict_train) #mistake
  metrics_train <- metrics(data_metrics_train, truth, pred)
  
  data_metrics_test <- data.frame(truth=ddf_test[,i], pred=predict_test)
  metrics_test <- metrics(data_metrics_test, truth, pred)
  
  # plotting prediction
  library(patchwork)
  
  gg_test <- ggplot(data_metrics_test, aes(x=truth, y=pred))+
    geom_point()+
    labs(title = "Density of Data Point for Testing Data",
         subtitle = paste("Linear Model: Rsq = ", round(metrics_test$.estimate[2], 3), sep = ""), 
         x = paste("predicted values for ", i, sep = ""),
         y = paste("observed values for ", i, sep = "")) +
    theme(aspect.ratio=1)+
    geom_abline(intercept = 0, slope = 1, color = "#fc5e03")
  
  gg_train <- ggplot(data_metrics_train, aes(x=truth, y=pred))+
    geom_point()+
    labs(title = "Density of Data Point for Training Data",
         subtitle = paste("Linear Model: Rsq = ", round(metrics_train$.estimate[2], 3), sep = ""),
         x = paste("predicted values for ", i, sep = ""),
         y = paste("observed values for ", i, sep = "")) +
    theme(aspect.ratio=1)+
    geom_abline(intercept = 0, slope = 1, color = "#fc5e03")
  
  # arrange them into one file to export
  setwd(output)
  ggarrange(gg_test, gg_train, ncol = 2, nrow = 1,
            labels = c("A", "B"), common.legend = TRUE, legend = "right")
  ggsave(paste("./mechanistics/phylogenetic diversity/lapi_ID/training and predicting/linear model/", i, "_lm_prediction_pasc.png", sep = ""), width = 8, height = 5)
  setwd(input)
} # end i loop

#### Random Forest ----------------------------------------------------------------------------------------

for (i in metrics) {
  f <- as.formula(paste(i,  "~ intertegular_distance + proboscis_ratio + fore_wing_ratio + corbicula_ratio + landscape", sep = ""))
  train <- train(form = f,
                 data = ddf_train, 
                 method = "rf",
                 trControl = trainControl(method = "cv", number = 10))
  
  ## Get variable importance, and turn into a data frame
  var_imp <- varImp(train, scale=FALSE)$importance
  var_imp <- data.frame(variables=row.names(var_imp), importance=var_imp$Overall)
  
  ## Create a plot of variable importance
  var_imp %>%
    
    ## Sort the data by importance
    arrange(importance) %>%
    
    ## Create a ggplot object for aesthetic
    ggplot(aes(x=reorder(variables, importance), y=importance)) + 
    
    ## Plot the bar graph
    geom_bar(stat='identity') + 
    
    ## Flip the graph to make a horizontal bar plot
    coord_flip() + 
    
    ## Add x-axis label
    xlab('Variables') +
    
    ## Add a title
    labs(title='Random forest variable importance') + 
    
    ## Some layout for the plot
    theme_minimal() + 
    theme(axis.text = element_text(size = 10), 
          axis.title = element_text(size = 15), 
          plot.title = element_text(size = 20), 
    )
  setwd(output)
  ggsave(paste("./mechanistics/phylogenetic diversity/lapi_ID/training and predicting/random forest/", i, "_rf_var_imp_pasc.png", sep = ""), width = 8, height = 8)
  setwd(input)
  
  predict_train <- predict(
    ## lm object
    object=train, 
    ## Data to use for predictions
    newdata=ddf_train)
  
  predict_test <- predict(
    ## lm object
    object=train, 
    ## Data to use for predictions
    newdata=ddf_test)
  
  data_metrics_train <- data.frame(truth=ddf_train[, i], pred=predict_train)
  metrics_train <- metrics(data_metrics_train, truth, pred)
  
  data_metrics_test <- data.frame(truth=ddf_test[, i], pred=predict_test)
  metrics_test <- metrics(data_metrics_test, truth, pred)
  
  # plotting prediction
  library(patchwork)
  
  gg_test <- ggplot(data_metrics_test, aes(x=truth, y=pred))+
    geom_point()+
    labs(title = "Density of Data Point for Testing Data",
         subtitle = paste("Random Forest: Rsq = ", round(metrics_test$.estimate[2], 3), sep = ""),
         x = paste("predicted values for ", i, sep = ""),
         y = paste("observed values for ", i, sep = "")) +
    theme(aspect.ratio=1)+
    geom_abline(intercept = 0, slope = 1, color = "#fc5e03")
  
  gg_train <- ggplot(data_metrics_train, aes(x=truth, y=pred))+
    geom_point()+ 
    labs(title = "Density of Data Point for Training Data",
         subtitle = paste("Random Forest: Rsq = ", round(metrics_train$.estimate[2], 3), sep = ""),
         x = paste("predicted values for ", i, sep = ""),
         y = paste("observed values for ", i, sep = "")) +
    theme(aspect.ratio=1)+
    geom_abline(intercept = 0, slope = 1, color = "#fc5e03")
  
  # arrange them into one file to export
  setwd(output)
  ggarrange(gg_test, gg_train, ncol = 2, nrow = 1,
            labels = c("A", "B"), common.legend = TRUE, legend = "right")
  ggsave(paste("./mechanistics/phylogenetic diversity/lapi_ID/training and predicting/random forest/", i, "_rf_prediction_pasc.png", sep = ""), width = 8, height = 5)
  setwd(input)
  
} # end i loop
