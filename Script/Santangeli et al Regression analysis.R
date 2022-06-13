## ------------------------------------------------------------------------
## 'The effects of protected areas on the ecological niches of birds and mammals'
## ------------------------------------------------------------------------

# Santangeli et al.

## ------------------------------------------------------------------------
# 'R script to reproduce the analyses'
## ------------------------------------------------------------------------

# Analysis performed with R (v. R 4.1.0) and R studio (v. 1.4.1103)
# Authors: Stefano Mammola & Andrea Santangeli
# Location: Helsinki, April-July 2021

# clean the workspace -----------------------------------------------------

rm(list=ls())

# Loading R package -------------------------------------------------------

library("ape")        # Analyses of Phylogenetics and Evolution
library("caper")      # Comparative Analyses of Phylogenetics and Evolution in R
library("dplyr")      # A Grammar of Data Manipulation
library("geiger")     # Analysis of Evolutionary Diversification
library("grid")
library("ggstance")   # Horizontal 'ggplot2' Components
library("ggtree")     # an R package for visualization of tree and annotation data
library("gridExtra")  # Miscellaneous Functions for "Grid" Graphics
library("ggplot2")    # Create Elegant Data Visualisations Using the Grammar of Graphics
library("glmmTMB")    # Generalized Linear Mixed Models using Template Model Builder
library("phytools")   # Phylogenetic Tools for Comparative Biology (and Other Things)
library("png")        # For silouhettes
library("parameters") # Processing of Model Parameters
library("psych")      # Procedures for Psychological, Psychometric, and Personality research
library("patchwork")  # The Composer of Plots

# Custom theme for ggplot2 ------------------------------------------------

theme_ggplot <- theme(
  legend.position = "none",
  axis.title = element_text(size = 12),
  axis.text.x = element_text(size = 11),
  axis.text.y = element_text(size = 11),
  panel.grid = element_blank(),
  plot.caption = element_text(size = 10, color = "gray50"),
  plot.title = element_text(face="bold", size=12)
)

custom.label <- c("Body mass", "Habitat\nspecialization", "Diet\nspecialization", "Carnivore", "Red list status\n[threatened]")

# Loading silhouettes for figures ------------------------------------------

# Taken from PhyloPic (http://phylopic.org/) with open license.

Silu_bird   <- png::readPNG("Silhouettes/Silu_turdus.png") 
Silu_mamm   <- png::readPNG("Silhouettes/Silu_vulpes.png")

# Loading the databases ---------------------------------------------------

# Read the hypervolume results that includes the traits:
db <- read.csv("Data/hypervolume_metrics.csv", header = TRUE, sep = ",", dec = ".", row.names = NULL, as.is = FALSE)

# Read Trait matrix
trait <- read.csv("Data/traits.csv", header = TRUE, sep = ",", dec = ".", row.names = NULL, as.is = FALSE)

db <- dplyr::left_join(db, trait, by = c("species" = "Species")) 

#rename columns
colnames(db)[1]  <- "Species"
colnames(db)[18] <- "activity"
colnames(db)[19] <- "mass"
colnames(db)[22] <- "red"
colnames(db)[23] <- "hab"
colnames(db)[24] <- "diet"
colnames(db)[27] <- "lat"
colnames(db)[28] <- "Name_phylo"

#checking the strcuture
str(db)

#converting numeric to factor
db$red <- as.factor(db$red)

#Renaming the variable Name_phylo
db$Name_phylo <- gsub(' ', '_', db$Name_phylo)

# Loading Phylogenetic trees ----------------------------------------------

mammTREE <- read.tree("Data/mammTree.tre")
birdTREE <- read.tree("Data/birdTree.tre")

# Analysis with birds -----------------------------------------------------

birds <- db[db$Group == "bird",]

#### Data exploration

#Checking values in relation to axes 
dev.off()
plot(birds$naxes,birds$beta_total)

birds <- birds[birds$naxes>4,] ; birds <- na.omit(birds)

## checking outliers
dev.off() ; par(mfrow = c(4,2), mar = c(rep(2,4)))
dotchart(birds$beta_repl, main= "Beta repl") #ok
dotchart(birds$beta_diff, main= "Beta diff") #outlier to be omitted
dotchart(birds$hab, main= "habitat") #ok
dotchart(birds$diet, main= "diet") #maybe an outlier. Trasform?
dotchart(birds$Diet_vertebrates, main= "Diet Vertebrate") #log trasform!
dotchart(birds$mass, main = "mass") #log transform!
dotchart(birds$Gen_Length, main = "Gen_length") #log trasform!
dotchart(birds$Diet_vertebrates, main= "Diet_vert") #log trasform!

# transforming variables
birds$mass_log   <-  log(birds$mass+1)
birds$gen_log    <-  log(birds$Gen_Length+1)
birds$diet_asin  <-  asin(birds$diet)

#double check
dev.off() ; par(mfrow = c(2,2),mar = c(rep(2,4)))
dotchart(birds$mass_log, main= "mass_log") #ok
dotchart(birds$gen_log , main= "Gen_log ") #ok
dotchart(birds$diet_asin , main= "diet_asin") #ok

# Check collinearity
psych::pairs.panels(birds[,c("hab","diet_asin","gen_log","mass_log","Diet_vertebrates")]) 
#mass and generation are collinear! We select only mass

pdf(file = "Figures/Figure Collinearity bird.pdf", width = 7, height = 5)

psych::pairs.panels(birds[,c("hab","diet_asin","gen_log","mass_log","Diet_vertebrates")]) 

dev.off()

#Checking association with categorical
boxplot(mass_log ~ red, data = birds)
boxplot(hab ~ red, data = birds)
boxplot(diet_asin ~ red, data = birds)
boxplot(lat ~ red, data = birds)
#Ok

# Check relationships
psych::pairs.panels(birds[,c("beta_repl","beta_diff","hab","diet_asin","gen_log","mass_log","Diet_vertebrates")]) 

# GLMM ----------------------------------------

##Setting baseline for year
birds <- within(birds, year <- relevel(year, ref = "2018s"))

# scale and center all continuous variables before analyses
birds$mass_log          <-  scale(birds$mass_log, center = TRUE, scale = TRUE)
birds$hab               <-  scale(birds$hab, center = TRUE, scale = TRUE)
birds$diet_asin         <-  scale(birds$diet_asin, center = TRUE, scale = TRUE)
birds$Diet_vertebrates  <-  scale(birds$Diet_vertebrates, center = TRUE, scale = TRUE)
 
#creating a new dataframe
birds_model <- data.frame(
                 species  = birds$Species,
                 name_phylo = birds$Name_phylo,
                 order = birds$Order,
                 family = birds$Family,
                 year    = birds$year,
                 beta_repl = birds$beta_repl,
                 beta_diff = birds$beta_diff,
                 red = birds$red,
                 mass = birds$mass_log,
                 mass_log = birds$mass_log,
                 diet = birds$diet_asin,
                 vol = birds$delta_vol,
                 hab = birds$hab,
                 Diet_v = birds$Diet_vertebrates,
                 naxes = birds$naxes)

birds_model <- na.omit(birds_model)

# Model for replacement
bird_M1 <- glmmTMB(beta_repl ~ mass + hab + diet + Diet_v + red + (1 | year) + (1 | species),
                   data = birds_model, beta_family(link = "logit"))

#Checking the model
parameters::model_parameters(bird_M1)

# Fixed Effects 

# Parameter   | Coefficient |   SE |         95% CI |     z |      p
# ------------------------------------------------------------------
# (Intercept) |       -0.30 | 0.08 | [-0.46, -0.13] | -3.54 | < .001
# mass        |        0.14 | 0.08 | [-0.01,  0.30] |  1.83 | 0.067 
# hab         |        0.09 | 0.07 | [-0.05,  0.24] |  1.26 | 0.206 
# diet        |   -5.04e-03 | 0.08 | [-0.16,  0.15] | -0.07 | 0.948 
# Diet_v      |       -0.13 | 0.07 | [-0.27,  0.02] | -1.73 | 0.084 
# red [1]     |        0.05 | 0.24 | [-0.42,  0.53] |  0.22 | 0.826 

# # Random Effects 
# 
# Parameter               | Coefficient
# -------------------------------------
# SD (Intercept: year)    |        0.01
# SD (Intercept: species) |        0.42
# SD (Residual)           |        1.78

performance::check_model(bird_M1)

pdf(file = "Figures/Validation M1.pdf", width = 8, height = 6)
performance::check_model(bird_M1)
dev.off()

#Checking for phylogenetic signal in the residuals
res <- residuals(bird_M1)
birds_test <- data.frame(birds_model, 
                      residuals = residuals(bird_M1)) 

birds_test <- birds_test %>% group_by(name_phylo) %>% summarise_at(vars(c("residuals")), funs(mean(., na.rm=TRUE)))
birds_test <- data.frame(birds_test)

#checking match of names in the phylogenetic tree (and dropping unused tips)
rownames(birds_test) <- make.names(birds_test$name_phylo, unique = TRUE)

(missing.tip <- birdTREE$tip.label[birdTREE$tip.label %in% birds_test$name_phylo == FALSE]) # 7 species to drop

length(birds_test$name_phylo)
sort(birdTREE$tip.label)

birdTREE2 <- drop.tip(birdTREE, missing.tip)

unique(birdTREE2$tip.label %in% birds_test$name_phylo) # all good now

#Checking signal in the residuals
BirdTreeList <- birdTREE2$tip.label      # get order of the sp in tree
birds_test   <- birds_test[order(match(rownames(birds_test), BirdTreeList)),] # reorder the DF

phylosig(birdTREE2, birds_test$residuals, method = "lambda", test = TRUE)

# Phylogenetic signal lambda : 6.68102e-05 
# logL(lambda) : 62.674 
# LR(lambda=0) : -0.003054 
# P-value (based on LR test) : 1 

# No phylo signal in the residuals!

mod <- parameters::model_parameters(bird_M1) 

mod1 <- data.frame(Variable  = mod$Parameter[2:6],
                   Estimate  = mod$Coefficient[2:6],
                   CI_low    = mod$CI_low[2:6],
                   CI_high   = mod$CI_high[2:6]) ; rm(mod)

mod1$Variable <- custom.label
mod1$Variable <- factor(mod1$Variable, rev(custom.label)) 

(p3 <- ggplot2::ggplot(data = mod1, aes(Variable,Estimate)) +
    geom_hline(lty = 3, size = 0.7, col = "grey50", yintercept = 0) +
    geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = 0, col = "grey10") +
    geom_point(size = 2, pch = 21, col = "grey10", fill = "grey20") +
    geom_text(aes(Variable, Estimate), 
              label = round(mod1$Estimate,2), 
              vjust = -1, size = 3) +
    labs(title = "C - Habitat shift",
         y = "Estimate [95% Confidence interval]",
         x = NULL) + 
    ylim(-1.2,1.2)+
    coord_flip() +
    annotation_custom(grid::rasterGrob(Silu_bird),
                      xmin = unit(0.7, "native"), xmax = unit(1.2,"native"),
                      ymin = unit(0.7,"npc"),  ymax = unit(1,"npc"))+ 
    theme_bw() + 
    theme_ggplot)

# Model for expantion/contraction

#removing the outlier
birds_model2  <- birds_model[birds_model$beta_diff < 0.2,]
birds_model2  <- na.omit(birds_model2)

#Doesn't work!
bird_M2 <- glmmTMB(beta_diff ~ mass + hab + diet + Diet_v + red + (1 | year) + (1 | species), 
                   data = birds_model2, beta_family(link = "logit"))
#Error messagge: there are 5 value = 0 in the response variable.
range(db$beta_diff) 
table(db$beta_diff)

#replacing the value with a very small number
birds_model2$beta_diff <- ifelse(birds_model2$beta_diff == 0, 
                                 birds_model2$beta_diff + 0.000000000000000000001,
                                 birds_model2$beta_diff)

#Re-fitting the model
bird_M2 <- glmmTMB(beta_diff ~ mass + hab + diet + Diet_v + red + (1 | year) + (1 | species), 
                   data = birds_model2, beta_family(link = "logit"))

#Checking the model
parameters::model_parameters(bird_M2)

# Parameter   | Coefficient |   SE |         95% CI |      z |      p
# -------------------------------------------------------------------
# (Intercept) |       -3.98 | 0.14 | [-4.26, -3.71] | -28.17 | < .001
# mass        |       -0.04 | 0.10 | [-0.24,  0.15] |  -0.45 | 0.656 
# hab         |        0.09 | 0.10 | [-0.10,  0.27] |   0.89 | 0.375 
# diet        |       -0.02 | 0.10 | [-0.21,  0.18] |  -0.18 | 0.859 
# Diet_v      |        0.04 | 0.08 | [-0.13,  0.20] |   0.44 | 0.663 
# red [1]     |       -0.02 | 0.32 | [-0.65,  0.61] |  -0.07 | 0.945

# # Random Effects 
# 
# Parameter               | Coefficient
# -------------------------------------
# SD (Intercept: year)    |        0.15
# SD (Intercept: species) |        0.76
# SD (Residual)           |        5.68

performance::check_model(bird_M2)

pdf(file = "Figures/Validation M2.pdf", width = 8, height = 6)
performance::check_model(bird_M2)
dev.off()

#Checking for phylogenetic signal in the residuals
res <- residuals(bird_M2)
birds_test <- data.frame(birds_model2, 
                         residuals = residuals(bird_M2)) 

birds_test <- birds_test %>% group_by(name_phylo) %>% summarise_at(vars(c("residuals")), funs(mean(., na.rm=TRUE)))

birds_test <- data.frame(birds_test)

unique(birdTREE2$tip.label %in% birds_test$name_phylo) # all good now

#Checking the signal in the residuals
BirdTreeList <- birdTREE2$tip.label      # get order of the sp in tree
birds_test <- birds_test[order(match(rownames(birds_test), BirdTreeList)),] # reorder the DF

phylosig(birdTREE2, birds_test$residuals, method = "lambda", test = TRUE)

# Phylogenetic signal lambda : 6.68102e-05 
# logL(lambda) : 291.166 
# LR(lambda=0) : -0.00459575 
# P-value (based on LR test) : 1 

# No phylo signal in the residuals! 

mod <- parameters::model_parameters(bird_M2) 

mod2 <- data.frame(Variable  = mod$Parameter[2:6],
                   Estimate  = mod$Coefficient[2:6],
                   CI_low    = mod$CI_low[2:6],
                   CI_high   = mod$CI_high[2:6]) ; rm(mod)

mod2$Variable <- custom.label
mod2$Variable <- factor(mod1$Variable, rev(custom.label)) 

(p1 <- ggplot2::ggplot(data = mod2, aes(Variable,Estimate)) +
    geom_hline(lty = 3, size = 0.7, col = "grey50", yintercept = 0) +
    geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = 0, col = "grey10") +
    geom_point(size = 2, pch = 21, col = "grey10", fill = "grey20") +
    geom_text(aes(Variable, Estimate), 
              label = round(mod2$Estimate,2), 
              vjust = -1, size = 3) +
    labs(title = "A - Niche expansion",
         y = NULL,
         x = NULL) + 
    ylim(-1.2,1.2)+
    coord_flip() +
    annotation_custom(grid::rasterGrob(Silu_bird),
                      xmin = unit(0.7, "native"), xmax = unit(1.2,"native"),
                      ymin = unit(0.7,"npc"),  ymax = unit(1,"npc"))+ 
    theme_bw() + 
    theme_ggplot)

# Analysis with mammals -----------------------------------------------------

mamm <- db[db$Group == "mammal",] ; droplevels(mamm)

#### Data exploration

#Checking values in relation to axes 
dev.off()
plot(mamm$naxes,mamm$beta_total)

mamm <- mamm[mamm$naxes>4,]

## checking outliers
dev.off() ; par(mfrow=c(4,2),mar=c(rep(2,4)))
dotchart(mamm$beta_repl, main= "Beta repl") #outlier
dotchart(mamm$beta_diff, main= "Beta diff") #outlier 
dotchart(mamm$hab, main= "habitat") #maybe an outlier. Trasform?
dotchart(mamm$diet, main= "diet") #maybe an outlier. Trasform?
dotchart(mamm$mass, main= "mass") #log transform!
dotchart(mamm$Diet_vertebrates, main= "Carnivore") #ok
dotchart(mamm$Gen_Length, main= "Gen_length") #ok

# trasforming variables
mamm$mass_log   <-  log(mamm$mass+1)
mamm$diet_asin  <-  asin(mamm$diet)
mamm$hab_asin   <-  asin(mamm$hab)

#double check
dev.off() ; par(mfrow=c(2,2),mar=c(rep(2,4)))
dotchart(mamm$mass_log, main= "mass_log") #ok
dotchart(mamm$hab_asin , main= "hab_asin ") #ok
dotchart(mamm$diet_asin , main= "diet_asin") #ok

# Check collinearity
psych::pairs.panels(mamm[,c("hab_asin","diet_asin","mass_log","Gen_Length","Diet_vertebrates")]) 
#mass and generation are collinear! We select only mass

# saving the figure
pdf(file = "Figures/Figure Collinearity mammals.pdf", width = 7, height = 5)
psych::pairs.panels(mamm[,c("hab_asin","diet_asin","mass_log","Gen_Length","Diet_vertebrates")]) 
dev.off()

#Checking association with categorical
boxplot(mass_log ~ red, data = mamm)
boxplot(hab_asin ~ red, data = mamm)
boxplot(diet_asin ~ red, data = mamm)
boxplot(lat ~ red, data = mamm)
#Ok

# Check relationships
psych::pairs.panels(mamm[,c("beta_repl","beta_diff","hab_asin","diet_asin","mass_log","Gen_Length","Diet_vertebrates")]) 

# GLMM ----------------------------------------

##Setting baseline for year
mamm <- within(mamm, year <- relevel(year, ref = "2018s"))

# scale and center all continuous variables before analyses
mamm$mass_log  <- scale(mamm$mass_log, center = TRUE, scale = TRUE)
mamm$hab       <- scale(mamm$hab_asin, center = TRUE, scale = TRUE)
mamm$diet_asin <- scale(mamm$diet_asin, center = TRUE, scale = TRUE)
mamm$Diet_vertebrates <- scale(mamm$Diet_vertebrates, center = TRUE, scale = TRUE)
 
#creating a new dataframe
mamm_model <- data.frame(species = mamm$Species,
                         name_phylo = mamm$Name_phylo,
                         vol = mamm$delta_vol,
                         order = mamm$Order,
                         family = mamm$Family,
                         year    = mamm$year,
                         beta_repl = mamm$beta_repl,
                         beta_diff = mamm$beta_diff,
                         red = mamm$red,
                         mass = mamm$mass_log,
                         diet =mamm$diet_asin,
                         hab = mamm$hab,
                         Diet_vertebrates = mamm$Diet_vertebrates)

mamm_model <- na.omit(mamm_model)

mamm_model1 <- mamm_model$beta_repl

# Model for replacement
mamm_M1 <- glmmTMB(beta_repl ~ mass + hab + diet + Diet_vertebrates + red + (1|year) + (1|species), data = mamm_model, beta_family(link = "logit"))

#Checking the model
parameters::model_parameters(mamm_M1)

# # Fixed Effects 
# 
# Parameter        | Coefficient |   SE |         95% CI |     z |     p
# ----------------------------------------------------------------------
# (Intercept)      |       -0.33 | 0.11 | [-0.55, -0.10] | -2.87 | 0.004
# mass             |       -0.09 | 0.12 | [-0.32,  0.15] | -0.73 | 0.468
# hab              |       -0.11 | 0.10 | [-0.31,  0.10] | -1.00 | 0.316
# diet             |        0.04 | 0.10 | [-0.16,  0.25] |  0.40 | 0.687
# Diet_vertebrates |       -0.11 | 0.12 | [-0.35,  0.13] | -0.90 | 0.371
# red [1]          |       -0.04 | 0.29 | [-0.61,  0.52] | -0.14 | 0.886

# # Random Effects 
# 
# Parameter               | Coefficient
# -------------------------------------
# SD (Intercept: year)    |        0.04
# SD (Intercept: species) |        0.37
# SD (Residual)           |        4.39

performance::check_model(mamm_M1)

pdf(file = "Figures/Validation M3.pdf", width = 8, height = 6)
performance::check_model(mamm_M1)
dev.off()

#Checking for phylogenetic signal in the residuals
res <- residuals(mamm_M1)
mamm_test <- data.frame(mamm_model, 
                        residuals = residuals(mamm_M1)) 

mamm_test <- mamm_test %>% group_by(name_phylo) %>% summarise_at(vars(c("residuals")), funs(mean(., na.rm=TRUE)))

mamm_test <- data.frame(mamm_test)
rownames(mamm_test) <- make.names(mamm_test$name_phylo, unique = TRUE)

# Checking match of names in the phylogenetic tree 
unique(mammTREE$tip.label %in% mamm_test$name_phylo) # all good now

MammTreeList <- mammTREE$tip.label      # get order of the sp in tree
mamm_test    <- mamm_test[order(match(rownames(mamm_test), MammTreeList)),] # reorder the DF

phylosig(mammTREE, mamm_test$residuals, method = "lambda", test = TRUE)

# Phylogenetic signal lambda : 6.94123e-05 
# logL(lambda) : 34.2737 
# LR(lambda=0) : -0.000826767 
# P-value (based on LR test) : 1 

# No phylo signal in the residuals! 

mod <- parameters::model_parameters(mamm_M1) 

mod3 <- data.frame(Variable  = mod$Parameter[2:6],
                   Estimate  = mod$Coefficient[2:6],
                   CI_low    = mod$CI_low[2:6],
                   CI_high   = mod$CI_high[2:6]) ; rm(mod)

mod3$Variable <- custom.label
mod3$Variable <- factor(mod3$Variable, rev(custom.label)) 

(p4 <- ggplot2::ggplot(data = mod3, aes(Variable,Estimate)) +
    geom_hline(lty = 3, size = 0.7, col = "grey50", yintercept = 0) +
    geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = 0, col = "grey10") +
    geom_point(size = 2, pch = 21, col = "grey10", fill = "grey20") +
    geom_text(aes(Variable, Estimate), 
              label = round(mod3$Estimate,2), 
              vjust = -1, size = 3) +
    labs(title = "D - Habitat shift",
         y = "Estimate [95% Confidence interval]",
         x = NULL) + 
    ylim(-1.2,1.2)+
    coord_flip() +
    annotation_custom(grid::rasterGrob(Silu_mamm),
                      xmin = unit(0.7, "native"), xmax = unit(1.2,"native"),
                      ymin = unit(0.7,"npc"),  ymax = unit(1,"npc"))+ 
    theme_bw() + 
    theme_ggplot + theme(axis.text.y=element_blank()))

# Adding jitter for modelling
mamm_model$beta_diff <- ifelse(mamm_model$beta_diff == 0, 
                               mamm_model$beta_diff + 0.000000000000000000001,
                               mamm_model$beta_diff)

mamm_M2 <- glmmTMB(beta_diff ~ mass + hab + diet + Diet_vertebrates + red +(1|year) + (1|species), data = mamm_model, beta_family(link = "logit"))

#Checking the model
parameters::model_parameters(mamm_M2)

# # Fixed Effects 
# 
# Parameter        | Coefficient |   SE |         95% CI |      z |      p
# ------------------------------------------------------------------------
# (Intercept)      |       -3.12 | 0.10 | [-3.32, -2.92] | -30.68 | < .001
# mass             |       -0.13 | 0.09 | [-0.32,  0.06] |  -1.36 | 0.173 
# hab              |       -0.02 | 0.09 | [-0.20,  0.16] |  -0.23 | 0.816 
# diet             |        0.12 | 0.09 | [-0.05,  0.30] |   1.40 | 0.162 
# Diet_vertebrates |        0.02 | 0.10 | [-0.17,  0.21] |   0.23 | 0.818 
# red [1]          |       -0.48 | 0.27 | [-1.02,  0.06] |  -1.74 | 0.081 

# # 
# # Random Effects 
# 
# Parameter               | Coefficient
# -------------------------------------
# SD (Intercept: year)    |    1.40e-09
# SD (Intercept: species) |    1.06e-04
# SD (Residual)           |        6.23

pdf(file = "Figures/Validation M4.pdf", width = 8, height = 6)
performance::check_model(mamm_M2)
dev.off()

#Checking for phylogenetic signal in the residuals
res <- residuals(mamm_M2)
mamm_test <- data.frame(mamm_model, 
                        residuals = residuals(mamm_M2)) 

mamm_test <- mamm_test %>% group_by(name_phylo) %>% summarise_at(vars(c("residuals")), funs(mean(., na.rm=TRUE)))

mamm_test <- data.frame(mamm_test)
rownames(mamm_test) <- make.names(mamm_test$name_phylo, unique = TRUE)

MammTreeList <- mammTREE$tip.label      # get order of the sp in tree
mamm_test <- mamm_test[order(match(rownames(mamm_test), MammTreeList)),] # reorder the DF

phylosig(mammTREE, mamm_test$residuals, method = "lambda", test = TRUE)

# Phylogenetic signal lambda : 6.94123e-05 
# logL(lambda) : 34.2737 
# LR(lambda=0) : -0.000826767 
# P-value (based on LR test) : 1 

# No phylo signal in the residuals!

mod <- parameters::model_parameters(mamm_M2) 

mod4 <- data.frame(Variable  = mod$Parameter[2:6],
                   Estimate  = mod$Coefficient[2:6],
                   CI_low    = mod$CI_low[2:6],
                   CI_high   = mod$CI_high[2:6]) ; rm(mod)

mod4$Variable <- custom.label
mod4$Variable <- factor(mod4$Variable, rev(custom.label)) 

(p2 <- ggplot2::ggplot(data = mod4, aes(Variable,Estimate)) +
    geom_hline(lty = 3, size = 0.7, col = "grey50", yintercept = 0) +
    geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = 0, col = "grey10") +
    geom_point(size = 2, pch = 21, col = "grey10", fill = "grey20") +
    geom_text(aes(Variable, Estimate), 
              label = round(mod4$Estimate,2), 
              vjust = -1, size = 3) +
    labs(title = "B - Niche expansion",
         y = NULL,
         x = NULL) + 
    ylim(-1.2,1.2)+
    coord_flip() +
    annotation_custom(grid::rasterGrob(Silu_mamm),
                      xmin = unit(0.7, "native"), xmax = unit(1.2,"native"),
                      ymin = unit(0.7,"npc"),  ymax = unit(1,"npc"))+ 
    theme_bw() + 
    theme_ggplot + theme(axis.text.y=element_blank()))

# Arranging in a plot ------------------------------------------------------

plots <- list(p1, p2, p3, p4) ; grobs <- list() ; widths <- list()

for (i in 1:length(plots)){
  grobs[[i]] <- ggplotGrob(plots[[i]])
  widths[[i]] <- grobs[[i]]$widths[2:5]
}

maxwidth <- do.call(grid::unit.pmax, widths)

for (i in 1:length(grobs)){
  grobs[[i]]$widths[2:5] <- as.list(maxwidth)
}

pdf(file = "Figures/Figure 5.pdf", width = 12, height = 10)
do.call("grid.arrange", c(grobs, nrow = 2, ncol = 2))
dev.off()

# Plotting the value on the phylogeny -------------------------------------

# Code partly rearranged from: https://4va.github.io/biodatasci/r-ggtree.html
mamm_tree_plot <- mamm %>% group_by(Name_phylo) %>% summarise_at(vars(c("beta_total","beta_repl","beta_diff")), funs(mean(., na.rm=TRUE)))

mamm_tree_plot <- data.frame(mamm_tree_plot)
rownames(mamm_tree_plot) <- mamm_tree_plot$Name_phylo

MammTreeList <- mammTREE$tip.label      # get order of the sp in tree
mamm_tree_plot <- mamm_tree_plot [order(match(rownames(mamm_tree_plot), MammTreeList)),] # reorder the DF

mammTREE$tip.label <- paste("  ", mammTREE$tip.label,sep='')
mammTREE$tip.label <- gsub('_', ' ', mammTREE$tip.label)

#renaming Castor fiber to castor sp.
mammTREE$tip.label[19] <- "  Castor sp."

length(mamm_tree_plot$beta_repl)
length(mammTREE$tip.label)

d1 <- data.frame(id=rep(mammTREE$tip.label,2),
                 value = c(mamm_tree_plot$beta_repl,mamm_tree_plot$beta_diff),
                 category = c(rep("Shift",length(mammTREE$tip.label)),
                              rep("Expansion",length(mammTREE$tip.label))))

p <-  ggtree(mammTREE, branch.length="none") + 
      geom_tippoint() + 
      geom_tiplab(fontface="italic") + 
      xlim_tree(26) 

p2 <- facet_plot(p, panel = 'Niche differentiation', data = d1, 
                 geom = geom_barh,
                 mapping = aes(x = value, fill = as.factor(category)), 
                 stat='identity' ) + scale_fill_manual(values =  c("grey60", "grey10")) +
                 labs(fill = NULL)
                 
p2 + theme_tree2()

#Save the plot
pdf(file = "Figures/Figure 4.pdf", width = 8, height = 4)
p2 + theme_tree2()
dev.off()

### Birds

birds_tree_plot <- birds %>% group_by(Name_phylo) %>% summarise_at(vars(c("beta_total","beta_repl","beta_diff")), funs(mean(., na.rm=TRUE)))

birds_tree_plot <- data.frame(birds_tree_plot)
rownames(birds_tree_plot) <- birds_tree_plot$Name_phylo

BirdTreeList <- birdTREE2$tip.label      # get order of the sp in tree
birds_tree_plot <- birds_tree_plot[order(match(rownames(birds_tree_plot), BirdTreeList)),] # reorder the DF

birdTREE2$tip.label <- paste("  ", birdTREE2$tip.label,sep='')
birdTREE2$tip.label <- gsub('_', ' ', birdTREE2$tip.label)

d2 <- data.frame(id=rep(birdTREE2$tip.label,2),
                 value = c(birds_tree_plot$beta_repl,birds_tree_plot$beta_diff),
                 category = c(rep("Shift",length(birdTREE2$tip.label)),
                              rep("Expansion",length(birdTREE2$tip.label))))

p3 <-  ggtree(birdTREE2, branch.length="none") + 
  geom_tippoint() + 
  geom_tiplab(size=2,fontface="italic") + 
  xlim_tree(26) 

p4 <- facet_plot(p3, panel = 'Niche differentiation', data = d2, 
                 geom = geom_barh, 
                 mapping = aes(x = value, fill = as.factor(category)), 
                 stat='identity' ) + scale_fill_manual(values =  c("grey60", "grey10")) +
  labs(fill = NULL)

p4 + theme_tree2()

pdf(file = "Figures/Figure 3.pdf", width = 10, height = 10)
p4 + theme_tree2()
dev.off()

# How many species for each hypothesis ---------------------------------

#Threshold of 0.1
threshold <- 0.1

sum(nrow(mamm_tree_plot[mamm_tree_plot$beta_repl < threshold & mamm_tree_plot$beta_diff < threshold,])) # HP 0
sum(nrow(mamm_tree_plot[mamm_tree_plot$beta_diff >= threshold,])) # HP 1
sum(nrow(mamm_tree_plot[mamm_tree_plot$beta_repl >= threshold,])) # HP 2
sum(nrow(mamm_tree_plot[mamm_tree_plot$beta_repl >= threshold & mamm_tree_plot$beta_diff >= threshold,])) # HP 3

sum(nrow(birds_tree_plot[birds_tree_plot$beta_repl < threshold & birds_tree_plot$beta_diff < threshold,])) # HP 0
sum(nrow(birds_tree_plot[birds_tree_plot$beta_diff >= threshold,])) # HP 1
sum(nrow(birds_tree_plot[birds_tree_plot$beta_repl >= threshold,])) # HP 2
sum(nrow(birds_tree_plot[birds_tree_plot$beta_repl >= threshold & birds_tree_plot$beta_diff >= threshold,])) # HP 3

HP0_mamm <- mamm_tree_plot[mamm_tree_plot$beta_repl < threshold & mamm_tree_plot$beta_diff < threshold,]$Name_phylo
HP2_mamm <- mamm_tree_plot[mamm_tree_plot$beta_repl >= threshold,]$Name_phylo

HP0_birds <- birds_tree_plot[birds_tree_plot$beta_repl < threshold & birds_tree_plot$beta_diff < threshold,]$Name_phylo
HP2_birds <- birds_tree_plot[birds_tree_plot$beta_repl >= threshold,]$Name_phylo

# Bind all 4 df and reorganize them
All_HP <- data.frame(species = c(HP0_birds, HP0_mamm, HP2_birds, HP2_mamm),
                     HP = c(rep("HP0_birds",length(HP0_birds)),
                            rep("HP0_mamm",length(HP0_mamm)),
                            rep("HP2_birds",length(HP2_birds)),
                            rep("HP2_mamm",length(HP2_mamm))))

trait <- data.frame(Name_phylo = db$Name_phylo, Red = db$red)

trait <- distinct(trait, Name_phylo, .keep_all = TRUE)

# now I have all data with red list status:
All_HP2 <- left_join(All_HP, trait, 
                     by = c("species" = "Name_phylo")) ; names(All_HP2)

All_HP2$HP  <- as.factor(All_HP2$HP)
All_HP2$Red <- as.factor(All_HP2$Red)

# Summarizing data

Data4hp <- All_HP2 %>% group_by(HP, Red) %>% tally()
Data4hp <- as.data.frame(Data4hp)

### Further manipulation
HP <- rep(c("HP0_birds", "HP1_birds","HP2_birds","HP3_birds", "HP0_mamm", "HP1_mamm","HP2_mamm","HP3_mamm"),2)
Red <- c(rep(0,8), rep(1,8))
HP <- data.frame(HP, Red)
HP$HP <- as.factor(HP$HP)
HP$Red <- as.factor(HP$Red)
HP2 <- left_join(HP, Data4hp, by = c("HP", "Red"))

# Turn Nas in the n to Zeros:
HP2 <- HP2 %>% mutate_at(vars("n"), ~replace(., is.na(.), 0)) 

# Now I have the DF for the figures:
birdsDF <- HP2[c(1:4, 9:12) ,]  ; birdsDF <- gdata::drop.levels(birdsDF)
mammDF  <- HP2[c(5:8, 13:16) ,] ; mammDF  <- gdata::drop.levels(mammDF) 

#Rename levels
levels(birdsDF$HP) <- levels(mammDF$HP) <- c("HP0\nNo change", "HP1\nExpansion", "HP2\nShift", "HP3\nExpansion & shift")  # relevel
birdsDF$Red <- relevel(birdsDF$Red, "1")    # order levels, threatened first
mammDF$Red  <- relevel(mammDF$Red, "1")     # order levels, threatened first

#creating a fake low value for graphical reasons
birdsDF[c(2,4),3] <- 0.5
mammDF[c(2,4),3] <- 0.1

# Make plots:
(p7 <- ggplot(birdsDF, aes(x = HP, y = n)) +
    geom_col(aes(fill = Red), width = 0.7) +
    labs(x = NULL, 
         y = "Count of species",
         title = "A")+
    scale_fill_manual(values = c("grey60", "grey10")) +
    annotation_custom(grid::rasterGrob(Silu_bird),
                      xmin = 3.5, xmax = 4.5,
                      ymin = 65,  ymax = 75)
  + theme_bw() + theme_ggplot)

(p8 <-  ggplot(mammDF, aes(x = HP, y = n)) +
    geom_col(aes(fill = Red), width = 0.7) +
    labs(x = NULL, 
         y = NULL,
         title = "B")+
    
    annotation_custom(grid::rasterGrob(Silu_mamm),
                      xmin = 3.5, xmax = 4.5,
                      ymin = 15,  ymax = 17)+
    
    scale_fill_manual(values = c("grey60", "grey10"), labels = c("Threatened", "Not threatened"), name =  NULL)
  + theme_bw() + theme_ggplot + theme(legend.position = c(0.2, 0.7))
)

# Arranging in a plot ------------------------------------------------------

plots <- list(p7, p8)
grobs <- list()
widths <- list()

for (i in 1:length(plots)){
  grobs[[i]] <- ggplotGrob(plots[[i]])
  widths[[i]] <- grobs[[i]]$widths[2:5]
}

maxwidth <- do.call(grid::unit.pmax, widths)

for (i in 1:length(grobs)){
  grobs[[i]]$widths[2:5] <- as.list(maxwidth)
}

pdf(file = "Figures/Figure S6_thresold 1.pdf", width = 10, height = 3)

do.call("grid.arrange", c(grobs, nrow = 1, ncol = 2))

dev.off()

# Theshold 0.2
threshold <- 0.2

sum(nrow(mamm_tree_plot[mamm_tree_plot$beta_repl < threshold & mamm_tree_plot$beta_diff < threshold,])) # HP 0
sum(nrow(mamm_tree_plot[mamm_tree_plot$beta_diff >= threshold,])) # HP 1
sum(nrow(mamm_tree_plot[mamm_tree_plot$beta_repl >= threshold,])) # HP 2
sum(nrow(mamm_tree_plot[mamm_tree_plot$beta_repl >= threshold & mamm_tree_plot$beta_diff >= threshold,])) # HP 3

sum(nrow(birds_tree_plot[birds_tree_plot$beta_repl < threshold & birds_tree_plot$beta_diff < threshold,])) # HP 0
sum(nrow(birds_tree_plot[birds_tree_plot$beta_diff >= threshold,])) # HP 1
sum(nrow(birds_tree_plot[birds_tree_plot$beta_repl >= threshold,])) # HP 2
sum(nrow(birds_tree_plot[birds_tree_plot$beta_repl >= threshold & birds_tree_plot$beta_diff >= threshold,])) # HP 3

HP0_mamm <- mamm_tree_plot[mamm_tree_plot$beta_repl < threshold & mamm_tree_plot$beta_diff < threshold,]$Name_phylo
HP2_mamm <- mamm_tree_plot[mamm_tree_plot$beta_repl >= threshold,]$Name_phylo

HP0_birds <- birds_tree_plot[birds_tree_plot$beta_repl < threshold & birds_tree_plot$beta_diff < threshold,]$Name_phylo
HP2_birds <- birds_tree_plot[birds_tree_plot$beta_repl >= threshold,]$Name_phylo

# Bind all 4 df and reorganize them
All_HP <- data.frame(species = c(HP0_birds, HP0_mamm, HP2_birds, HP2_mamm),
                     HP = c(rep("HP0_birds",length(HP0_birds)),
                            rep("HP0_mamm",length(HP0_mamm)),
                            rep("HP2_birds",length(HP2_birds)),
                            rep("HP2_mamm",length(HP2_mamm))))

trait <- data.frame(Name_phylo = db$Name_phylo, Red = db$red)

trait <- distinct(trait, Name_phylo, .keep_all = TRUE)

# now I have all data with red list status:
All_HP2 <- left_join(All_HP, trait, 
                     by = c("species" = "Name_phylo")) ; names(All_HP2)

All_HP2$HP  <- as.factor(All_HP2$HP)
All_HP2$Red <- as.factor(All_HP2$Red)

# Summarizing data

Data4hp <- All_HP2 %>% group_by(HP, Red) %>% tally()
Data4hp <- as.data.frame(Data4hp)

### Further manipulation
HP <- rep(c("HP0_birds", "HP1_birds","HP2_birds","HP3_birds", "HP0_mamm", "HP1_mamm","HP2_mamm","HP3_mamm"),2)
Red <- c(rep(0,8), rep(1,8))
HP <- data.frame(HP, Red)
HP$HP <- as.factor(HP$HP)
HP$Red <- as.factor(HP$Red)
HP2 <- left_join(HP, Data4hp, by = c("HP", "Red"))

# Turn Nas in the n to Zeros:
HP2 <- HP2 %>% mutate_at(vars("n"), ~replace(., is.na(.), 0)) 

#           HP Red  n
# 1  HP0_birds   0 13
# 2  HP1_birds   0  0
# 3  HP2_birds   0 69
# 4  HP3_birds   0  0
# 5   HP0_mamm   0  1
# 6   HP1_mamm   0  0
# 7   HP2_mamm   0 16
# 8   HP3_mamm   0  0
# 9  HP0_birds   1  2
# 10 HP1_birds   1  0
# 11 HP2_birds   1  8
# 12 HP3_birds   1  0
# 13  HP0_mamm   1  1
# 14  HP1_mamm   1  0
# 15  HP2_mamm   1  3
# 16  HP3_mamm   1  0

# Now I have the DF for the figures:
birdsDF <- HP2[c(1:4, 9:12) ,]  ; birdsDF <- gdata::drop.levels(birdsDF)
mammDF  <- HP2[c(5:8, 13:16) ,] ; mammDF  <- gdata::drop.levels(mammDF) 

#Rename levels
levels(birdsDF$HP) <- levels(mammDF$HP) <- c("HP0\nNo change", "HP1\nExpansion", "HP2\nShift", "HP3\nExpansion & shift")  # relevel
birdsDF$Red <- relevel(birdsDF$Red, "1")    # order levels, threatened first
mammDF$Red  <- relevel(mammDF$Red, "1")     # order levels, threatened first

#creating a fake low value for graphical reasons
birdsDF[c(2,4),3] <- 0.5
mammDF[c(2,4),3] <- 0.1

# Make plots:
(p7 <- ggplot(birdsDF, aes(x = HP, y = n)) +
    geom_col(aes(fill = Red), width = 0.7) +
    labs(x = NULL, 
         y = "Count of species",
         title = "A")+
    scale_fill_manual(values = c("grey60", "grey10")) +
    annotation_custom(grid::rasterGrob(Silu_bird),
                      xmin = 3.5, xmax = 4.5,
                      ymin = 65,  ymax = 75)
  + theme_bw() + theme_ggplot)

(p8 <-  ggplot(mammDF, aes(x = HP, y = n)) +
    geom_col(aes(fill = Red), width = 0.7) +
    labs(x = NULL, 
         y = NULL,
         title = "B")+
    
    annotation_custom(grid::rasterGrob(Silu_mamm),
                      xmin = 3.5, xmax = 4.5,
                      ymin = 15,  ymax = 17)+
    
    scale_fill_manual(values = c("grey60", "grey10"), labels = c("Threatened", "Not threatened"), name =  NULL)
  + theme_bw() + theme_ggplot + theme(legend.position = c(0.2, 0.7))
)

# Arranging in a plot ------------------------------------------------------

plots <- list(p7, p8)
grobs <- list()
widths <- list()

for (i in 1:length(plots)){
  grobs[[i]] <- ggplotGrob(plots[[i]])
  widths[[i]] <- grobs[[i]]$widths[2:5]
}

maxwidth <- do.call(grid::unit.pmax, widths)

for (i in 1:length(grobs)){
  grobs[[i]]$widths[2:5] <- as.list(maxwidth)
}

pdf(file = "Figures/Figure 2.pdf", width = 10, height = 3)

do.call("grid.arrange", c(grobs, nrow = 1, ncol = 2))

dev.off()

# Theshold 0.3
threshold <- 0.3

sum(nrow(mamm_tree_plot[mamm_tree_plot$beta_repl < threshold & mamm_tree_plot$beta_diff < threshold,])) # HP 0
sum(nrow(mamm_tree_plot[mamm_tree_plot$beta_diff >= threshold,])) # HP 1
sum(nrow(mamm_tree_plot[mamm_tree_plot$beta_repl >= threshold,])) # HP 2
sum(nrow(mamm_tree_plot[mamm_tree_plot$beta_repl >= threshold & mamm_tree_plot$beta_diff >= threshold,])) # HP 3

sum(nrow(birds_tree_plot[birds_tree_plot$beta_repl < threshold & birds_tree_plot$beta_diff < threshold,])) # HP 0
sum(nrow(birds_tree_plot[birds_tree_plot$beta_diff >= threshold,])) # HP 1
sum(nrow(birds_tree_plot[birds_tree_plot$beta_repl >= threshold,])) # HP 2
sum(nrow(birds_tree_plot[birds_tree_plot$beta_repl >= threshold & birds_tree_plot$beta_diff >= threshold,])) # HP 3

HP0_mamm <- mamm_tree_plot[mamm_tree_plot$beta_repl < threshold & mamm_tree_plot$beta_diff < threshold,]$Name_phylo
HP2_mamm <- mamm_tree_plot[mamm_tree_plot$beta_repl >= threshold,]$Name_phylo

HP0_birds <- birds_tree_plot[birds_tree_plot$beta_repl < threshold & birds_tree_plot$beta_diff < threshold,]$Name_phylo
HP2_birds <- birds_tree_plot[birds_tree_plot$beta_repl >= threshold,]$Name_phylo

# Bind all 4 df and reorganize them
All_HP <- data.frame(species = c(HP0_birds, HP0_mamm, HP2_birds, HP2_mamm),
                     HP = c(rep("HP0_birds",length(HP0_birds)),
                            rep("HP0_mamm",length(HP0_mamm)),
                            rep("HP2_birds",length(HP2_birds)),
                            rep("HP2_mamm",length(HP2_mamm))))

trait <- data.frame(Name_phylo = db$Name_phylo, Red = db$red)

trait <- distinct(trait, Name_phylo, .keep_all = TRUE)

# now I have all data with red list status:
All_HP2 <- left_join(All_HP, trait, 
                     by = c("species" = "Name_phylo")) ; names(All_HP2)

All_HP2$HP  <- as.factor(All_HP2$HP)
All_HP2$Red <- as.factor(All_HP2$Red)

# Summarizing data

Data4hp <- All_HP2 %>% group_by(HP, Red) %>% tally()
Data4hp <- as.data.frame(Data4hp)

### Further manipulation
HP <- rep(c("HP0_birds", "HP1_birds","HP2_birds","HP3_birds", "HP0_mamm", "HP1_mamm","HP2_mamm","HP3_mamm"),2)
Red <- c(rep(0,8), rep(1,8))
HP <- data.frame(HP, Red)
HP$HP <- as.factor(HP$HP)
HP$Red <- as.factor(HP$Red)
HP2 <- left_join(HP, Data4hp, by = c("HP", "Red"))

# Turn Nas in the n to Zeros:
HP2 <- HP2 %>% mutate_at(vars("n"), ~replace(., is.na(.), 0)) 

# Now I have the DF for the figures:
birdsDF <- HP2[c(1:4, 9:12) ,]  ; birdsDF <- gdata::drop.levels(birdsDF)
mammDF  <- HP2[c(5:8, 13:16) ,] ; mammDF  <- gdata::drop.levels(mammDF) 

#Rename levels
levels(birdsDF$HP) <- levels(mammDF$HP) <- c("HP0\nNo change", "HP1\nExpansion", "HP2\nShift", "HP3\nExpansion & shift")  # relevel
birdsDF$Red <- relevel(birdsDF$Red, "1")    # order levels, threatened first
mammDF$Red  <- relevel(mammDF$Red, "1")     # order levels, threatened first

#creating a fake low value for graphical reasons
birdsDF[c(2,4),3] <- 0.5
mammDF[c(2,4),3] <- 0.1

# Make plots:
(p7 <- ggplot(birdsDF, aes(x = HP, y = n)) +
    geom_col(aes(fill = Red), width = 0.7) +
    labs(x = NULL, 
         y = "Count of species",
         title = "A")+
    scale_fill_manual(values = c("grey60", "grey10")) +
    annotation_custom(grid::rasterGrob(Silu_bird),
                      xmin = 3.5, xmax = 4.5,
                      ymin = 45,  ymax = 55)
  + theme_bw() + theme_ggplot)

(p8 <-  ggplot(mammDF, aes(x = HP, y = n)) +
    geom_col(aes(fill = Red), width = 0.7) +
    labs(x = NULL, 
         y = NULL,
         title = "B")+
    
    annotation_custom(grid::rasterGrob(Silu_mamm),
                      xmin = 3.5, xmax = 4.5,
                      ymin = 13,  ymax = 15)+
    
    scale_fill_manual(values = c("grey60", "grey10"), labels = c("Threatened", "Not threatened"), name =  NULL)
  + theme_bw() + theme_ggplot + theme(legend.position = c(0.2, 0.7))
)

# Arranging in a plot ------------------------------------------------------

plots <- list(p7, p8)
grobs <- list()
widths <- list()

for (i in 1:length(plots)){
  grobs[[i]] <- ggplotGrob(plots[[i]])
  widths[[i]] <- grobs[[i]]$widths[2:5]
}

maxwidth <- do.call(grid::unit.pmax, widths)

for (i in 1:length(grobs)){
  grobs[[i]]$widths[2:5] <- as.list(maxwidth)
}

pdf(file = "Figures/Figure S7_thresold 3.pdf", width = 10, height = 3)

do.call("grid.arrange", c(grobs, nrow = 1, ncol = 2))

dev.off()

# Checking alternative thresholds --------------------------------------

# Checking how species would fall in the different hypotheses 
# using alternative thresholds

threshold <- seq(from = 0, to = 1, by = 0.01)

HP0 <- c()
HP1 <- c()
HP2 <- c()
HP3 <- c()

for(i in 1:length(threshold)) {
  
  HP0 <- append(HP0, sum(nrow(birds_tree_plot[birds_tree_plot$beta_repl < threshold[i] & birds_tree_plot$beta_diff < threshold[i],]))) # HP 0
  HP1 <- append(HP1, sum(nrow(birds_tree_plot[birds_tree_plot$beta_diff >= threshold[i],]))) # HP 1
  HP2 <- append(HP2, sum(nrow(birds_tree_plot[birds_tree_plot$beta_repl >= threshold[i],]))) # HP 2
  HP3 <- append(HP3, sum(nrow(birds_tree_plot[birds_tree_plot$beta_repl >= threshold[i] & birds_tree_plot$beta_diff >= threshold[i],]))) # HP 3
  
}

HP <- data.frame(threshold = rep(threshold,4),
                 n         = c(HP0,HP1,HP2,HP3),
                 ID        = c(rep("HP0", length(threshold)),
                               rep("HP1", length(threshold)),
                               rep("HP2", length(threshold)),
                               rep("HP3", length(threshold))))

(p9 <- ggplot(HP, aes(x = threshold, y = n)) + facet_grid(~ ID) +
    geom_point(col="grey20", fill = "grey30") +
    scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
    labs(x = "Threshold", 
         y = "Number of species") + theme_bw() + theme_ggplot)

# and plot

pdf(file = "Figures/Figure S5.pdf", width = 10, height = 4)

p9

dev.off()

# Checking the difference in volume and beta withing and outside ----------

(p1 <- ggplot(data = db[db$Group == "bird" & db$naxes > 5,],  #if only 7 axes are shown, we lose year 2000s
              aes(x = year, y = delta_vol)) +
   geom_point(aes(x = year, y = delta_vol), 
              position = position_jitter(width = 0.15), size = 1, alpha = 0.2) +
   geom_boxplot(width = 0.4, outlier.shape = NA, alpha = 0.4, fill="grey50") +
   geom_hline(yintercept=0,col="blue",lty="dotted")+
   labs(title = "A", subtitle = "Birds",
        y = "Delta volume\n(Volume in protected - unprotected)", x = NULL) +
   guides(fill = "none", color = "none") +
   theme_classic() + theme_ggplot
)

(p2 <- ggplot(data = db[db$Group == "mammal" & db$naxes > 6,], 
              aes(x = year, y = delta_vol)) +
    geom_point(aes(x = year, y = delta_vol), 
               position = position_jitter(width = 0.15), size = 1, alpha = 0.2) +
    geom_boxplot(width = 0.4, outlier.shape = NA, alpha = 0.4, fill="grey50") +
    geom_hline(yintercept=0,col="blue",lty="dotted")+
    labs(title = "B", subtitle = "Mammal", 
         y = "Delta volume\n(Volume in protected - unprotected)", x = NULL) +
    
    guides(fill = "none", color = "none") +
    theme_classic() + theme_ggplot
)

## Supplementary Figures S4
plots <- list(p1, p2) ; grobs <- list() ; widths <- list()

for (i in 1:length(plots)){
  grobs[[i]] <- ggplotGrob(plots[[i]])
  widths[[i]] <- grobs[[i]]$widths[2:5]
}

maxwidth <- do.call(grid::unit.pmax, widths)

for (i in 1:length(grobs)){
  grobs[[i]]$widths[2:5] <- as.list(maxwidth)
}

pdf("Figures/Figure S4.pdf", width = 10, height = 4)
do.call("grid.arrange", c(grobs, nrow = 1, ncol = 2))
dev.off()

# End 