## ------------------------------------------------------------------------
## 'Effects of protected areas on ecological niche properties of terrestrial vertebrates'
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

# Loading silhouettes for figures ------------------------------------------

# Taken from PhyloPic (http://phylopic.org/) with open license.

Silu_bird   <- png::readPNG("Silhouettes/Silu_turdus.png") 
Silu_mamm   <- png::readPNG("Silhouettes/Silu_vulpes.png")

# Loading the databases ---------------------------------------------------

# Read the hypervolume results that includes the traits:
db <- read.csv("Data/hypervolume_metrics_TRAITS_15_7_21.csv", header = TRUE, sep = ";", dec = ",",row.names = NULL, as.is = FALSE)

#rename columns
colnames(db)[1]  <- "Name_phylo"
colnames(db)[16] <- "mass"
colnames(db)[19] <- "red"
colnames(db)[20] <- "hab"
colnames(db)[21] <- "diet"
colnames(db)[24] <- "lat"

#checking the strcuture
str(db)

#converting numeric to factor
db$red <- as.factor(db$red)

#Renaming the variable Name_phylo
db$Name_phylo <- gsub(' ', '_', db$Name_phylo)

# Loading Phylogenetic trees ----------------------------------------------

mammTREE <- read.tree("Data/mammTree_15_7_21.tre")
birdTREE <- read.tree("Data/birdTree_15_7_21.tre")

# Checking the difference in volume and beta withing and outside ----------

plot <- db[db$naxes == 7,]

(p1 <- ggplot(data = plot[plot$Group == "bird",], 
              aes(x = year, y = delta_vol)) +
    geom_point(aes(x = year, y = delta_vol), 
               position = position_jitter(width = 0.15), size = 1, alpha = 0.2) +
    geom_boxplot(width = 0.4, outlier.shape = NA, alpha = 0.4, fill="grey50") +
      geom_hline(yintercept=0,col="blue",fill="gray60",lty="dotted")+
      labs(title = "A", subtitle = "Birds",
             y = "Delta volume\n(Volume in protected - unprotected)", x = NULL) +
         guides(fill = FALSE, color = FALSE) +
      theme_classic() + theme_ggplot
)

(p2 <- ggplot(data = plot[plot$Group == "mammal",], 
              aes(x = year, y = delta_vol)) +
    geom_point(aes(x = year, y = delta_vol), 
               position = position_jitter(width = 0.15), size = 1, alpha = 0.2) +
    geom_boxplot(width = 0.4, outlier.shape = NA, alpha = 0.4, fill="grey50") +
    geom_hline(yintercept=0,col="blue",fill="gray60",lty="dotted")+
    labs(title = "B", subtitle = "Mammal", 
         y = "Delta volume\n(Volume in protected - unprotected)", x = NULL) +
    
   guides(fill = FALSE, color = FALSE) +
    theme_classic() + theme_ggplot
)

## Supplementary Figures S3
plots <- list(p1, p2) ; grobs <- list() ; widths <- list()

for (i in 1:length(plots)){
  grobs[[i]] <- ggplotGrob(plots[[i]])
  widths[[i]] <- grobs[[i]]$widths[2:5]
}

maxwidth <- do.call(grid::unit.pmax, widths)

for (i in 1:length(grobs)){
  grobs[[i]]$widths[2:5] <- as.list(maxwidth)
}

pdf("Figures/Figure_S3.pdf", width = 10, height = 4)
do.call("grid.arrange", c(grobs, nrow = 1, ncol = 2))
dev.off()

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
dotchart(birds$mass, main= "mass") #log transform!
dotchart(birds$lat, main= "latitude") #ok
dotchart(birds$Gen_Length, main= "Gen_length") #log trasform!

# trasforming variables
birds$mass_log   <-  log(birds$mass+1)
birds$Gen_log    <-  log(birds$Gen_Length+1)
birds$diet_asin  <-  asin(birds$diet)

#double check
dev.off() ; par(mfrow = c(2,2),mar = c(rep(2,4)))
dotchart(birds$mass_log, main= "mass_log") #ok
dotchart(birds$Gen_log , main= "Gen_log ") #ok
dotchart(birds$diet_asin , main= "diet_asin") #ok

# Check collinearity
psych::pairs.panels(birds[,colnames(birds)[c(20,24,27:29)]]) 
#mass and generation are collinear! We select only mass

#Checking association with categorical
boxplot(mass_log ~ red, data = birds)
boxplot(hab ~ red, data = birds)
boxplot(diet_asin ~ red, data = birds)
boxplot(lat ~ red, data = birds)
#Ok

# Check relationships
psych::pairs.panels(birds[,colnames(birds)[c(6,7,20,24,27,29)]]) 

# GLMM ----------------------------------------

##Setting baseline for year
birds <- within(birds, year <- relevel(year, ref = "2018s"))

# scale and center all continuous variables before analyses
birds$mass_log2  <-  scale(birds$mass_log, center = TRUE, scale = TRUE)
birds$hab        <-  scale(birds$hab, center = TRUE, scale = TRUE)
birds$diet_asin  <-  scale(birds$diet_asin, center = TRUE, scale = TRUE)
birds$lat        <-  scale(birds$lat, center = TRUE, scale = TRUE)
 
#creating a new dataframe
birds_model <- data.frame(species = birds$species,
                 name_phylo = birds$Name_phylo,
                 order = birds$Order,
                 family = birds$Family,
                 year    = birds$year,
                 beta_repl = birds$beta_repl,
                 beta_diff = birds$beta_diff,
                 red = birds$red,
                 mass = birds$mass_log2,
                 mass_log = birds$mass_log,
                 diet = birds$diet_asin,
                 vol = birds$delta_vol,
                 hab = birds$hab,
                 lat = birds$lat,
                 naxes = birds$naxes)

birds_model <- na.omit(birds_model)

# Model for replacement
bird_M1 <- glmmTMB(beta_repl ~ mass + hab + diet + lat + red + (1 | year) + (1 | species),
                   data = birds_model, beta_family(link = "logit"))

#Checking the model
parameters::model_parameters(bird_M1)
performance::check_model(bird_M1)

#Checking for phylogenetic signal in the residuals
res <- residuals(bird_M1)
birds_test <- data.frame(birds_model, 
                      residuals = residuals(bird_M1)) 

birds_test <- birds_test %>% group_by(name_phylo) %>% summarise_at(vars(c(13)), funs(mean(., na.rm=TRUE)))

birds_test <- data.frame(birds_test)
rownames(birds_test) <- make.names(birds_test$name_phylo, unique = TRUE)

#checking match of names in the phylogenetic tree (and dropping unused tips)
rownames(birds_test) <- make.names(birds_test$name_phylo, unique = TRUE)

(missing.tip <- birdTREE$tip.label[birdTREE$tip.label %in% birds_test$name_phylo == FALSE]) # 7 species to drop

birdTREE2 <- drop.tip(birdTREE, missing.tip)

unique(birdTREE2$tip.label %in% birds_test$name_phylo) # all good now

#Checking signal in the residuals
BirdTreeList <- birdTREE2$tip.label      # get order of the sp in tree
birds_test   <- birds_test[order(match(rownames(birds_test), BirdTreeList)),] # reorder the DF

phylosig(birdTREE2, birds_test$residuals, method = "lambda", test = TRUE)

# Phylogenetic signal lambda : 6.67426e-05 
# logL(lambda) : 74.5521 
# LR(lambda=0) : -0.00462835 
# P-value (based on LR test) : 1 

# No phylo signal in the residuals!

(p3 <- sjPlot::plot_model(bird_M1, sort.est = FALSE, se = TRUE,col="black",
                   vline.color ="grey70",
                   title = "C - Habitat shift",
                   show.values = TRUE, value.offset = .3,
                   axis.labels = c("Red list status\n[threatened]","Mean\nlatitude","Diet\nspecialization",
                                   "Habitat\nspecialization","Body mass")) +
    
                  annotation_custom(grid::rasterGrob(Silu_bird),
                                    xmin = unit(1, "native"), xmax = unit(1.8,"native"),
                                    ymin = unit(0.3,"npc"),  ymax = unit(1,"npc"))+ 
  
               theme_bw() + theme_ggplot)

# Model for expantion/contraction

#removing the outlier
birds_model2  <- birds_model[birds_model$beta_diff < 0.2,]
birds_model2  <- na.omit(birds_model2)

#Doesn't work!
bird_M2 <- glmmTMB(beta_diff ~ mass + hab + diet + lat + red + (1 | year) + (1 | species), 
                   data = birds_model2, beta_family(link = "logit"))

#Error messagge: there are 5 value = 0 in the response variable.

range(db$beta_diff) 
table(db$beta_diff)

#replacing the value with a very small number
birds_model2$beta_diff <- ifelse(birds_model2$beta_diff == 0, 
                                 birds_model2$beta_diff + 0.000000000000000000001,
                                 birds_model2$beta_diff)

#Re-fitting the model
bird_M2 <- glmmTMB(beta_diff ~ mass + hab + diet + lat + red + (1 | year) + (1 | species), 
                   data = birds_model2, beta_family(link = "logit"))

#Checking the model
performance::check_model(bird_M2)
parameters::model_parameters(bird_M2)

#Checking for phylogenetic signal in the residuals
res <- residuals(bird_M2)
birds_test <- data.frame(birds_model2, 
                         residuals = residuals(bird_M2)) 

birds_test <- birds_test %>% group_by(name_phylo) %>% summarise_at(vars(c(13)), funs(mean(., na.rm=TRUE)))

birds_test <- data.frame(birds_test)

unique(birdTREE2$tip.label %in% birds_test$name_phylo) # all good now

#Checking the signal in the residuals
BirdTreeList <- birdTREE2$tip.label      # get order of the sp in tree
birds_test <- birds_test[order(match(rownames(birds_test), BirdTreeList)),] # reorder the DF

phylosig(birdTREE2, birds_test$residuals, method = "lambda", test = TRUE)

# hylogenetic signal lambda : 6.67426e-05 
# logL(lambda) : 340.327 
# LR(lambda=0) : -0.00158656 
# P-value (based on LR test) : 1 

# No phylo signal in the residuals! 

(p1 <- sjPlot::plot_model(bird_M2, sort.est = FALSE, se = TRUE, col="black",
                   vline.color ="grey70",
                   title = "A - Niche expansion",
                   show.values = TRUE, value.offset = .3, axis.title = " ",
                   axis.labels = c("Red list status\n[threatened]","Mean\nlatitude","Diet\nspecialization",
                                   "Habitat\nspecialization","Body mass")) +
                    
                     annotation_custom(grid::rasterGrob(Silu_bird),
                      xmin = unit(1, "native"), xmax = unit(1.8,"native"),
                      ymin = unit(0.3,"npc"),  ymax = unit(1,"npc"))+ 
  
                   theme_bw() + theme_ggplot)

## Check mass result --- Figure S4

birds_model2$expansion <- ifelse(birds_model2$vol > 0,"Contraction outside protected areas", "Expansion outside protected areas")

s4<- ggplot(birds_model2, aes(x=mass, fill=expansion, color=expansion)) +
  geom_density(alpha =.5) +  labs(x = "Logarithm of body mass (scaled)" , y = "Density")+
  scale_fill_manual(values =  c( "orange","grey20")) +
  scale_color_manual(values =  c("orange","grey20")) + theme_bw() + theme_ggplot + theme(legend.title = element_blank(),legend.position = c(0.5, 0.7))

pdf("Figures/Figure_S4.pdf", width = 8, height = 4)
s4
dev.off()

# Analysis with mammals -----------------------------------------------------

mamm <- db[db$Group == "mammal",]

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
dotchart(mamm$lat, main= "latitude") #ok
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
psych::pairs.panels(mamm[,colnames(mamm)[c(17,24,27:29)]]) 
#mass and generation are collinear! We select only mass

#Checking association with categorical
boxplot(mass_log ~ red, data = mamm)
boxplot(hab_asin ~ red, data = mamm)
boxplot(diet_asin ~ red, data = mamm)
boxplot(lat ~ red, data = mamm)
#Ok

# Check relationships
psych::pairs.panels(mamm[,colnames(mamm)[c(6,7,24,27:29)]]) 

# GLMM ----------------------------------------

##Setting baseline for year
mamm <- within(mamm, year <- relevel(year, ref = "2018s"))

# scale and center all continuous variables before analyses
mamm$mass_log  <-  scale(mamm$mass_log, center = TRUE, scale = TRUE)
mamm$hab       <-  scale(mamm$hab_asin, center = TRUE, scale = TRUE)
mamm$diet_asin <-  scale(mamm$diet_asin, center = TRUE, scale = TRUE)
mamm$lat       <-  scale(mamm$lat, center = TRUE, scale = TRUE)
 
#creating a new dataframe
mamm_model <- data.frame(species = mamm$species,
                         name_phylo = mamm$Name_phylo,
                         order = mamm$Order,
                         family = mamm$Family,
                         year    = mamm$year,
                         beta_repl = mamm$beta_repl,
                         beta_diff = mamm$beta_diff,
                         red = mamm$red,
                         mass = mamm$mass_log,
                         diet =mamm$diet_asin,
                         hab = mamm$hab,
                         lat = mamm$lat)

mamm_model <- na.omit(mamm_model)

mamm_model1 <- mamm_model$beta_repl

# Model for replacement
mamm_M1 <- glmmTMB(beta_repl ~ mass + hab + diet + lat + red + (1|year) + (1|species), data = mamm_model, beta_family(link = "logit"))

#Checking the model
parameters::model_parameters(mamm_M1)
performance::check_model(mamm_M1)

#Checking for phylogenetic signal in the residuals
res <- residuals(mamm_M1)
mamm_test <- data.frame(mamm_model, 
                        residuals = residuals(mamm_M1)) 

mamm_test <- mamm_test %>% group_by(name_phylo) %>% summarise_at(vars(c(12)), funs(mean(., na.rm=TRUE)))

mamm_test <- data.frame(mamm_test)
rownames(mamm_test) <- make.names(mamm_test$name_phylo, unique = TRUE)

# Checking match of names in the phylogenetic tree 
unique(mammTREE$tip.label %in% mamm_test$name_phylo) # all good now

MammTreeList <- mammTREE$tip.label      # get order of the sp in tree
mamm_test    <- mamm_test[order(match(rownames(mamm_test), MammTreeList)),] # reorder the DF

phylosig(mammTREE, mamm_test$residuals, method = "lambda", test = TRUE)

# Phylogenetic signal lambda : 6.95863e-05 
# logL(lambda) : 14.2156 
# LR(lambda=0) : -0.00100615 
# P-value (based on LR test) : 1 

# No phylo signal in the residuals! 

(p4 <- sjPlot::plot_model(mamm_M1, sort.est = FALSE, se = TRUE, col = "black",
                   vline.color ="grey70",
                   title = "D - Habitat shift",
                   show.values = TRUE, value.offset = .3,
                   axis.labels = c(rep(" ",5))) + theme_bw() + 
    
                    annotation_custom(grid::rasterGrob(Silu_mamm),
                      xmin = unit(1, "native"), xmax = unit(1.8,"native"),
                      ymin = unit(0.4,"npc"),  ymax = unit(1,"npc"))+ 
    
                   theme_ggplot)

#removing the outlier
mamm_model2 <- mamm_model[mamm_model$beta_diff < 0.12,]

mamm_M2 <- glmmTMB(beta_diff ~ mass + hab + diet + lat + red +(1|year) + (1|species), data = mamm_model2, beta_family(link = "logit"))

#Checking the model
parameters::model_parameters(mamm_M2)
performance::check_model(mamm_M2)

#Checking for phylogenetic signal in the residuals
res <- residuals(mamm_M2)
mamm_test <- data.frame(mamm_model2, 
                        residuals = residuals(mamm_M2)) 

mamm_test <- mamm_test %>% group_by(name_phylo) %>% summarise_at(vars(c(12)), funs(mean(., na.rm=TRUE)))

mamm_test <- data.frame(mamm_test)
rownames(mamm_test) <- make.names(mamm_test$name_phylo, unique = TRUE)

MammTreeList <- mammTREE$tip.label      # get order of the sp in tree
mamm_test <- mamm_test[order(match(rownames(mamm_test), MammTreeList)),] # reorder the DF

phylosig(mammTREE, mamm_test$residuals, method = "lambda", test = TRUE)

# Phylogenetic signal lambda : 6.95863e-05 
# logL(lambda) : 62.7059 
# LR(lambda=0) : -0.0011535 
# P-value (based on LR test) : 1  

# No phylo signal in the residuals!

(p2 <- sjPlot::plot_model(mamm_M2, sort.est = FALSE, se = TRUE, col= "black",
                   vline.color ="grey70",
                   title = "B - Niche expansion", axis.title = " ",
                   show.values = TRUE, value.offset = .3,
                   axis.labels = c(rep(" ",5))) +
    
                   annotation_custom(grid::rasterGrob(Silu_mamm),
                      xmin = unit(1, "native"), xmax = unit(1.8,"native"),
                      ymin = unit(0.4,"npc"),  ymax = unit(1,"npc"))+ 
                  
                    theme_bw() + theme_ggplot)

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

mamm_tree_plot <- mamm %>% group_by(Name_phylo) %>% summarise_at(vars(c(4,5,6)), funs(mean(., na.rm=TRUE)))

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

birds_tree_plot <- birds %>% group_by(Name_phylo) %>% summarise_at(vars(c(4,5,6)), funs(mean(., na.rm=TRUE)))

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
                      ymin = 70,  ymax = 80)
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

pdf(file = "Figures/Figure_S5.pdf", width = 10, height = 4)

p9

dev.off()

