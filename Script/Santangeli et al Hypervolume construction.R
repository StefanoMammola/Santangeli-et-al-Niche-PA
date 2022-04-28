## ------------------------------------------------------------------------
## 'The effects of protected areas on the ecological niches of birds and mammals'
## ------------------------------------------------------------------------

# Santangeli et al.

## ------------------------------------------------------------------------
# 'R script to reproduce the analyses'
## ------------------------------------------------------------------------

# Analysis performed with R (v. R 4.0.3) and R studio (v. 1.4.1103)
# Authors: Stefano Mammola
# Location: Helsinki, April 2021

# clean the workspace -----------------------------------------------------

rm(list=ls())

# Loading R package -------------------------------------------------------

library("BAT")
library("dplyr")
library("ggplot2")
library("gridExtra")
library("hypervolume")
library("pbmcapply")
library("tidyr")

# Loading the databases ---------------------------------------------------

#Database is just an example with 10 species; the full database is not provided because it exceed the size limit for GitHub. 

db2  <- read.csv("Data/hypervolume_example_database.csv", sep=",", header = TRUE, as.is = TRUE)

# Parallels ---------------------------------------------------------------

# Settings for parallelisation
cores <- detectCores()-2

# Selecting only matched areas --------------------------------------------

db <- db2 %>% drop_na(matched) ; db <- select(db, -c(matched))

# Preparing the database --------------------------------------------------

# converting chr to factors
db$scientificName    <- as.factor(db$scientificName)
db$PA                <- as.factor(db$PA)
db$CLC_yr            <- as.factor(db$CLC_yr)

# renaming PA
levels(db$PA) <- c("unprotected","protected")

# reducing the number of axes
db <- db %>% mutate(artif = pr_artif) ; db <- select(db, -c(pr_artif))
db <- db %>% mutate(shrub = pr_shrub_herb) ; db <- select(db, -c(pr_shrub_herb))
db <- db %>% mutate(water = pr_wetland + pr_water) ; db <- select(db, -c(pr_wetland, pr_water))
db <- db %>% mutate(agro  = pr_arable + pr_pasture + pr_heter_agric) ; db <- select(db, -c(pr_arable, pr_pasture, pr_heter_agric))
db <- db %>% mutate(forest_broad = pr_forest_broadl) ; db <- select(db, -c(pr_forest_broadl))
db <- db %>% mutate(forest_conif = pr_forest_conif) ; db <- select(db, -c(pr_forest_conif))
db <- db %>% mutate(forest_mix = pr_forest_mixed) ; db <- select(db, -c(pr_forest_mixed))

# we are removing others & open

#  Selecting only counts > 0 & species present in both areas -------------------

db <- db[db$count>0,] ; db$scientificName <- droplevels(db$scientificName)

db_dist <- distinct(db, scientificName, CLC_yr, PA) 
db_count <- db_dist %>% group_by(scientificName, CLC_yr) %>% tally()

db_count <- db_count[db_count$n > 1,]

db <- left_join(db,db_count, by = c("scientificName", "CLC_yr"))
db <- db %>% drop_na(n) ; db <- select(db, -c(n))
db <- droplevels(db)

# Filtering out species occurring in less than 5 sites/year/PA ------------

# change count to occurrence:
db_count <- db %>%
  mutate(Occurrence = case_when(
    count > 0 ~ 1,
  ))

# calc the sum of counts per sp and route and yr
db_count <- db_count %>%
  group_by(scientificName, CLC_yr, PA) %>%
  summarise_at(vars(Occurrence), sum)

db_count2 <- db_count[db_count$Occurrence > 4, ]

db_count2 <- db_count2 %>%
  group_by(scientificName, CLC_yr) %>%
  tally()

db <- left_join(db, db_count2, by = c("scientificName", "CLC_yr"))
db <- db[db$n > 1, ] ; db <- select(db, -c(n))
db <- na.omit(db)

db$scientificName <- droplevels(as.factor(db$scientificName))

# tide up the environment -------------------------------------------------

rm(db_count,db_count2,db_dist,db2)

# Generating the hypervolumes ---------------------------------------------

str(db)

# env variable from 11 : 17
n_min = 11
n_max = 17

# empty data.frame
result_beta <-  data.frame(species           = NULL, 
                            year             = NULL,
                            naxes            = NULL,
                            beta_total       = NULL,
                            beta_repl        = NULL,
                            beta_diff        = NULL,
                            vol_protected    = NULL,
                            vol_unprotected  = NULL,
                            disp_protected   = NULL,
                            disp_unprotected = NULL,
                            delta_vol        = NULL,
                            delta_disp       = NULL)

# Number of iteration = number of species
n_iter <- unique(nlevels(db$scientificName))

# Starting the loop
for(i in 1 : n_iter) {
  
  message(paste("-- Modelling species", i ,"out of", n_iter))
  message(paste("--", unique(db$scientificName)[i]), "--")
  
  db_i <- db[db$scientificName == unique(db$scientificName)[i],]
   
  # Creating a list of observation by year in protected areas
  db_i_protected <- db_i[db_i$PA == "protected",]
  
  db_i_protected <- lapply(1:nlevels(db_i_protected$CLC_yr), 
                           function(s) db_i_protected[db_i_protected$CLC_yr == levels(db_i_protected$CLC_yr)[s], ])
  
  # Creating a list of observation by year in unprotected areas
  db_i_unprotected <- db_i[db_i$PA == "unprotected",]
  
  db_i_unprotected <- lapply(1:nlevels(db_i_unprotected$CLC_yr), 
                             function(s) db_i_unprotected[db_i_unprotected$CLC_yr == levels(db_i_unprotected$CLC_yr)[s], ])

  # procedure to remove hypervolume axis variable when one of the habitat is 0
  remove <- c()  
 
  for(h in 1:length(db_i_protected)) {
    
    db_h_protected <- db_i_protected[[h]]
    db_h_unprotected <- db_i_unprotected[[h]]
    
    omit <- ifelse(colSums(db_h_protected[,n_min : n_max]) == 0, colnames(db_h_protected[,n_min : n_max]), NA)
    omit <- as.character(na.omit(as.factor(omit)))
    
    omit2 <- ifelse(colSums(db_h_unprotected[,n_min : n_max]) == 0, colnames(db_h_protected[,n_min : n_max]), NA)
    omit2 <- as.character(na.omit(as.factor(omit2)))
    
    omit <- unique(c(omit,omit2))
    
    if(length(omit)==7) {
      
    remove <- c(remove,h)
  
    } else {
      
      db_i_protected[[h]]   <- select(db_h_protected, -omit)
      db_i_unprotected[[h]] <- select(db_h_protected, -omit)
      
    }
    
  }
  
  # Remove null hypervolumes if needed
  if(length(remove) > 0) {
    db_i_protected   <- db_i_protected[-c(remove)] ; db_i_unprotected <- db_i_unprotected[-c(remove)]    
  }
    
  #estimating bandwidth
  bw  <- lapply(1:length(db_i_protected), function(x) estimate_bandwidth(db_i_protected[[x]][,n_min : ncol(db_i_protected[[x]])]))
  bw2 <- lapply(1:length(db_i_unprotected), function(x) estimate_bandwidth(db_i_unprotected[[x]][,n_min : ncol(db_i_protected[[x]])]))
  

  #Creating the hypervolumes ---> protected areas
  message("----- Creating hypervolumes for protected areas")
  
  hv_protected <- do.call(hypervolume_join, pbmclapply(1:length(db_i_protected), 
                                                       function(x) hypervolume_gaussian(db_i_protected[[x]][,n_min : ncol(db_i_protected[[x]])],
                                                                                        weight = db_i_protected[[x]]$count/sum(db_i_protected[[x]]$count), #relative abundance
                                                                                        verbose = FALSE,
                                                                                        kde.bandwidth = bw[[x]],
                                                                                        name = as.character(unique(db_i_protected[[x]]$CLC_yr))),
                                                       mc.cores = cores,
                                                       mc.style = "ETA"))
  
  
  # Creating the hypervolumes ---> unprotected areas
  message("----- Creating hypervolumes for unprotected areas")
  
  hv_unprotected <- do.call(hypervolume_join, pbmclapply(1:length(db_i_unprotected), 
                                                         function(x) hypervolume_gaussian(db_i_unprotected[[x]][,n_min : ncol(db_i_protected[[x]])],
                                                                                          weight = db_i_unprotected[[x]]$count/sum(db_i_unprotected[[x]]$count), #relative abundance
                                                                                          verbose = FALSE,
                                                                                          kde.bandwidth = bw2[[x]],
                                                                                          name = as.character(unique(db_i_unprotected[[x]]$CLC_yr))),
                                                         mc.cores = cores,
                                                         mc.style = "ETA"))
  
  
  # Calculating beta
  message("----- Calculating hypervolume metrics")
  
  beta_tot         <- c()
  beta_repl        <- c()
  beta_diff        <- c()
  vol_protected    <- c()
  vol_unprotected  <- c()
  disp_protected   <- c()
  disp_unprotected <- c()
  delta_vol        <- c()
  delta_disp       <- c()
  year             <- c()
  naxes            <- c()
  
  for(j in 1:length(hv_protected@HVList)) {
    
    beta <- BAT::kernel.beta(hypervolume_join(hv_protected@HVList[[j]],hv_unprotected@HVList[[j]]))
    
    beta_tot         <- c(beta_tot,  beta[[1]])
    beta_repl        <- c(beta_repl, beta[[2]])
    beta_diff        <- c(beta_diff, beta[[3]])
    vol_protected    <- c(vol_protected,   hv_protected@HVList[[j]]@Volume*100000)
    vol_unprotected  <- c(vol_unprotected, hv_unprotected@HVList[[j]]@Volume*100000)
    disp_protected   <- c(disp_protected,   BAT::kernel.dispersion(hv_protected[[j]],frac=0.02))
    disp_unprotected <- c(disp_unprotected, BAT::kernel.dispersion(hv_unprotected[[j]],frac=0.02))
    delta_vol        <- c(delta_vol, (vol_protected[j]-vol_unprotected[j]))
    delta_disp       <- c(delta_disp, (disp_protected[j]-disp_unprotected[j]))
    
    year             <- c(year,hv_protected@HVList[[j]]@Name)
    naxes            <- c(naxes,hv_protected@HVList[[j]]@Dimensionality)
    
  }
  
  result_beta2 <-  data.frame(species = rep(unique(db$scientificName)[i],length(hv_protected@HVList)), 
                              year         = year,
                              naxes        = naxes,
                              beta_total   = beta_tot,
                              beta_repl    = beta_repl,
                              beta_diff    = beta_diff,
                              vol_protected    = vol_protected,
                              vol_unprotected  = vol_unprotected,
                              disp_protected   = disp_protected, 
                              disp_unprotected = disp_unprotected,
                              delta_vol    = delta_vol,
                              delta_disp   = delta_disp)
  
  result_beta <- rbind(result_beta,result_beta2)
  
}

warnings() # these rows have to be omitted.

# Saving the results ------------------------------------------------------

write.csv(result_beta, file="hypervolume_metrics_15_7_21.csv", row.names = FALSE)
