## Authors: LTER pipeline team
## Purpose: build a general function to analyze popler data sets and generate desired outputs
## Last update: February 1, 2019

## install popler
#install.packages("devtools")
#devtools::install_github("AldoCompagnoni/popler", build_vignettes=T, force=T)
#install.packages("dplyr")
library(popler)
library(tidyverse)
library(rstan)

## let's use the heron data set as a guinea pig
herons_metadat <- pplr_browse(proj_metadata_key==88)
herons <- pplr_get_data(herons_metadat,cov_unpack=T)
herons_metadat$lat_lter

bigfun <- function(popler_proj_key,rarity_threshold){
  
  ## extract popler project data and metadata
  metadat <- pplr_browse(proj_metadata_key==popler_proj_key, full_tbl = T)
  ## diagnose the data type
  type <- metadat$datatype
  if(type=="individual" | type=="basal_area"){return("Non-desired data type")}
  ## get data and combine spatial rep info
  n_spat_levels <- metadat$n_spat_levs
  dat <- pplr_get_data(metadat) 
  dat <- dat %>% mutate(ran_effect = ifelse(n_spat_levels==1,spatial_replication_level_1,
                          ifelse(n_spat_levels==2,interaction(spatial_replication_level_1,spatial_replication_level_2),
                                 ifelse(n_spat_levels==3,interaction(spatial_replication_level_1,spatial_replication_level_2,spatial_replication_level_3),
                                        interaction(spatial_replication_level_1,spatial_replication_level_2,spatial_replication_level_3,spatial_replication_level_4)))))

  ## filter out NAs and very rare species -- still need to do

  ## number of species we are left with
  spp <- na.omit(unique(dat$sppcode))
  n_spp <- length(spp)

  ## prep data for analysis
  datalist<-list(
    n=nrow(dat),
    nyear=length(unique(as.factor(dat$year))),
    nrep=length(unique(as.factor(dat$ran_effect))),
    year=as.numeric(as.factor(as.character(dat$year))),
    rep=as.numeric(as.factor(as.character(dat$ran_effect))),
    count=dat$abundance_observation)
  
  ## send to the appropriate Stan model, given data type
  stan_model <- ifelse(type=="count","count_model.stan",
                       ifelse(type=="biomass","biomass_model.stan",
                              ifelse(type=="cover","cover_model.stan","density_model.stan")))
  abund_fit<-stan(file=stan_model,data=datalist,iter=10000,chains=3)
  
  
  
  ## collect climate covariates
  latlong_DD <- c(metadat$lat_lter, metadat$long_lter)
  years <- metadat$studystartyr:metadat$studyendyr
  
  ## pull out items of interest from Stan output - could embed climate change simulation here
  
  ## package outputs into data frame or list
  
}

## pseudo-code for Stan model
#1. Get year-specific lambdas
#2. Fit spline for climate covariate
#3. Derive expected value
#4. Calculate model diagnostics




# trash

dat %>% 
  group_by(sppcode) %>% 
  select(abundance_observation) %>% 
  filter(abundance_observation > 0) %>% 
  summarise(n())

test <- dat %>% 
  group_by(sppcode) %>% 
  filter(abundance_observation > 0,
         !is.na(abundance_observation))

filter(test,sppcode=="WI")
spp <- na.omit(unique(test$sppcode))
