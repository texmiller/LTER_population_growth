## Authors: LTER pipeline team
## Purpose: build a general function to analyze popler data sets and generate desired outputs
## Last update: February 13, 2019

## install popler
#install.packages("devtools")
#devtools::install_github("AldoCompagnoni/popler", build_vignettes=T, force=T)
#install.packages("dplyr")
library(popler)
library(tidyverse)
library(R2jags)
library(rstan)
library(bayesplot)
library(rstanarm)

## let's use the heron data set (88) as a guinea pig
k <- 88

pipeline <- function(k){
  
  ## extract popler project data and metadata
  metadat <- pplr_browse(proj_metadata_key==as.integer(k), full_tbl = T)
  ## diagnose the data type
  type <- metadat$datatype
  ## break out of function if datatype is individual or basal cover
  if(type=="individual" | type=="basal_cover"){return("Non-desired data type")}
  
  ## get data and combine spatial rep info
  n_spat_levels <- metadat$n_spat_levs
  dat <- pplr_get_data(metadat,cov_unpack = T) %>% 
    as.data.frame %>% 
    mutate(n_spat_levels = n_spat_levels) %>% 
    mutate(ran_effect = ifelse(n_spat_levels==1,spatial_replication_level_1,
                               ifelse(n_spat_levels==2,interaction(spatial_replication_level_1,spatial_replication_level_2),
                                      ifelse(n_spat_levels==3,interaction(spatial_replication_level_1,spatial_replication_level_2,spatial_replication_level_3),
                                             interaction(spatial_replication_level_1,spatial_replication_level_2,spatial_replication_level_3,spatial_replication_level_4))))) %>% 
    filter(!is.na(abundance_observation))
  ## filter out NAs and very rare species -- we made a decision to use only the data provided by PIs-- we are not
  ## assumming that NAs are zero
  
  ## if study is experimental, use only the control group
  ## check with Aldo what to specify here
  if(metadat$studytype=="exp"){}
  
  ## keep track of what year*ran_effect levels were lost by na.omit
  summ <- dat %>% 
    group_by(year,ran_effect) %>% 
    summarise(n(),sum(is.na(abundance_observation)))

  ## prep data for analysis
  newdat <- dat %>% 
    select(year,ran_effect,sppcode,abundance_observation) %>% 
    drop_na()
  
  datalist<-list(
    n=nrow(newdat),
    nyear=length(unique(as.factor(newdat$year))),
    nrep=length(unique(as.factor(newdat$ran_effect))),
    nsp=length(unique(as.factor(newdat$sppcode))),
    year=as.numeric(as.factor(as.character(newdat$year))),
    rep=as.numeric(as.factor(as.character(newdat$ran_effect))),
    sp=as.numeric(as.factor(as.character(newdat$sppcode))),
    count=newdat$abundance_observation)
  
  inits<-function(){list(a = matrix(rnorm(datalist$nyear*datalist$nsp,0,1),nrow=datalist$nyear,ncol=datalist$nsp),
                         sigma.rep = rlnorm(1),sigma.od = rlnorm(1))}
  parameters<-c("a","sigma.rep","sigma.od","r","fit","fit.new")
    
  ## send to the appropriate JAGS model, given data type
  jags_model <- ifelse(type=="count","count_model_jags.txt",
                       ifelse(type=="biomass","biomass_model_jags.txt",
                              ifelse(type=="cover","cover_model_jags.txt","density_model_jags.txt")))
  abund_fit<-jags(data=datalist,
                  inits=inits,
                  parameters.to.save=parameters,
                  model.file="count_model_jags.txt",
                  n.thin=5,
                  n.chains=3,
                  n.burnin=1000,
                  n.iter=5000,
                  working.directory=getwd())
  
  ## Bayesian p-value
  bayes.p <- mean(abund_fit$BUGSoutput$sims.list$fit > abund_fit$BUGSoutput$sims.list$fit.new)

  ## collect posterior mean and CI growth rates, convert to year * species matrix
  r.indices <- which(str_sub(rownames(abund_fit$BUGSoutput$summary),1,1)=="r")
  r.out <- tibble(r_mean = abund_fit$BUGSoutput$summary[r.indices,"mean"],
                  r_lowCI = abund_fit$BUGSoutput$summary[r.indices,"2.5%"],
                  r_highCI = abund_fit$BUGSoutput$summary[r.indices,"97.5%"],
                  year = rep(1:(datalist$nyear-1),times=datalist$nsp),
                  sp = rep(1:datalist$nsp,each=(datalist$nyear-1)),
                  r_Rhats = abund_fit$BUGSoutput$summary[r.indices,"Rhat"])
  
  r.plot <- ggplot(r.out)+
    geom_point(aes(x=as.factor(year),y=r_mean))+
    geom_errorbar(aes(x=as.factor(year),ymin=r_lowCI, ymax=r_highCI))+
    facet_wrap(~sp)

  return(list(metadat=metadat,
              r.out=r.out,
              r.plot=r.plot,
              bayes.p=bayes.p))
}

test <- pipeline(88)

r.out %>% filter(sp!=8) %>% 
ggplot()+
  geom_point(aes(x=as.factor(year),y=r_mean))+
  geom_errorbar(aes(x=as.factor(year),ymin=r_lowCI, ymax=r_highCI))+
  facet_wrap(~sp)













## pseudo-code for Stan model
#1. Get year-specific lambdas
#2. Fit spline for climate covariate
#3. Derive expected value
#4. Calculate model diagnostics


traceplot(abund_fit)
print(abund_fit)

mcmc_areas(
  posterior, 
  pars = c("cyl", "drat", "am", "sigma"),
  prob = 0.8, # 80% intervals
  prob_outer = 0.99, # 99%
  point_est = "mean"
)



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


## collect climate covariates
latlong_DD <- c(metadat$lat_lter, metadat$long_lter)
years <- metadat$studystartyr:metadat$studyendyr

## pull out items of interest from Stan output - could embed climate change simulation here

