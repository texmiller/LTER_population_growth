## Authors: LTER pipeline team
## Purpose: build a general function to analyze popler data sets and generate desired outputs
## Last update: February 26, 2019

## install popler and packages
install.packages("devtools")
devtools::install_github("AldoCompagnoni/popler", build_vignettes=T, force=T)
install.packages("dplyr")
library(popler)
library(knitr)
library(dplyr)
library(tidyverse)
library(devtools) #needed to download prism from github
library(reshape2) ##melting dataframes
library(raster) ##working with raster data
library(sp) ##manipulationg spatial data
#install_github(repo = "prism", username = "ropensci")
library(prism) ##prism data access
library(segmented) #for peicewise regression
library(rstan)
source("API_workaround.R")

## let's use the heron data set (88) as a guinea pig
#k <- 88
obs_studies <- pplr_browse(studytype=="obs" & datatype!="individual" & datatype!="basal_cover",full_tbl = T)
#write_csv(obs_studies %>% 
#            select(proj_metadata_key,title,metalink),"obs_studies.csv")

## we think we want to drop study 300 (end year = 2030)

## all of our data types
obs_count <- pplr_browse(studytype=="obs" & datatype == "count")$proj_metadata_key
obs_cover <- pplr_browse(studytype=="obs" & datatype == "cover")$proj_metadata_key
obs_density <- pplr_browse(studytype=="obs" & datatype == "density")$proj_metadata_key
obs_biomass <- pplr_browse(studytype=="obs" & datatype == "biomass")$proj_metadata_key

## climate data, already subsetted for lat/longs of the popler studies (see Bene's script)
prism <- read_csv("prismdata.csv") 


pipeline <- function(k){

  ## extract popler project data and metadata
  command <- paste0('pplr_browse( proj_metadata_key==',k,', full_tbl = T)')
  metadat<-parse(n=1, text=command) %>% eval
  ## diagnose the data type
  type <- metadat$datatype
  ## break out of function if datatype is individual or basal cover
  if(type=="individual" | type=="basal_cover"){return("Non-desired data type")}
  
  ## get data and combine spatial rep info
  n_spat_levels <- metadat$n_spat_levs
  #dat <- pplr_get_data(metadat,cov_unpack = T) %>% ## API problems
  dat <- query_get(conn, efficienty_query( k )) %>%   
    as.data.frame %>% 
    mutate(n_spat_levels = n_spat_levels) %>% 
    mutate(ran_effect = ifelse(n_spat_levels==1,spatial_replication_level_1,
                               ifelse(n_spat_levels==2,interaction(spatial_replication_level_1,spatial_replication_level_2),
                                      ifelse(n_spat_levels==3,interaction(spatial_replication_level_1,spatial_replication_level_2,spatial_replication_level_3),
                                             interaction(spatial_replication_level_1,spatial_replication_level_2,spatial_replication_level_3,spatial_replication_level_4))))) %>% 
    mutate(abundance_observation = count_observation) %>%  # Aldo's current code in the API bypass calls all abundance obs "count_observation"
    filter(!is.na(abundance_observation))
  ## filter out NAs and very rare species -- we made a decision to use only the data provided by PIs-- we are not
  ## assumming that NAs are zero
  
  ## sometimes sppcode does not exist, so we have to use species
  if("sppcode"%in%colnames(dat)==FALSE)
  {
    colnames(dat)[colnames(dat)=="species"]<-"sppcode"
  }
  
  ## find species that were never observed in one or more years, sometimes sppcode does not exist, so we have to use species
  drop_rare <- dat %>% 
    group_by(sppcode,year) %>% 
    summarise(total_obs_per_year = sum(abundance_observation))%>% 
    filter(total_obs_per_year == 0)%>% 
    summarise(unique(sppcode))
  ## drop rare spp and prep data for analysis
  newdat <- dat[(dat$sppcode%in%drop_rare$sppcode)==FALSE,] 
  
  ## if abudance_obs is structured, sum over structure
  ## check with Aldo
  
  #cover ranges from 0-100 and sometimes ranges from 0-1.. with a few entries above 1. I still have to go through all of the cover dataset (26ish) to see what else is there.
  if(type=="cover")
  {
    q<-quantile(newdat$abundance_observation)
    if(q[4]<1)
    {
      newdat$abundance_observation<-newdat$abundance_observation*100
      newdat<-newdat[newdat$abundance_observation<=100,]
    }
  }
  
  #prep data for stan
  datalist<-list(
    n=nrow(newdat),
    nyear=length(unique(as.factor(newdat$year))),
    nrep=length(unique(as.factor(newdat$ran_effect))),
    nsp=length(unique(as.factor(newdat$sppcode))),
    year=as.numeric(as.factor(as.character(newdat$year))),
    rep=as.numeric(as.factor(as.character(newdat$ran_effect))),
    sp=as.numeric(as.factor(as.character(newdat$sppcode))),
    count=newdat$abundance_observation)
  
  ## send to the appropriate JAGS model, given data type
  stan_model <- ifelse(type=="count","Bene\\count_model.stan",
                       ifelse(type=="density","biomass_density_model.stan",
                              ifelse(type=="biomass","biomass_density_model.stan","cover_model.stan")))

## Need to round up count and cover because some entries have decimals  
if(type=="count" | type=="cover")
{
  datalist$count<-round(datalist$count)
}

  abund_fit<-stan(file=stan_model,data=datalist,iter=5000,chains=3,warmup=500)
  
  ## Bayesian p-value
  newy<-rstan::extract(abund_fit,"newy")[[1]]
  new_mean<-apply(newy,2,mean)
  bayes.p<-length(new_mean[new_mean>mean(newdat$abundance_observation)])/length(new_mean)

  ## collect posterior mean and CI growth rates, convert to year * species matrix
  r.out<-summary(abund_fit,"r")[[1]][,c(1,4,8,10)]
  year<-rep(as.numeric(levels(as.factor(as.character(newdat$year))))[-1],1,each=datalist$nsp)
  sp<-rep(levels(as.factor(as.character(newdat$sppcode))),(datalist$nyear-1))
  r.out<-cbind(r.out,year,sp)
  r.out<-data.frame(r.out)
  for(j in 1:5)
  {
    r.out[,j]<-as.numeric(as.character(r.out[,j]))
  }
  colnames(r.out)<-c("r_mean","r_lowCI","r_highCI","r_Rhats","year","sp")
  
  r.plot <- ggplot(r.out)+
    geom_point(aes(x=as.factor(year),y=r_mean))+
    geom_errorbar(aes(x=as.factor(year),ymin=r_lowCI, ymax=r_highCI))+
    facet_wrap(~sp)+
    ggtitle(metadat$title)
  
  ## climate covariate
  study_years <- metadat$studystartyr:metadat$studyendyr
  study_site <- metadat$lterid
  site_lat <- round(metadat$lat_lter,2)
  site_long <- round(metadat$lng_lter,2)
  
  census_month <- 5 #read_csv("census_months.csv") %>% 
    #filter(proj_metadata_key == k) %>% 
    #summary(unique(Census_month))
  
  climate <- prism %>% 
    filter(lter == study_site,
           #lat == site_lat,
           #lon == site_long,
           year %in% study_years) %>% 
    spread(key = variable, value = value)
  
  %>% 
    ##df with years and months for ppt and temp
    ## convert calendar year to transition year
    mutate(climate_year = ifelse(month >=census_month, year+1, year),
           # Compute potential evapotranspiration (PET) and climatic water balance (BAL)
           PET = thornthwaite(tmean,site_lat),
           BAL = ppt-PET) %>% 
    filter(climate_year>=min(study_years))
  
  ## climate data should start 12 months before first abundance observation
  spei12 <- spei(climate[,'BAL'], 12)
  
  return(list(metadat=metadat,
              r.out=r.out,
              r.plot=r.plot,
              bayes.p=bayes.p))
}



pplr_dictionary(full_tbl = T)

pplr_dictionary(samplefreq)

plot(as.numeric(obs_studies$samplefreq))

thing <- pipeline(k=88) 

## plot r estimates
thing$r.plot

## look at distribution of R-hat values
hist(thing$r.out$r_Rhats)

## Bayesian p-value
thing$bayes.p




