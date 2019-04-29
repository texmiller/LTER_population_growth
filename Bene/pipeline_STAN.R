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
library(rstan)
library(SPEI)
library(mgcv)
source("API_workaround.R")

## Here are all the studies we will use (observational only, excluding individual and basal cover)
obs_studies <- pplr_browse(studytype=="obs" & datatype!="individual" & datatype!="basal_cover",full_tbl = T)
#write_csv(obs_studies %>% 
#            select(proj_metadata_key,title,metalink),"obs_studies.csv")
## we think we want to drop study 300 (end year = 2030)

## proj_metadata_key for all data types
obs_count <- pplr_browse(studytype=="obs" & datatype == "count")$proj_metadata_key
obs_cover <- pplr_browse(studytype=="obs" & datatype == "cover")$proj_metadata_key
obs_density <- pplr_browse(studytype=="obs" & datatype == "density")$proj_metadata_key
obs_biomass <- pplr_browse(studytype=="obs" & datatype == "biomass")$proj_metadata_key

## load climate data, already subsetted for lat/longs of the popler studies (see Bene's script)
prism <- read_csv("prismdata.csv") 
## census month info
census_months <- read_csv("census_months.csv")

## Specify the proj_metadata_key of the focal data set 
## let's use the heron data set (88) as a guinea pig
k <- 88

  ## extract popler project data and metadata
  command <- paste0('pplr_browse( proj_metadata_key==',k,', full_tbl = T)')
  metadat<-parse(n=1, text=command) %>% eval
  ## diagnose the data type
  type <- metadat$datatype
  ## break out of function if datatype is individual or basal cover
  ##if(type=="individual" | type=="basal_cover"){return("Non-desired data type")}
  
  ## get data and combine spatial rep info
  n_spat_levels <- metadat$n_spat_levs
  #dat <- pplr_get_data(metadat,cov_unpack = T) %>% ## API problems, use workaround
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
  if("sppcode"%in%colnames(dat)==FALSE) {colnames(dat)[colnames(dat)=="species"]<-"sppcode"}
  
  ## find species that were never observed in one or more years, and drop
  drop_rare <- dat %>% 
    group_by(sppcode,year) %>% 
    summarise(total_obs_per_year = sum(abundance_observation))%>% 
    filter(total_obs_per_year == 0)%>% 
    summarise(unique(sppcode))
  ## drop rare spp and prep data for analysis
  newdat <- dat[(dat$sppcode%in%drop_rare$sppcode)==FALSE,] 
  
  ## if abudance_obs is structured, sum over structure
  ## NEED TO DO -- check with Aldo
  
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
  
  temp<-newdat %>% 
    group_by(sppcode,year) %>%  
    summarise(n())
  
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
                       ifelse(type=="density","Bene\\biomass_density_model.stan",
                              ifelse(type=="biomass","Bene\\biomass_density_model.stan",
                                     "Bene\\cover_model.stan")))

## Need to round up count and cover because some entries have decimals  
if(type=="count" | type=="cover"){datalist$count<-round(datalist$count)}

  abund_fit<-stan(file=stan_model,data=datalist,iter=5000,chains=3,warmup=500)
  
  ## Bayesian p-value
  newy<-rstan::extract(abund_fit,"newy")[[1]]
  new_mean<-apply(newy,2,mean)
  bayes.p<-length(new_mean[new_mean>mean(newdat$abundance_observation)])/length(new_mean)
  
  ## climate covariate
  study_site <- metadat$lterid
  site_lat <- round(metadat$lat_lter,2)
  site_long <- round(metadat$lng_lter,2)
  study_years <- sort(unique(newdat$year)) #metadat$studystartyr:metadat$studyendyr
  census_month <- 5 #read_csv("census_months.csv") %>% 
  #filter(proj_metadata_key == k) %>% 
  #summary(unique(Census_month))
  
  climate <- prism %>% 
    dplyr::select(-X1) %>% 
    filter(LTERlocation.id == study_site,
           year %in% min(study_years):max(study_years)) %>% 
    group_by_at(vars(-value)) %>%  # group by everything other than the value column. 
    mutate(row_id=1:n()) %>% ungroup() %>%  # build group index
    spread(key = variable, value = value)%>% 
    dplyr::select(-row_id) %>%
    ##thornthwaite function does not like NAs
    filter(!is.na(tmean),
           !is.na(ppt))%>% 
    distinct() %>% 
    ##df with years and months for ppt and temp
    ## convert calendar year to transition year
    mutate(climate_year = ifelse(month >=census_month, year+1, year),
           # Compute potential evapotranspiration (PET) and climatic water balance (BAL)
           PET = thornthwaite(tmean,site_lat),
           BAL = ppt-PET) %>% 
    filter(climate_year>min(study_years) & climate_year<=max(study_years)) 
  ## climate data should start 12 months before first abundance observation
  spei12 <- spei(climate[,'BAL'], 12)
  
  ## climate data for prediction
  climate_pred <- prism %>% 
    dplyr::select(-X1) %>% 
    filter(LTERlocation.id == study_site) %>% 
    group_by_at(vars(-value)) %>%  # group by everything other than the value column. 
    mutate(row_id=1:n()) %>% ungroup() %>%  # build group index
    spread(key = variable, value = value)%>% 
    dplyr::select(-row_id) %>%
    ##thornthwaite function does not like NAs
    filter(!is.na(tmean),
           !is.na(ppt))%>% 
    distinct() %>% 
    ##df with years and months for ppt and temp
    ## convert calendar year to transition year
    mutate(climate_year = ifelse(month >=census_month, year+1, year),
           # Compute potential evapotranspiration (PET) and climatic water balance (BAL)
           PET = thornthwaite(tmean,site_lat),
           BAL = ppt-PET) %>% 
    filter(climate_year<=max(year))
    climate_pred <- climate_pred[-(1:(census_month-1)),]
    spei12_pred <- spei(climate_pred[,'BAL'], 12)
    
  ## pull out relevant quantities from Stan output

  lambda <- rstan::extract(abund_fit,"lambda")[[1]][,(study_years[2:length(study_years)]-study_years[1:(length(study_years)-1)])==1,]
  #atest<-a[,(study_years[2:length(study_years)]-study_years[1:(length(study_years)-1)])==1,]
  year<-as.numeric(levels(as.factor(as.character(newdat$year))))[-1][(study_years[2:length(study_years)]-study_years[1:(length(study_years)-1)])==1]
  sp<-rep(levels(as.factor(as.character(newdat$sppcode))),(datalist$nyear-1))
  
  for(i in 1:abund_fit@sim$iter){
    lambda_i <- lambda[i,,]
    df_i <- data.frame(lambda_i)
    colnames(df_i) <- levels(as.factor(as.character(newdat$sppcode)))
    df_i$year <- year
    df_i <- df_i %>% 
      gather(BC:TH,key="species",value="lambda") %>% 
      mutate(r = log(lambda),
             species=as.factor(species))
    
    lambda_clim_i <- full_join(df_i,
                        tibble(SPEI = spei12$fitted[seq(from=12,to=length(spei12$fitted),by=12)],
                               year = (min(study_years):max(study_years))[-1]),
                        by="year")%>% 
      filter(!is.na(lambda))
    
    gam_fit_i <- gam(r ~ species + s(SPEI, by=species), data=lambda_clim_i)
    
    gam_dat<-plot(gam_fit_i)
    
    plot(gam_dat[[1]]$x)
  }
  
  
  set.seed(0)
  dat <- gamSim(1, n=400, dist="normal", scale=2) 
  b   <- gam(y~s(x0), data=dat) 
  plot(b, seWithMean=TRUE)
  
  
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
  
  ## Histogram of R-hats
  Rhat.plot <- ggplot(r.out)+
    geom_histogram(aes(x=r_Rhats))
    
  r.plot <- ggplot(r.out)+
    geom_point(aes(x=as.factor(year),y=r_mean))+
    geom_errorbar(aes(x=as.factor(year),ymin=r_lowCI, ymax=r_highCI))+
    facet_wrap(~sp)+
    ggtitle(metadat$title)+ 
    theme(axis.text.x = element_text(angle = 270, hjust = 1))
  


 
  
  
  
  r_clim %>% 
    ggplot()+
    geom_point(aes(x=SPEI,y=r_mean))+
    facet_grid(~sp)
  
  test <- r_clim %>% filter(sp=="CE")
  test_gam <- gam(r_mean ~ s(SPEI), data=test)
  plot(test_gam)
  points(test$SPEI,test$r_mean)
  summary(test_gam)
  
  predict.gam(test_gam,newdata=data.frame(seq(-1.5,1.5,0.1)))

  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
pipeline <- function(k){  
  return(list(metadat=metadat,
              r_clim=r_clim,
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





