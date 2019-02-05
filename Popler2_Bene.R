
install.packages("devtools")
devtools::install_github("AldoCompagnoni/popler", build_vignettes=T, force=T)
install.packages("dplyr")
library(popler)
library(knitr)
library(dplyr)
library(devtools) #needed to download prism from github
library(reshape2) ##melting dataframes
library(raster) ##working with raster data
library(sp) ##manipulationg spatial data
install_github(repo = "prism", username = "ropensci")
library(prism) ##prism data access
library(segmented) #for peicewise regression
library(rstan)

datatype<-pplr_dictionary(datatype)[[1]]
#"individual"  "count"       "cover"       "biomass"     "density"     "basal_cover"
# count: data in abundance_observation model1
#biomass: data in abundance_observation model2
#cover: data in abundance_observation model3
#density: data in abundance_observation model??? weird densities model3
#individual: no basal_cover + observational
#basal_cover: no basal_cover + observational

obs_counts <- pplr_browse(datatype == datatype[2] & treatment == "observational")
str(obs_counts)
project<-unique(obs_counts$proj_metadata_key)

pass <- pplr_get_data(proj_metadata_key==project[1])
pass$abundance_observation<-round(as.numeric(pass$abundance_observation))  
pass$spatial_replication_level_1<-as.factor(pass$spatial_replication_level_1)  
pass$spatial_replication_level_2<-as.factor(pass$spatial_replication_level_2)
pass$spatial_replication_level_3<-as.factor(pass$spatial_replication_level_3)  
pass$spatial_replication_level_4<-as.factor(pass$spatial_replication_level_4)
rep<-paste(pass$spatial_replication_level_1,"_",pass$spatial_replication_level_2,"_",pass$spatial_replication_level_3,"_",pass$spatial_replication_level_4)
pass<-cbind(pass,rep)
sp<-unique(pass$species)
tt<-table(pass$species)
tt<-tt[tt>9]
sp<-dimnames(tt)[[1]]
sp<-sp[sp!="NA"]

LAMBDA<-c()
for(k in 1:length(sp))
{
  pass2<-pass[pass$species==sp[k],]
  pass3<-pass2[is.na(pass2$abundance_observation)==FALSE,]
  
  datalist<-list(
    n=nrow(pass3),
    nyear=length(unique(as.factor(pass3$year))),
    nrep=length(unique(as.factor(pass3$rep))),
    year=as.numeric(as.factor(as.character(pass3$year))),
    rep=as.numeric(as.factor(as.character(pass3$rep))),
    count=pass3$abundance_observation
  )
  data_fit<-stan(file="Model1.stan",data=datalist,iter=10000,chains=3)
  
  test2<-extract(data_fit,"a")[[1]]
  lamb<-exp(apply(test2,2,quantile,0.5))[2:length(apply(test2,2,quantile,0.5))]/exp(apply(test2,2,quantile,0.5))[1:(length(apply(test2,2,quantile,0.5))-1)]
  lambl<-exp(apply(test2,2,quantile,0.025))[2:length(apply(test2,2,quantile,0.025))]/exp(apply(test2,2,quantile,0.025))[1:(length(apply(test2,2,quantile,0.025))-1)]
  lambh<-exp(apply(test2,2,quantile,0.975))[2:length(apply(test2,2,quantile,0.975))]/exp(apply(test2,2,quantile,0.975))[1:(length(apply(test2,2,quantile,0.975))-1)]
  LAMBDA<-rbind(LAMBDA,cbind(rep(project[1],length(lamb)),rep(sp[k],length(lamb)),levels(as.factor(as.character(pass3$year)))[-1],lamb,lambl,lambh))
}
LAMBDA<-cbind(LAMBDA,rep(unique(obs_counts$lat_lter[obs_counts$proj_metadata_key==project[1]]),nrow(LAMBDA)),rep(unique(obs_counts$lng_lter[obs_counts$proj_metadata_key==project[1]]),nrow(LAMBDA)))
colnames(LAMBDA)<-c("Project","Species","Year","lamb","lambl","lambh","Lat","Lon")

for(i in 2:length(project))
{
  pass <- pplr_get_data(proj_metadata_key==project[i])
  pass$structure_type_1<-as.numeric(pass$structure_type_1)
  pass$spatial_replication_level_1<-as.factor(pass$spatial_replication_level_1)  
  pass$spatial_replication_level_2<-as.factor(pass$spatial_replication_level_2)
  pass$spatial_replication_level_3<-as.factor(pass$spatial_replication_level_3)  
  pass$spatial_replication_level_4<-as.factor(pass$spatial_replication_level_4)
  rep<-paste(pass$spatial_replication_level_1,"_",pass$spatial_replication_level_2,"_",pass$spatial_replication_level_3,"_",pass$spatial_replication_level_4)
  pass<-cbind(pass,rep)
  sp<-unique(pass$species)
  tt<-table(pass$species)
  tt<-tt[tt>9]
  sp<-dimnames(tt)[[1]]
  sp<-sp[sp!="NA"]
  
  for(k in 1:length(sp))
  {
    pass2<-pass[pass$species==sp[k],]
    pass3<-pass2[is.na(pass2$abundance_observation)==FALSE,]
    datalist<-list(
      n=nrow(pass3),
      nyear=length(unique(as.factor(pass3$year))),
      nrep=length(unique(as.factor(pass3$rep))),
      year=as.numeric(as.factor(as.character(pass3$year))),
      rep=as.numeric(as.factor(as.character(pass3$rep))),
      count=round(pass3$abundance_observation) #sometimes counts are slightly off
    )
    data_fit<-stan(file="Model1.stan",data=datalist,iter=10000,chains=3)
    
    test2<-extract(data_fit,"a")[[1]]
    lamb<-exp(apply(test2,2,quantile,0.5))[2:length(apply(test2,2,quantile,0.5))]/exp(apply(test2,2,quantile,0.5))[1:(length(apply(test2,2,quantile,0.5))-1)]
    lambl<-exp(apply(test2,2,quantile,0.025))[2:length(apply(test2,2,quantile,0.025))]/exp(apply(test2,2,quantile,0.025))[1:(length(apply(test2,2,quantile,0.025))-1)]
    lambh<-exp(apply(test2,2,quantile,0.975))[2:length(apply(test2,2,quantile,0.975))]/exp(apply(test2,2,quantile,0.975))[1:(length(apply(test2,2,quantile,0.975))-1)]
    LAMBDA<-rbind(LAMBDA,cbind(rep(project[i],length(lamb)),rep(sp[k],length(lamb)),levels(as.factor(as.character(pass3$year)))[-1],lamb,lambl,lambh,rep(unique(obs_counts$lat_lter[obs_counts$proj_metadata_key==project[i]]),length(lamb)),rep(unique(obs_counts$lng_lter[obs_counts$proj_metadata_key==project[i]]),length(lamb))))
  }
}

LAMBDA<-data.frame(LAMBDA)
LAMBDA[,1]<-as.numeric(as.character(LAMBDA[,1]))
LAMBDA[,2]<-as.character(LAMBDA[,2])
LAMBDA[,3]<-as.numeric(as.character(LAMBDA[,3]))
LAMBDA[,4]<-as.numeric(as.character(LAMBDA[,4]))
LAMBDA[,5]<-as.numeric(as.character(LAMBDA[,5]))
LAMBDA[,6]<-as.numeric(as.character(LAMBDA[,6]))
LAMBDA[,7]<-as.numeric(as.character(LAMBDA[,7]))
LAMBDA[,8]<-as.numeric(as.character(LAMBDA[,8]))



LAMBDA_saved<-LAMBDA
LAMBDA<-LAMBDA[is.na(LAMBDA$Year)==FALSE,]
LAMBDA$Lat<-round(LAMBDA$Lat,1)
LAMBDA$Lon<-round(LAMBDA$Lon,1)
options(prism.path = "/Volumes/collnell/CAT/data/prism/prism_temp")
get_prism_annual(type = 'ppt', years=1990:2016, keepZip = TRUE)
ls_prism_data(name=TRUE)
year<-c(1990:2016)
ppt<-c()
for(i in 1:nrow(LAMBDA))
{
  new_file<-c(1:27)[year==LAMBDA$Year[i]] ##change to corresponding file numbers
  RS <- prism_stack(ls_prism_data()[new_file,1]) ##raster file
  df <- data.frame(rasterToPoints(RS)) ##creates a dataframe of points
  annual.df <- melt(df, c("x", "y"))
  names(annual.df)[1:2] <- c("lon", "lat") #rename columns
  annual.df$lat<-round(annual.df$lat,1)
  annual.df$lon<-round(annual.df$lon,1)
  ppt<-c(ppt,mean(annual.df$value[annual.df$lat==LAMBDA$Lat[i] & annual.df$lon==LAMBDA$Lon[i]]))
}
LAMBDA<-cbind(LAMBDA,ppt)

get_prism_annual(type = 'tmean', years=1990:2016, keepZip = TRUE)
ls_prism_data(name=TRUE)
year<-c(1990:2016)
tmean<-c()
for(i in 1:nrow(LAMBDA))
{
  new_file<-c(29:55)[year==LAMBDA$Year[i]] ##change to corresponding file numbers
  RS <- prism_stack(ls_prism_data()[new_file,1]) ##raster file
  df <- data.frame(rasterToPoints(RS)) ##creates a dataframe of points
  annual.df <- melt(df, c("x", "y"))
  names(annual.df)[1:2] <- c("lon", "lat") #rename columns
  annual.df$lat<-round(annual.df$lat,1)
  annual.df$lon<-round(annual.df$lon,1)
  tmean<-c(tmean,mean(annual.df$value[annual.df$lat==LAMBDA$Lat[i] & annual.df$lon==LAMBDA$Lon[i]]))
}
LAMBDA<-cbind(LAMBDA,tmean)


# Need the stan code for spline


