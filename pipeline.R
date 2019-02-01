## Authors: LTER pipeline team
## Purpose: build a general function to analyze popler data sets and generate desired outputs
## Last update: February 1, 2019

## install popler
#install.packages("devtools")
#devtools::install_github("AldoCompagnoni/popler", build_vignettes=T, force=T)
#install.packages("dplyr")
library(popler)

bigfun <- function(popler_obj){
  
  ## extract popler project metadata
  
  ## collect climate covariates
  
  ## diagnose the data type
  
  ## prep data for analysis
  
  ## send to the appropriate Stan model, given data type
  
  ## pull out items of interest from Stan output - could embed climate change simulation here
  
  ## package outputs into data frame or list
  
}

## pseudo-code for Stan model
#1. Get year-specific lambdas
#2. Fit spline for climate covariate
#3. Derive expected value
#4. Calculate model diagnostics
