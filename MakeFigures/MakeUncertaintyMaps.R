
# Purpose --------------------------------

# This script does predictions from top models for sl/disp responses for all watershed/seasons

#Inputs: wash_envcov_final.rds, pigsums, XXX
#Output:

# Setup --------------------------------

#load libraries
library(job)
library(tictoc)
library(sf)
library(terra)
library(exactextractr)
library(dplyr)
library(sp)
library(raster)
library(tidyr)
library(dismo)
library(gbm)
library(doParallel)
library(foreach)

## set working directories -------
#local:
home<-"/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/StatPigMvmt/Pipeline_R2"
filestr<-"TestRuns"
outdir<-file.path(home,"4_Outputs",filestr)
objdir<-file.path(home,"2_Data","Objects")
pmat_path=file.path(objdir,"UncPredMats")

## grab uncertainty matrices -------

list.files(path=pmat_path,full.names=TRUE)





