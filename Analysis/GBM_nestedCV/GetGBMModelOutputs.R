home="/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/StatPigMvmt/Pipeline_R2"

##Identifying optimal parameters for gradient boosted models

#created using tutorials:
#https://bradleyboehmke.github.io/HOML/gbm.html
#https://rspatial.org/raster/sdm/9_sdm_brt.html

# Set variables ----------------

#vars that need to be changed from run to run
repname="10DEC24_Runs"
split_type="Region"

# Process Outline ----------------

#input: pigsums.csv

# Setup ----------------

#Set dirs
setwd(home)
#set path where gbm functions live
gbm_funcdir=file.path(home,"1_Scripts","Analysis","GBM_nestedCV","GBM_Functions", fsep = .Platform$file.sep)
#set objdir
objdir=file.path(home,"2_Data","Objects", fsep = .Platform$file.sep)

#Load libraries
require(tidyverse)
require(dplyr)
require(dismo)
require(gbm)
require(gpboost)
require(job)
require(SHAPforxgboost)
require(caret)
require(ggplot2)
require(glmnet)

#read data
pigsums=readRDS(file.path(objdir,"dailyPigSums.rds", fsep = .Platform$file.sep))

# Format Data ----------------

#set column used for region split
colnames(pigsums)[3]<-"Region"

#all covars need to be either numeric or factor
pigsums[,which(colnames(pigsums)=="season")]<-as.factor(pigsums[,which(colnames(pigsums)=="season")])
pigsums[,which(colnames(pigsums)=="period")]<-as.factor(pigsums[,which(colnames(pigsums)=="period")])
pigsums[,which(colnames(pigsums)=="sex")]<-as.factor(pigsums[,which(colnames(pigsums)=="sex")])

#response as integers for poisson
pigsums$sl_mean<-as.integer(pigsums$sl_mean)
pigsums$displ_mean<-as.integer(pigsums$displ_mean)

#check for any NAs
#number of NAs between sl/displ may differ
pigsums_sl<-pigsums[!is.na(pigsums$sl_mean),]
pigsums_displ<-pigsums[!is.na(pigsums$displ_mean),]

#make cutoff set for sigma models
pigsums_sigmasl=pigsums[pigsums$not_na_sl>30,]
pigsums_sigmadisp=pigsums[pigsums$not_na_displ>30,]

# Source functions
source(file.path(gbm_funcdir,"GBM_FunctionSourcer.R", fsep = .Platform$file.sep))

# Hyperparameter optimization and Variable selection ---------
repname="21FEB24_Runs"
split_types=c("Region","Random")
responses=colnames(pigsums)[59:62]
response_strings=c("sl","sigma.sl","disp","sigma.disp")
distributions=c("poisson","gaussian","poisson","gaussian")
pigsums_list=list(pigsums_sl,pigsums_sigmasl,pigsums_displ,pigsums_sigmadisp)

#Loop through each split type
for(s in 1:2){
split_type=split_types[s]

#Loop through each of four responses
for(r in 1:4){
response=responses[r]
response_str=response_strings[r]
distribution=distributions[r]
pigsums=pigsums_list[[r]]

  filestr=paste("4_Outputs",repname,split_type,sep="/")
  if(!dir.exists(file.path(home,"4_Outputs",repname, fsep = .Platform$file.sep))){dir.create(file.path(home,file.path("4_Outputs",repname, fsep = .Platform$file.sep)))}
  if(!dir.exists(file.path(home,filestr))){dir.create(file.path(home,filestr))}
  path=file.path(home,filestr,fsep = .Platform$file.sep)
MakeAllGBMOutputs(path,split_type,pigsums,response,response_str,distribution)

}
}