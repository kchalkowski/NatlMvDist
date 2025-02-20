
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

## read data -------

#read in pigsums dataset
pigsums<-readRDS(file.path(objdir,"dailyPigSums.rds"))
pigswsite<-readRDS(file.path(objdir,"geolocsnatl_wDispl.rds"))

#Read in formatted wastershed shapefile
wash<-readRDS(file.path(objdir,"wash_envcov_final.rds"))

#Read in X sel table/X vec list
X_sel=readRDS(file.path(objdir,"X_sel.rds"))
X_vec_list=readRDS(file.path(objdir,"X_vec_list.rds"))

#read in optimal hyperparameter sets
sl.opt.params.kfold=read.csv(file.path(outdir,filestr,X_sel[1,3],"sl_bestmodelparams.csv"))
sigma.sl.opt.params.kfold=read.csv(file.path(outdir,filestr,X_sel[2,3],"sigma.sl_bestmodelparams.csv"))
disp.opt.params.kfold=read.csv(file.path(outdir,filestr,X_sel[3,3],"disp_bestmodelparams.csv"))
sigma.disp.opt.params.kfold=read.csv(file.path(outdir,filestr,X_sel[4,3],"sigma.disp_bestmodelparams.csv"))

#get helper shapefiles
usplot <- st_read(dsn = file.path(objdir,"usmapplot_best.shp"), layer = "usmapplot_best")

## format data -------

#format pigsums data
#set column used for region split
colnames(pigsums)[1]<-"Region"

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

#format helper shapefiles for plotting
usplot<-st_as_sf(usplot)
usplot <- st_cast(usplot, "MULTIPOLYGON")
usplot=st_transform(usplot,crs=st_crs(wash))

## Get Preds function ------------------
GetPreds<-function(pigsums2,washq,reps,opt.params,X.vec,response,distribution){
  #1-initialize matrix
  predmat=matrix(nrow=nrow(washq),ncol=reps)
  #2-loop
  #for(i in 1:reps){
  foreach(i=1:reps) %dopar% {
    print(i)
    gbm=dismo::gbm.fixed(data=pigsums2, gbm.x=X.vec, gbm.y=which(colnames(pigsums2)==response),
                         learning.rate=opt.params$learning.rate, 
                         tree.complexity=opt.params$tree.complexity, 
                         n.trees=opt.params$n.trees,
                         bag.fraction=opt.params$bag.fraction,
                         family=distribution) 
    #3-take m/f for each prediction
    #displacement per quarter male
    pred.m <- gbm::predict.gbm(gbm, washq, n.trees=opt.params$n.trees, "response")
    
    #displacement per quarter female
    pred.f <- gbm::predict.gbm(gbm, washq, n.trees=opt.params$n.trees, "response")
    
    #combine quarterly estimates into matrix to get means
    pred=rowMeans(cbind(pred.m,pred.f))
    
    pred
    #5-put average in matrix
    #predmat[,i]=pred  
    #pred
    
  }
  #6-output matrix
  #return(predmat)
  #return(predmat)
  
}

## Run function ------------------
### Washq1 ------------------

cl <- parallel::makeCluster(4)
doParallel::registerDoParallel(cl)
plist=GetPreds(pigsums2,washq1,1000,sl.opt.params.kfold_region,X_vec.start,"sl_","poisson")
pmat.q1.sl=sapply(plist,as.matrix)
parallel::stopCluster(cl)
write.csv(pmat.q1.sl,paste0(outdir,"/PredMats_22APR23/pmat_q1_sl.csv"))

cl <- parallel::makeCluster(4)
doParallel::registerDoParallel(cl)
plist=GetPreds(cutoff_sl,washq1,1000,sigma.sl.opt.params.kfold,X_vec.start,"sigma_sl","gaussian")
pmat.q1.sigmasl=sapply(plist,as.matrix)
parallel::stopCluster(cl)
write.csv(pmat.q1.sigmasl,paste0(outdir,"/PredMats_22APR23/pmat_q1_sigmasl.csv"))


cl <- parallel::makeCluster(4)
doParallel::registerDoParallel(cl)
plist=GetPreds(pigsums2,washq1,1000,disp.opt.params.kfold,X.vec.disp.01drop,"displacement","poisson")
pmat.q1.disp=sapply(plist,as.matrix)
parallel::stopCluster(cl)
write.csv(pmat.q1.disp,paste0(outdir,"/PredMats_22APR23/pmat_q1_disp.csv"))

cl <- parallel::makeCluster(4)
doParallel::registerDoParallel(cl)
plist=GetPreds(cutoff_sl,washq1,1000,sigma.disp.opt.params.kfold_region,X_vec.start,"sigma_disp","gaussian")
pmat.q1.sigmadisp=sapply(plist,as.matrix)
parallel::stopCluster(cl)
write.csv(pmat.q1.sigmadisp,paste0(outdir,"/PredMats_22APR23/pmat_q1_sigmadisp.csv"))

cl <- parallel::makeCluster(4)
doParallel::registerDoParallel(cl)
plist=GetPreds(pigsums2,washq1,1000,tenavg.opt.params.kfold_region,X_vec.start,"tenavg","poisson")
pmat.q1.tenavg=sapply(plist,as.matrix)
parallel::stopCluster(cl)
write.csv(pmat.q1.tenavg,paste0(outdir,"/PredMats_22APR23/pmat_q1_tenavg.csv"))

### Washq2 ------------------

cl <- parallel::makeCluster(4)
doParallel::registerDoParallel(cl)
plist=GetPreds(pigsums2,washq2,1000,sl.opt.params.kfold_region,X_vec.start,"sl_","poisson")
pmat.q2.sl=sapply(plist,as.matrix)
parallel::stopCluster(cl)
write.csv(pmat.q2.sl,paste0(outdir,"/PredMats_22APR23/pmat_q2_sl.csv"))

cl <- parallel::makeCluster(4)
doParallel::registerDoParallel(cl)
plist=GetPreds(cutoff_sl,washq2,1000,sigma.sl.opt.params.kfold,X_vec.start,"sigma_sl","gaussian")
pmat.q2.sigmasl=sapply(plist,as.matrix)
parallel::stopCluster(cl)
write.csv(pmat.q2.sigmasl,paste0(outdir,"/PredMats_22APR23/pmat_q2_sigmasl.csv"))

cl <- parallel::makeCluster(4)
doParallel::registerDoParallel(cl)
plist=GetPreds(pigsums2,washq2,1000,disp.opt.params.kfold,X.vec.disp.01drop,"displacement","poisson")
pmat.q2.disp=sapply(plist,as.matrix)
parallel::stopCluster(cl)
write.csv(pmat.q2.disp,paste0(outdir,"/PredMats_22APR23/pmat_q2_disp.csv"))

cl <- parallel::makeCluster(4)
doParallel::registerDoParallel(cl)
plist=GetPreds(cutoff_sl,washq2,1000,sigma.disp.opt.params.kfold_region,X_vec.start,"sigma_disp","gaussian")
pmat.q2.sigmadisp=sapply(plist,as.matrix)
parallel::stopCluster(cl)
write.csv(pmat.q2.sigmadisp,paste0(outdir,"/PredMats_22APR23/pmat_q2_sigmadisp.csv"))

### Washq3 ------------------

cl <- parallel::makeCluster(4)
doParallel::registerDoParallel(cl)
plist=GetPreds(pigsums2,washq3,1000,sl.opt.params.kfold_region,X_vec.start,"sl_","poisson")
pmat.q3.sl=sapply(plist,as.matrix)
parallel::stopCluster(cl)
write.csv(pmat.q3.sl,paste0(outdir,"/PredMats_22APR23/pmat_q3_sl.csv"))

cl <- parallel::makeCluster(4)
doParallel::registerDoParallel(cl)
plist=GetPreds(cutoff_sl,washq3,1000,sigma.sl.opt.params.kfold,X_vec.start,"sigma_sl","gaussian")
pmat.q3.sigmasl=sapply(plist,as.matrix)
parallel::stopCluster(cl)
write.csv(pmat.q3.sigmasl,paste0(outdir,"/PredMats_22APR23/pmat_q3_sigmasl.csv"))

cl <- parallel::makeCluster(4)
doParallel::registerDoParallel(cl)
plist=GetPreds(pigsums2,washq3,1000,disp.opt.params.kfold,X.vec.disp.01drop,"displacement","poisson")
pmat.q3.disp=sapply(plist,as.matrix)
parallel::stopCluster(cl)
write.csv(pmat.q3.disp,paste0(outdir,"/PredMats_22APR23/pmat_q3_disp.csv"))

cl <- parallel::makeCluster(4)
doParallel::registerDoParallel(cl)
plist=GetPreds(cutoff_sl,washq3,1000,sigma.disp.opt.params.kfold_region,X_vec.start,"sigma_disp","gaussian")
pmat.q3.sigmadisp=sapply(plist,as.matrix)
parallel::stopCluster(cl)
write.csv(pmat.q3.sigmadisp,paste0(outdir,"/PredMats_22APR23/pmat_q3_sigmadisp.csv"))

### Washq4 ------------------

cl <- parallel::makeCluster(4)
doParallel::registerDoParallel(cl)
plist=GetPreds(pigsums2,washq4,1000,sl.opt.params.kfold_region,X_vec.start,"sl_","poisson")
pmat.q4.sl=sapply(plist,as.matrix)
parallel::stopCluster(cl)
write.csv(pmat.q4.sl,paste0(outdir,"/PredMats_22APR23/pmat_q4_sl.csv"))

cl <- parallel::makeCluster(4)
doParallel::registerDoParallel(cl)
plist=GetPreds(cutoff_sl,washq4,1000,sigma.sl.opt.params.kfold,X_vec.start,"sigma_sl","gaussian")
pmat.q4.sigmasl=sapply(plist,as.matrix)
parallel::stopCluster(cl)
write.csv(pmat.q4.sigmasl,paste0(outdir,"/PredMats_22APR23/pmat_q4_sigmasl.csv"))

cl <- parallel::makeCluster(4)
doParallel::registerDoParallel(cl)
plist=GetPreds(pigsums2,washq4,1000,disp.opt.params.kfold,X.vec.disp.01drop,"displacement","poisson")
pmat.q4.disp=sapply(plist,as.matrix)
parallel::stopCluster(cl)
write.csv(pmat.q4.disp,paste0(outdir,"/PredMats_22APR23/pmat_q4_disp.csv"))

cl <- parallel::makeCluster(4)
doParallel::registerDoParallel(cl)
plist=GetPreds(cutoff_sl,washq4,1000,sigma.disp.opt.params.kfold_region,X_vec.start,"sigma_disp","gaussian")
pmat.q4.sigmadisp=sapply(plist,as.matrix)
parallel::stopCluster(cl)
write.csv(pmat.q4.sigmadisp,paste0(outdir,"/PredMats_22APR23/pmat_q4_sigmadisp.csv"))

