
# Purpose --------------------------------

# This script does predictions from top models for sl/disp responses for all watershed/seasons

#Inputs: wash, pigsums
#Output:

# Vars --------------------------------
reps=2
clustnum=4
filestr<-"04APR25_Runs"

# Setup --------------------------------

#load libraries
library(sf)
library(terra)
library(dplyr)
library(sp)
library(raster)
library(tidyr)
library(dismo)
library(gbm)
library(doParallel)
library(foreach)
library(abind)

## set working directories -------
#local:
home<-"/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/2_Projects/StatPigMvmt/Pipeline_R2"
#home="/home/kayleigh.chalkowski/RR_GBM"
outdir<-file.path(home,"4_Outputs")
objdir<-file.path(home,"2_Data","Objects")
indir<-file.path(home,"2_Data","Input","Env_Cov")

if(!dir.exists(file.path(outdir,filestr))){dir.create(file.path(outdir,filestr))}
if(!dir.exists(file.path(outdir,filestr,"FigTab"))){dir.create(file.path(outdir,filestr,"FigTab"))}
if(!dir.exists(file.path(outdir,filestr,"UncPredMats"))){dir.create(file.path(outdir,filestr,"UncPredMats"))}

## read data -------

#read in pigsums dataset
pigsums<-readRDS(file.path(objdir,"dailyPigSums.rds"))
pigswsite<-readRDS(file.path(objdir,"geolocsnatl_wDispl.rds"))

#Read in formatted wastershed shapefile
wash<-readRDS(file.path(indir,"wash3.rds"))

#Xsel and Xvec
Xsel<-readRDS(file.path(objdir,"X_sel.rds"))
X_vec_list<-readRDS(file.path(objdir,"X_vec_list.rds"))

#hardcoded to gbm 01 drop for now
sl.opt.params.kfold=read.csv(file.path(outdir,filestr,X_sel[1,3],"GBM_01Drop","sl_bestmodelparams.csv"))
sigma.sl.opt.params.kfold=read.csv(file.path(outdir,filestr,X_sel[2,3],"sigma.sl_bestmodelparams.csv"))
disp.opt.params.kfold=read.csv(file.path(outdir,filestr,X_sel[3,3],"disp_bestmodelparams.csv"))
sigma.disp.opt.params.kfold=read.csv(file.path(outdir,filestr,X_sel[4,3],"GBM_01Drop","sigma.disp_bestmodelparams.csv"))

#X_vec_list=getXvec(model.sel.tbl.total,"Xlist")

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
#usplot<-st_as_sf(usplot)
#usplot <- st_cast(usplot, "MULTIPOLYGON")
#usplot=st_transform(usplot,crs=st_crs(wash))

# Make functions to get uncertainty matrices ------------------

# Run predictions over reps for single period/quarter/response
GetPreds<-function(p,q,pigsums2,ww,reps,opt.params,X.vec,response,distribution){
  allq=grep("_q",colnames(ww))
  
  q_vars=grep(paste0("_q",q),colnames(ww))
  wash_notemporal=ww[,-allq]
  washq=cbind(wash_notemporal,ww[,q_vars])
  
  #fix names to match pigsums
  colnames(washq)=gsub("_q[0-9]","",colnames(washq))
  means=grep("_mn",colnames(washq))
  colnames(washq)=gsub("_mn","",colnames(washq))
  colnames(washq)[means]<-paste0("mean_",colnames(washq)[means])
  
  vars=grep("_var",colnames(washq))
  colnames(washq)=gsub("_var","",colnames(washq))
  colnames(washq)[vars]<-paste0("var_",colnames(washq)[vars])
  
  washq$season=as.factor(paste0("q",q))
  washq$period=p
  
  washqf=washq
  washqm=washq
  washqf$sex=as.factor("Female")
  washqm$sex=as.factor("Male")
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
    pred.m <- gbm::predict.gbm(gbm, washqm, n.trees=opt.params$n.trees, "response")
    
    #displacement per quarter female
    pred.f <- gbm::predict.gbm(gbm, washqf, n.trees=opt.params$n.trees, "response")
    
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

#Make function to run GetPreds in parallel over each quarter/period
RunUncertaintyParallel<-function(clustnum,pigsums2,
                                 ww,reps,opt.params,
                                 X.vec,response,
                                 distribution,path.out,
                                 filename){
  for(q in 1:4){
  print(paste0("q: ",q))
  for(p in 1:4){
  print(paste0("p: ",p))
  cl <- parallel::makeCluster(clustnum)
  doParallel::registerDoParallel(cl)
  plist=GetPreds(p,q,pigsums2,ww,reps,opt.params,X.vec,response,distribution)
  pmat=sapply(plist,as.matrix)
  parallel::stopCluster(cl)
  saveRDS(pmat,file.path(path.out,paste0(filename,"_",p,".rds")))
  }
  
  pmat_p1=readRDS(file.path(path.out,paste0(filename,"_",1,".rds")))
  pmat_p2=readRDS(file.path(path.out,paste0(filename,"_",2,".rds")))
  pmat_p3=readRDS(file.path(path.out,paste0(filename,"_",3,".rds")))
  pmat_p4=readRDS(file.path(path.out,paste0(filename,"_",4,".rds")))
  
  allpmat=abind(pmat_p1,pmat_p2,pmat_p3,pmat_p4,along=3)
  allpmat2=(rowSums(allpmat, dims = 2))/4
  
  #save compiled
  saveRDS(allpmat2,file.path(path.out,paste0(filename,"_q",q,".rds")))
  
  #delete temp files
  file.remove(file.path(path.out,paste0(filename,"_",1,".rds")))
  file.remove(file.path(path.out,paste0(filename,"_",2,".rds")))
  file.remove(file.path(path.out,paste0(filename,"_",3,".rds")))
  file.remove(file.path(path.out,paste0(filename,"_",4,".rds")))
  
  }
}

#Make function to average the quarters to get single uncertainty matrix
AverageQuarters<-function(pmat.path,pmat.str){
  pmat.fs=list.files(pmat.path,full.names=TRUE)
  pmat.fs=pmat.fs[grep(pmat.str,pmat.fs)]
  pmat_list <- lapply(pmat.fs, readRDS)
  
  #should only be list of four, else, stop
  if(length(pmat_list)>4){stop("incorrect pmats pulled for averaging")}
  
  pmat_1=pmat_list[[1]]
  pmat_2=pmat_list[[2]]
  pmat_3=pmat_list[[3]]
  pmat_4=pmat_list[[4]]
  
  allqmat=abind(pmat_1,pmat_2,pmat_3,pmat_4,along=3)
  allqmat2=(rowSums(allqmat, dims = 2))/4
  
  #Saving final pmat file
  print(paste0("Saving seasonally-averaged ",pmat.str))
  saveRDS(allqmat2,file.path(pmat.path,paste0(pmat.str,".rds")))
  
}

# Run function for each response and q ------------------
if(!dir.exists(file.path(outdir,filestr,"UncPredMats"))){dir.create(file.path(outdir,filestr,"UncPredMats"))}

path.out=file.path(outdir,filestr,"UncPredMats")

### SL ------------------

RunUncertaintyParallel(clustnum,
                       pigsums_sl,wash,reps,
                       sl.opt.params.kfold,
                       X_vec_list$Xsl,
                       "sl_mean","poisson",
                       path.out,"sl_pmat")

AverageQuarters(path.out,"sl_pmat")

### Sigma SL ------------------

RunUncertaintyParallel(clustnum,
                       pigsums_sigmasl,wash,reps,
                       sigma.sl.opt.params.kfold,
                       X_vec_list$sigmasl,
                       "sl_disp","gaussian",
                       path.out,"sigmasl_pmat")

AverageQuarters(path.out,"sigmasl_pmat")

### Disp ------------------

RunUncertaintyParallel(clustnum,
                       pigsums_disp,wash,reps,
                       disp.opt.params.kfold,
                       X_vec_list$sigmadisp,
                       "displ_mean","poisson",
                       path.out,"disp_pmat")

AverageQuarters(path.out,"disp_pmat")

### Sigma Disp ------------------

RunUncertaintyParallel(clustnum,
                       pigsums_sigmadisp,wash,reps,
                       sigma.disp.opt.params.kfold,
                       X_vec_list$Xdisp,
                       "displ_disp","gaussian",
                       path.out,"sigmadisp_pmat")

AverageQuarters(path.out,"sigmadisp_pmat")



