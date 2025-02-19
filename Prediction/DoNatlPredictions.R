
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

## format data -------

#format pigsums data

#format watersheds shapefiles
#wash2<-st_as_sf(wash)
#wash2 <- st_cast(wash2, "POLYGON")

#Separate wash q's into m/f for separate analysis
washf=wash
washm=wash
washf$sex=as.factor("Female")
washm$sex=as.factor("Male")

# Run GBM models ----------------------------
gbm.sl=gbm.fixed(data=pigsums_sl, gbm.x=X_vec_list$Xsl, gbm.y=which(colnames(pigsums_sl)=="sl_mean"),
                 learning.rate=sl.opt.params.kfold$learning.rate, 
                 tree.complexity=sl.opt.params.kfold$tree.complexity, 
                 n.trees=sl.opt.params.kfold$n.trees,
                 bag.fraction=sl.opt.params.kfold$bag.fraction,
                 family="poisson") 
gbm.disp=gbm.fixed(data=pigsums_displ, gbm.x=X_vec_list$Xdisp, gbm.y=which(colnames(pigsums)=="displ_mean"),
                   learning.rate=disp.opt.params.kfold$learning.rate, 
                   tree.complexity=disp.opt.params.kfold$tree.complexity, 
                   n.trees=disp.opt.params.kfold$n.trees,
                   bag.fraction=disp.opt.params.kfold$bag.fraction,
                   family="poisson") 
gbm.sigma.sl=gbm.fixed(data=pigsums_sigmasl, gbm.x=X_vec_list$sigmasl, gbm.y=which(colnames(pigsums)=="sl_disp"),
                       learning.rate=sigma.sl.opt.params.kfold$learning.rate, 
                       tree.complexity=sigma.sl.opt.params.kfold$tree.complexity, 
                       n.trees=sigma.sl.opt.params.kfold$n.trees,
                       bag.fraction=sigma.sl.opt.params.kfold$bag.fraction,
                       family="gaussian") 
gbm.sigma.disp=gbm.fixed(data=pigsums_sigmadisp, gbm.x=X_vec_list$sigmadisp, gbm.y=which(colnames(pigsums)=="displ_disp"),
                         learning.rate=sigma.disp.opt.params.kfold$learning.rate, 
                         tree.complexity=sigma.disp.opt.params.kfold$tree.complexity, 
                         n.trees=sigma.disp.opt.params.kfold$n.trees,
                         bag.fraction=sigma.disp.opt.params.kfold$bag.fraction,
                         family="gaussian") 

# Predictions for male/female -----------------------

## Function for predictions ------------------
# get predictions for each response, each quarter, male and female
# average male/female predictions

## working -------
#fix colnames to match pigsums
#this will go in function


######

AvgWashPreds<-function(gbm.mod,wash,opt.params,response){
  #split by quarter
    #grep colnames with q in name
  q_vars=grep("_q",colnames(wash))
  q1_vars=grep("_q1",colnames(wash))
  q2_vars=grep("_q2",colnames(wash))
  q3_vars=grep("_q3",colnames(wash))
  q4_vars=grep("_q4",colnames(wash))
  q_list=list(q1_vars,q2_vars,q3_vars,q4_vars)
  predmat=matrix(nrow=nrow(wash),ncol=4)
  for(q in 1:4){
  #colnames(predmat)[q]<-paste(response,q,sep="_")
  wash_notemporal=wash[,-q_list[[q]]]
  washq=cbind(wash_notemporal,wash[,q_list[[q]]])
  
  #fix names to match pigsums
  colnames(washq)=gsub("_q[0-9]","",colnames(washq))
  means=grep("_mn",colnames(washq))
  colnames(washq)=gsub("_mn","",colnames(washq))
  colnames(washq)[means]<-paste0("mean_",colnames(washq)[means])
  
  vars=grep("_var",colnames(washq))
  colnames(washq)=gsub("_var","",colnames(washq))
  colnames(washq)[vars]<-paste0("var_",colnames(washq)[vars])
  
  washq$season=as.factor(paste0("q",q))

  #Loop through periods
  predmat_p=matrix(nrow=nrow(wash),ncol=4)
  for(p in 1:4){
  washq$period=as.factor(paste0("p",p))
    
  washqf=washq
  washqm=washq
  washqf$sex=as.factor("Female")
  washqm$sex=as.factor("Male")
  
  pred.q.f <- predict.gbm(gbm.mod, washqf, n.trees=opt.params$n.trees, "response")
  pred.q.m <- predict.gbm(gbm.mod, washqm, n.trees=opt.params$n.trees, "response")
  
  predmat_p[,p]=rowMeans(cbind(pred.q.f,pred.q.m))
  
  }
  
  predmat[,q]=rowMeans(predmat_p)
  
  }
   return(predmat)
}

#Get predictions
pred_sl=AvgWashPreds(gbm.sl,wash,sl.opt.params.kfold,"sl")
pred_sigmasl=AvgWashPreds(gbm.sigma.sl,wash,sigma.sl.opt.params.kfold,"sigma_sl")
pred_disp=AvgWashPreds(gbm.disp,wash,disp.opt.params.kfold,"disp")
pred_sigmadisp=AvgWashPreds(gbm.sigma.disp,wash,sigma.disp.opt.params.kfold,"sigma_disp")
  
# Add mean predictions to watersheds --------------

#add preds to washq1
washq1$pred.sl<-pred.sl.q1
washq1$pred.sigma.sl<-pred.sigma.sl.q1
washq1$pred.disp<-pred.disp.q1
washq1$pred.sigma.disp<-pred.sigma.disp.q1

#add preds to washq2
washq2$pred.sl<-pred.sl.q2
washq2$pred.sigma.sl<-pred.sigma.sl.q2
washq2$pred.disp<-pred.disp.q2
washq2$pred.sigma.disp<-pred.sigma.disp.q2

#add preds to washq3
washq3$pred.sl<-pred.sl.q3
washq3$pred.sigma.sl<-pred.sigma.sl.q3
washq3$pred.disp<-pred.disp.q3
washq3$pred.sigma.disp<-pred.sigma.disp.q3

#add preds to washq4
washq4$pred.sl<-pred.sl.q4
washq4$pred.sigma.sl<-pred.sigma.sl.q4
washq4$pred.disp<-pred.disp.q4
washq4$pred.sigma.disp<-pred.sigma.disp.q4

#save the watershed files with the predictions
#write out sf objects with preds
st_write(washq1,paste0(home,"/Data/washq1_preds21APR23.shp"),append=FALSE)
st_write(washq2,paste0(home,"/Data/washq2_preds21APR23.shp"),append=FALSE)
st_write(washq3,paste0(home,"/Data/washq3_preds21APR23.shp"),append=FALSE)
st_write(washq4,paste0(home,"/Data/washq4_preds21APR23.shp"),append=FALSE)

# Plot prediction maps ------------------------
#get helper shapefiles
usplot <- st_read(dsn = "./Data/Shapefiles_Mapping/usmapplot", layer = "usmapplot_best")
usplot<-st_as_sf(usplot)
usplot <- st_cast(usplot, "MULTIPOLYGON")
usplot=st_transform(usplot,crs=st_crs(washq1))

#crop with usplot
washq1=st_intersection(washq1,ncplot)  
washq2=st_intersection(washq2,usplot)  
washq3=st_intersection(washq3,usplot)  
washq4=st_intersection(washq4,usplot)  

## Washq1 -------------------------------------

#washq1 maps
washq1.sl=ggplot() + 
  geom_sf(data=ncplot, fill="black")+
  geom_sf(data = washq1, aes(fill = pred.sl),lwd=0)+scale_fill_viridis_c(name="Mean step length (m) Jan-Mar")+
  geom_sf(data=ncplot, fill=NA, color="#EBEBEB")+
  theme_map()
washq1.sl

washq1.sigma.sl=ggplot() + 
  geom_sf(data=usplot, fill="black")+
  geom_sf(data = washq1, aes(fill = pred.sigma.sl),lwd=0)+scale_fill_viridis_c(name="Mean step length dispersion Jan-Mar")+
  geom_sf(data=usplot, fill=NA, color="#EBEBEB")+
  theme_map()

washq1.disp=ggplot() + 
  geom_sf(data=usplot, fill="black")+
  geom_sf(data = washq1, aes(fill = pred.disp),lwd=0)+scale_fill_viridis_c(name="Mean displacement (m) Jan-Mar")+
  geom_sf(data=usplot, fill=NA, color="#EBEBEB")+
  theme_map()

washq1.sigma.disp=ggplot() + 
  geom_sf(data=usplot, fill="black")+
  geom_sf(data = washq1, aes(fill = pred.sigma.disp),lwd=0)+scale_fill_viridis_c(name="Mean displacement dispersion Jan-Mar")+
  geom_sf(data=usplot, fill=NA, color="#EBEBEB")+
  theme_map()

washq1.tenavg=ggplot() + 
  geom_sf(data=usplot, fill="black")+
  geom_sf(data = washq1, aes(fill = pred.tenavg),lwd=0)+scale_fill_viridis_c(name="Mean of 10% longest step lengths (m) Jan-Mar")+
  geom_sf(data=usplot, fill=NA, color="#EBEBEB")+
  theme_map()

## Washq2 -------------------------------------

#washq2 maps
washq2.sl=ggplot() + 
  geom_sf(data=usplot, fill="black")+
  geom_sf(data = washq2, aes(fill = pred.sl),lwd=0)+scale_fill_viridis_c(name="Mean step length (m) Apr-Jun")+
  geom_sf(data=usplot, fill=NA, color="#EBEBEB")+
  theme_map()

washq2.sigma.sl=ggplot() + 
  geom_sf(data=usplot, fill="black")+
  geom_sf(data = washq2, aes(fill = pred.sigma.sl),lwd=0)+scale_fill_viridis_c(name="Mean step length dispersion Apr-Jun")+
  geom_sf(data=usplot, fill=NA, color="#EBEBEB")+
  theme_map()

washq2.disp=ggplot() + 
  geom_sf(data=usplot, fill="black")+
  geom_sf(data = washq2, aes(fill = pred.disp),lwd=0)+scale_fill_viridis_c(name="Mean displacement (m) Apr-Jun")+
  geom_sf(data=usplot, fill=NA, color="#EBEBEB")+
  theme_map()

washq2.sigma.disp=ggplot() + 
  geom_sf(data=usplot, fill="black")+
  geom_sf(data = washq2, aes(fill = pred.sigma.disp),lwd=0)+scale_fill_viridis_c(name="Mean displacement dispersion Apr-Jun")+
  geom_sf(data=usplot, fill=NA, color="#EBEBEB")+
  theme_map()

washq2.tenavg=ggplot() + 
  geom_sf(data=usplot, fill="black")+
  geom_sf(data = washq2, aes(fill = pred.tenavg),lwd=0)+scale_fill_viridis_c(name="Mean of 10% longest step lengths (m) Apr Jun")+
  geom_sf(data=usplot, fill=NA, color="#EBEBEB")+
  theme_map()

## Washq3 -------------------------------------

#washq3 maps
washq3.sl=ggplot() + 
  geom_sf(data=usplot, fill="black")+
  geom_sf(data = washq3, aes(fill = pred.sl),lwd=0)+scale_fill_viridis_c(name="Mean step length (m) Jul-Sep")+
  geom_sf(data=usplot, fill=NA, color="#EBEBEB")+
  theme_map()

washq3.sigma.sl=ggplot() + 
  geom_sf(data=usplot, fill="black")+
  geom_sf(data = washq3, aes(fill = pred.sigma.sl),lwd=0)+scale_fill_viridis_c(name="Mean step length dispersion Jul-Sep")+
  geom_sf(data=usplot, fill=NA, color="#EBEBEB")+
  theme_map()

washq3.disp=ggplot() + 
  geom_sf(data=usplot, fill="black")+
  geom_sf(data = washq3, aes(fill = pred.disp),lwd=0)+scale_fill_viridis_c(name="Mean displacement (m) Jul-Sep")+
  geom_sf(data=usplot, fill=NA, color="#EBEBEB")+
  theme_map()

washq3.sigma.disp=ggplot() + 
  geom_sf(data=usplot, fill="black")+
  geom_sf(data = washq3, aes(fill = pred.sigma.disp),lwd=0)+scale_fill_viridis_c(name="Mean displacement dispersion Jul-Sep")+
  geom_sf(data=usplot, fill=NA, color="#EBEBEB")+
  theme_map()

washq3.tenavg=ggplot() + 
  geom_sf(data=usplot, fill="black")+
  geom_sf(data = washq3, aes(fill = pred.tenavg),lwd=0)+scale_fill_viridis_c(name="Mean of 10% longest step lengths (m) Jul-Sep")+
  geom_sf(data=usplot, fill=NA, color="#EBEBEB")+
  theme_map()

## Washq4 -------------------------------------

#washq4 maps
washq4.sl=ggplot() + 
  geom_sf(data=usplot, fill="black")+
  geom_sf(data = washq4, aes(fill = pred.sl),lwd=0)+scale_fill_viridis_c(name="Mean step length (m) Oct-Dec")+
  geom_sf(data=usplot, fill=NA, color="#EBEBEB")+
  theme_map()

washq4.sigma.sl=ggplot() + 
  geom_sf(data=usplot, fill="black")+
  geom_sf(data = washq4, aes(fill = pred.sigma.sl),lwd=0)+scale_fill_viridis_c(name="Mean step length dispersion Oct-Dec")+
  geom_sf(data=usplot, fill=NA, color="#EBEBEB")+
  theme_map()

washq4.disp=ggplot() + 
  geom_sf(data=usplot, fill="black")+
  geom_sf(data = washq4, aes(fill = pred.disp),lwd=0)+scale_fill_viridis_c(name="Mean displacement (m) Oct-Dec")+
  geom_sf(data=usplot, fill=NA, color="#EBEBEB")+
  theme_map()

washq4.sigma.disp=ggplot() + 
  geom_sf(data=usplot, fill="black")+
  geom_sf(data = washq4, aes(fill = pred.sigma.disp),lwd=0)+scale_fill_viridis_c(name="Mean displacement dispersion Oct-Dec")+
  geom_sf(data=usplot, fill=NA, color="#EBEBEB")+
  theme_map()

washq4.tenavg=ggplot() + 
  geom_sf(data=usplot, fill="black")+
  geom_sf(data = washq4, aes(fill = pred.tenavg),lwd=0)+scale_fill_viridis_c(name="Mean of 10% longest step lengths (m) Oct-Dec")+
  geom_sf(data=usplot, fill=NA, color="#EBEBEB")+
  theme_map()

## Make map grids -------------------

plot_grid(washq1.sl, washq2.sl, washq3.sl, washq4.sl, nrow = 2)
plot_grid(washq1.sigma.sl, washq2.sigma.sl, washq3.sigma.sl, washq4.sigma.sl, nrow = 2)
plot_grid(washq1.disp, washq2.disp, washq3.disp, washq4.disp, nrow = 2)
plot_grid(washq1.sigma.disp, washq2.sigma.disp, washq3.sigma.disp, washq4.sigma.disp, nrow = 2)
plot_grid(washq1.tenavg, washq2.tenavg, washq3.tenavg, washq4.tenavg, nrow = 2)

# Plot uncertainty maps ---------------------------------

#create maps that show the variance in the predictions across each watershed

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


