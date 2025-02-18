home="homepath"

##Identifying optimal parameters for gradient boosted models

#created using tutorials:
#https://bradleyboehmke.github.io/HOML/gbm.html
#https://rspatial.org/raster/sdm/9_sdm_brt.html

###################
## Set Variables ##
###################
#vars that may need to be changed from run to run
repname="10DEC24_Runs"
split_type="Region"

###################
##Process Outline##
###################

#input: pigsums.csv

################
##Script Setup##
################
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
require(job)

#read data
pigsums=readRDS(file.path(objdir,"dailyPigSums.rds", fsep = .Platform$file.sep))

###############
##Format Data##
###############
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

###################################
##Source Functions and Set Params##
###################################

source(file.path(gbm_funcdir,"/GBM_FunctionSourcer.R", fsep = .Platform$file.sep))

###############################
##--Hyperparameter optimization--##
###############################

#job::job({
  #resource needed parms

  filestr=paste("4_Outputs",repname,split_type,sep="/")
  if(!dir.exists(file.path("4_Outputs/",repname, fsep = .Platform$file.sep))){dir.create(file.path("4_Outputs/",repname, fsep = .Platform$file.sep))}
  if(!dir.exists(filestr)){dir.create(filestr)}
  
  #random=randomly splits the data into test/training sets according to ko/ki
  split_type="Region" #splits the data according to 'region', which is the study name
  
#for testing
ko_t=10 #outer k-fold cross validations
ki_t=3 #inner k-fold cross validations

#X vec is a vector of indices with parameters included in the model
X_vec.start=c(
  which(colnames(pigsums)=="sex"),
  which(colnames(pigsums)=="season"),
  which(colnames(pigsums)=="period"),
  which(colnames(pigsums)=="mean_tc"):
  which(colnames(pigsums)=="var_lc_24")
        )

#Run models
#res.sl=Run.GBM.Model(pigsums_sl,"sl_mean",X_vec.start,ko_t,ki_t,"poisson",split_type,ntreemax=8000)

#Export CV stats from first run
#set file.str and make directory for saving
#write.csv(res.sl[[2]],paste0(filestr,"/sl_meanRMSE.csv"))
#write.csv(res.sl[[3]],paste0(filestr,"/sl_bestmodelparams.csv"))
#write.csv(res.sl[[4]],paste0(filestr,"/sl_meanR2.csv"))
#saveRDS(res.sl[[5]],paste0(filestr,"/sl_preds.rds"))
#saveRDS(res.sl[[6]],paste0(filestr,"/sl_testobs.rds"))

#started res.sigma.sl run at aug 26/2024, 12pm
#res.sigma.sl=Run.GBM.Model(pigsums_sl,"sl_disp",X_vec.start,ko_t,ki_t,"gaussian",split_type,ntreemax=8000)

#write.csv(res.sigma.sl[[2]],paste0(filestr,"/sigma.sl_meanRMSE.csv"))
#write.csv(res.sigma.sl[[3]],paste0(filestr,"/sigma.sl_bestmodelparams.csv"))
#write.csv(res.sigma.sl[[4]],paste0(filestr,"/sigma.sl_meanR2.csv"))
#saveRDS(res.sigma.sl[[5]],paste0(filestr,"/sigma.sl_preds.rds"))
#saveRDS(res.sigma.sl[[6]],paste0(filestr,"/sigma.sl_testobs.rds"))

res.disp=Run.GBM.Model(pigsums_displ,"displ_mean",X_vec.start,ko_t,ki_t,"poisson",split_type,ntreemax=8000)

write.csv(res.disp[[2]],paste0(filestr,"/disp_meanRMSE.csv"))
write.csv(res.disp[[3]],paste0(filestr,"/disp_bestmodelparams.csv"))
write.csv(res.disp[[4]],paste0(filestr,"/disp_meanR2.csv"))
saveRDS(res.disp[[5]],paste0(filestr,"/disp_preds.rds"))
saveRDS(res.disp[[6]],paste0(filestr,"/disp_testobs.rds"))

res.sigma.disp=Run.GBM.Model(pigsums_displ,"displ_disp",X_vec.start,ko_t,ki_t,"gaussian",split_type,ntreemax=8000)

write.csv(res.sigma.disp[[2]],paste0(filestr,"/sigma.disp_meanRMSE.csv"))
write.csv(res.sigma.disp[[3]],paste0(filestr,"/sigma.disp_bestmodelparams.csv"))
write.csv(res.sigma.disp[[4]],paste0(filestr,"/sigma.disp_meanR2.csv"))
saveRDS(res.sigma.disp[[5]],paste0(filestr,"/sigma.disp_preds.rds"))
saveRDS(res.sigma.disp[[6]],paste0(filestr,"/sigma.disp_testobs.rds"))

#})

#Get x_vecs without unimportant variables
#Remove variables that had less than 0.1 influence
X_vec.2.sl=drop.unimportant(res.sl[[1]],X_vec.start,0.1)
X_vec.2.sigma.sl=drop.unimportant(res.sigma.sl[[1]],X_vec.start,0.1)
X_vec.2.disp=drop.unimportant(res.disp[[1]],X_vec.start,0.1)
X_vec.2.sigma.disp=drop.unimportant(res.sigma.disp[[1]],X_vec.start,0.1)
X_vec.2.tenavg=drop.unimportant(res.tenavg[[1]],X_vec.start,0.1)

#Save new vectors to file
write.csv(X_vec.2.sl,paste0(filestr,"/X_vec_01Drop/","Xsl.csv"))
write.csv(X_vec.2.sigma.sl,paste0(filestr,"/X_vec_01Drop/Xsigma.sl.csv"))
write.csv(X_vec.2.disp,paste0(filestr,"/X_vec_01Drop/Xdisp.csv"))
write.csv(X_vec.2.sigma.disp,paste0(filestr,"/X_vec_01Drop/Xsigma.disp.csv"))
write.csv(X_vec.2.tenavg,paste0(filestr,"/X_vec_01Drop/Xtenavg.csv"))

#Run models without dropped variables
#job::job({
res.sl=Run.GBM.Model(pigsums,"sl_",X_vec.2.sl,ko_t,ki_t,"poisson",split_type,ntreemax=8000)
res.sigma.sl=Run.GBM.Model(pigsums,"sigma_sl",X_vec.2.sigma.sl,ko_t,ki_t,"gaussian",split_type,ntreemax=8000)
res.disp=Run.GBM.Model(pigsums,"displacement",X_vec.2.disp,ko_t,ki_t,"poisson",split_type,ntreemax=8000)
res.sigma.disp=Run.GBM.Model(pigsums,"sigma_disp",X_vec.2.sigma.disp,ko_t,ki_t,"gaussian",split_type,ntreemax=8000)

#Export CV stats
write.csv(res.sl[[2]],paste0(filestr,"/GBM_01Drop/sl_meanRMSE.csv"))
write.csv(res.sl[[3]],paste0(filestr,"/GBM_01Drop/sl_bestmodelparams.csv"))
write.csv(res.sl[[4]],paste0(filestr,"/GBM_01Drop/sl_meanR2.csv"))
saveRDS(res.sl[[5]],paste0(filestr,"/GBM_01Drop/sl_preds.rds"))
saveRDS(res.sl[[6]],paste0(filestr,"/GBM_01Drop/sl_testobs.rds"))

write.csv(res.sigma.sl[[2]],paste0(filestr,"/GBM_01Drop/sigma.sl_meanRMSE.csv"))
write.csv(res.sigma.sl[[3]],paste0(filestr,"/GBM_01Drop/sigma.sl_bestmodelparams.csv"))
write.csv(res.sigma.sl[[4]],paste0(filestr,"/GBM_01Drop/sigma.sl_meanR2.csv"))
saveRDS(res.sigma.sl[[5]],paste0(filestr,"/GBM_01Drop/sigma.sl_preds.rds"))
saveRDS(res.sigma.sl[[6]],paste0(filestr,"/GBM_01Drop/sigma.sl_testobs.rds"))

write.csv(res.disp[[2]],paste0(filestr,"/GBM_01Drop/disp_meanRMSE.csv"))
write.csv(res.disp[[3]],paste0(filestr,"/GBM_01Drop/disp_bestmodelparams.csv"))
write.csv(res.disp[[4]],paste0(filestr,"/GBM_01Drop/disp_meanR2.csv"))
saveRDS(res.disp[[5]],paste0(filestr,"/GBM_01Drop/disp_preds.rds"))
saveRDS(res.disp[[6]],paste0(filestr,"/GBM_01Drop/disp_testobs.rds"))

write.csv(res.sigma.disp[[2]],paste0(filestr,"/GBM_01Drop/sigma.disp_meanRMSE.csv"))
write.csv(res.sigma.disp[[3]],paste0(filestr,"/GBM_01Drop/sigma.disp_bestmodelparams.csv"))
write.csv(res.sigma.disp[[4]],paste0(filestr,"/GBM_01Drop/sigma.disp_meanR2.csv"))
saveRDS(res.sigma.disp[[5]],paste0(filestr,"/GBM_01Drop/sigma.disp_preds.rds"))
saveRDS(res.sigma.disp[[6]],paste0(filestr,"/GBM_01Drop/sigma.disp_testobs.rds"))

#})

#Run lasso function for each model
#get list of remaining important variables to use
nrep=100
typemeasure="mae"

#some tidying of x, cant have characters
x=pigsums[,X_vec.start]
x$sexf=NA
x[x$sex=="Male",]$sexf=1
x[x$sex=="Female",]$sexf=0
x=x[,-2]

#class(x$sexf)
#make alpha search grid
grid=seq(0.01,1,by=0.01)

#job::job({
  source(paste0(home,"/Scripts_Polished/GBM_nestedCV/run_Lasso_function.R"))
  source(paste0(home,"/Scripts_Polished/GBM_nestedCV/Optimize_Alpha_Lasso.R"))
  
#Run lasso on each model
keepVars.sl=runLassofunc(x,pigsums$sl_,"poisson",grid,typemeasure,nrep)
X_vec.sl.lasso=which(colnames(pigsums)%in%names(keepVars.sl))

keepVars.sigma.sl=runLassofunc(x,pigsums$sigma_sl,"poisson",grid,typemeasure,nrep)
X_vec.sigma.sl.lasso=which(colnames(pigsums)%in%names(keepVars.sigma.sl))

keepVars.disp=runLassofunc(x,pigsums$displacement,"poisson",grid,typemeasure,nrep)
X_vec.disp.lasso=which(colnames(pigsums)%in%names(keepVars.disp))

keepVars.sigma.disp=runLassofunc(x,pigsums$sigma_disp,"poisson",grid,typemeasure,nrep)
X_vec.sigma.disp.lasso=which(colnames(pigsums)%in%names(keepVars.sigma.disp))

#Write out the post lasso X vecs
write.csv(X_vec.sl.lasso,paste0(filestr,"/X_vec_postlasso/Xsl_lasso.csv"))
write.csv(X_vec.sigma.sl.lasso,paste0(filestr,"/X_vec_postlasso/Xsigmasl_lasso.csv"))
write.csv(X_vec.disp.lasso,paste0(filestr,"/X_vec_postlasso/Xdisp_lasso.csv"))
write.csv(X_vec.sigma.disp.lasso,paste0(filestr,"/X_vec_postlasso/Xsigmadisp_lasso.csv"))

#})


#Run lasso models
#job::job({

  source(paste0(home,"/Scripts_Polished/GBM_nestedCV/K_Split.R"))
  source(paste0(home,"/Scripts_Polished/GBM_nestedCV/HyperparamModelFit.R"))
  source(paste0(home,"/Scripts_Polished/GBM_nestedCV/InnerCrossValidationLoop.R"))
  source(paste0(home,"/Scripts_Polished/GBM_nestedCV/OuterCrossValidationLoop.R"))
  source(paste0(home,"/Scripts_Polished/GBM_nestedCV/RunGBMModelFunction.R"))
  
  res.sl=Run.GBM.Model(pigsums,"sl_",X_vec.sl.lasso,ko_t,ki_t,"poisson",split_type,ntreemax=8000)
  
  #Export CV stats
  write.csv(res.sl[[2]],paste0(filestr,"/GBM_Lasso/sl_meanRMSE.csv"))
  write.csv(res.sl[[3]],paste0(filestr,"/GBM_Lasso/sl_bestmodelparams.csv"))
  write.csv(res.sl[[4]],paste0(filestr,"/GBM_Lasso/sl_meanR2.csv"))
  saveRDS(res.sl[[5]],paste0(filestr,"/GBM_Lasso/sl_preds.rds"))
  saveRDS(res.sl[[6]],paste0(filestr,"/GBM_Lasso/sl_testobs.rds"))
  
  res.sigma.sl=Run.GBM.Model(pigsums,"sigma_sl",X_vec.sigma.sl.lasso,ko_t,ki_t,"gaussian",split_type,ntreemax=8000)
  
  #Export CV stats
  write.csv(res.sigma.sl[[2]],paste0(filestr,"/GBM_Lasso/sigma.sl_meanRMSE.csv"))
  write.csv(res.sigma.sl[[3]],paste0(filestr,"/GBM_Lasso/sigma.sl_bestmodelparams.csv"))
  write.csv(res.sigma.sl[[4]],paste0(filestr,"/GBM_Lasso/sigma.sl_meanR2.csv"))
  saveRDS(res.sigma.sl[[5]],paste0(filestr,"/GBM_Lasso/sigma.sl_preds.rds"))
  saveRDS(res.sigma.sl[[6]],paste0(filestr,"/GBM_Lasso/sigma.sl_testobs.rds"))
  
  res.disp=Run.GBM.Model(pigsums,"displacement",X_vec.disp.lasso,ko_t,ki_t,"poisson",split_type,ntreemax=8000)
  
  #Export CV stats
  write.csv(res.disp[[2]],paste0(filestr,"/GBM_Lasso/disp_meanRMSE.csv"))
  write.csv(res.disp[[3]],paste0(filestr,"/GBM_Lasso/disp_bestmodelparams.csv"))
  write.csv(res.disp[[4]],paste0(filestr,"/GBM_Lasso/disp_meanR2.csv"))
  saveRDS(res.disp[[5]],paste0(filestr,"/GBM_Lasso/disp_preds.rds"))
  saveRDS(res.disp[[6]],paste0(filestr,"/GBM_Lasso/disp_testobs.rds"))
  
  res.sigma.disp=Run.GBM.Model(pigsums,"sigma_disp",X_vec.sigma.disp.lasso,ko_t,ki_t,"gaussian",split_type,ntreemax=8000)
  
  #Export CV stats
  write.csv(res.sigma.disp[[2]],paste0(filestr,"/GBM_Lasso/sigma.disp_meanRMSE.csv"))
  write.csv(res.sigma.disp[[3]],paste0(filestr,"/GBM_Lasso/sigma.disp_bestmodelparams.csv"))
  write.csv(res.sigma.disp[[4]],paste0(filestr,"/GBM_Lasso/sigma.disp_meanR2.csv"))
  saveRDS(res.sigma.disp[[5]],paste0(filestr,"/GBM_Lasso/sigma.disp_preds.rds"))
  saveRDS(res.sigma.disp[[6]],paste0(filestr,"/GBM_Lasso/sigma.disp_testobs.rds"))
  
#})

#job::job({
  pigsums$gbm.null=rep(1,nrow(pigsums))
  X_null=ncol(pigsums)
  
  res.sl=Run.GBM.Model(pigsums,"sl_",X_null,ko_t,ki_t,"poisson",split_type,ntreemax=8000)
  
  #Export CV stats
  write.csv(res.sl[[2]],paste0(filestr,"/GBM_Null/sl_meanRMSE.csv"))
  write.csv(res.sl[[3]],paste0(filestr,"/GBM_Null/sl_bestmodelparams.csv"))
  write.csv(res.sl[[4]],paste0(filestr,"/GBM_Null/sl_meanR2.csv"))
  
  res.sigma.sl=Run.GBM.Model(pigsums,"sigma_sl",X_null,ko_t,ki_t,"gaussian",split_type,ntreemax=8000)
  
  #Export CV stats
  write.csv(res.sigma.sl[[2]],paste0(filestr,"/GBM_Null/sigma.sl_meanRMSE.csv"))
  write.csv(res.sigma.sl[[3]],paste0(filestr,"/GBM_Null/sigma.sl_bestmodelparams.csv"))
  write.csv(res.sigma.sl[[4]],paste0(filestr,"/GBM_Null/sigma.sl_meanR2.csv"))
  
  res.disp=Run.GBM.Model(pigsums,"displacement",X_null,ko_t,ki_t,"poisson",split_type,ntreemax=8000)
  
  #Export CV stats
  write.csv(res.disp[[2]],paste0(filestr,"/GBM_Null/disp_meanRMSE.csv"))
  write.csv(res.disp[[3]],paste0(filestr,"/GBM_Null/disp_bestmodelparams.csv"))
  write.csv(res.disp[[4]],paste0(filestr,"/GBM_Null/disp_meanR2.csv"))
  
  res.sigma.disp=Run.GBM.Model(pigsums,"sigma_disp",X_null,ko_t,ki_t,"gaussian",split_type,ntreemax=8000)
  
  #Export CV stats
  write.csv(res.sigma.disp[[2]],paste0(filestr,"/GBM_Null/sigma.disp_meanRMSE.csv"))
  write.csv(res.sigma.disp[[3]],paste0(filestr,"/GBM_Null/sigma.disp_bestmodelparams.csv"))
  write.csv(res.sigma.disp[[4]],paste0(filestr,"/GBM_Null/sigma.disp_meanR2.csv"))

#})

