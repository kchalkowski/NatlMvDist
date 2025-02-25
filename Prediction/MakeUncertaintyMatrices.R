
# Purpose --------------------------------

# This script does predictions from top models for sl/disp responses for all watershed/seasons

#Inputs: wash, pigsums
#Output:

# Vars --------------------------------
reps=2
clustnum=4
filestr<-"TestRuns"

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
library(abind)

## set working directories -------
#local:
home<-"/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/StatPigMvmt/Pipeline_R2"
outdir<-file.path(home,"4_Outputs")
objdir<-file.path(home,"2_Data","Objects")
if(!dir.exists(file.path(outdir,filestr))){dir.create(file.path(outdir,filestr))}
if(!dir.exists(file.path(outdir,filestr,"FigTab"))){dir.create(file.path(outdir,filestr,"FigTab"))}
if(!dir.exists(file.path(outdir,filestr,"UncPredMats"))){dir.create(file.path(outdir,filestr,"UncPredMats"))}

## read data -------

#read in pigsums dataset
pigsums<-readRDS(file.path(objdir,"dailyPigSums.rds"))
pigswsite<-readRDS(file.path(objdir,"geolocsnatl_wDispl.rds"))

#Read in formatted wastershed shapefile
wash<-readRDS(file.path(objdir,"wash_envcov_final.rds"))

#get helper shapefiles
usplot <- st_read(dsn = file.path(objdir,"usmapplot_best.shp"), layer = "usmapplot_best")

#Need do this again bc this script runs on hpc
# Make model selection table ---------------------------
#Need do this again bc this script runs on hpc

#make table with all mean R2 and RMSEs
#row for each response (5)
#col for each 2 metrics (RMSE and R2) and 4 models (full, drop01, lassodrop, null)

#initiate tables, region and random cv methods
model.sel.tbl=as.data.frame(matrix(nrow=4,ncol=6))
rownames(model.sel.tbl)=c("step length", "step length sigma", "displacement", "disp. sigma")
colnames(model.sel.tbl)=c("full RMSE", "full r2", "drop 01 RMSE", "drop 01 r2", "null RMSE", "null R2")
model.sel.tbl.region=as.data.frame(matrix(nrow=4,ncol=6))
rownames(model.sel.tbl.region)=c("step length", "step length sigma", "displacement", "disp. sigma")
colnames(model.sel.tbl.region)=c("full RMSE", "full r2", "drop 01 RMSE", "drop 01 r2", "null RMSE", "null R2")

#Full model: get values from file
full.RMSE.sl=read.csv(file.path(outdir,filestr,"Random","sl_meanRMSE.csv"))
full.R2.sl=read.csv(file.path(outdir,filestr,"Random","sl_meanR2.csv"))
full.RMSE.sigma.sl=read.csv(file.path(outdir,filestr,"Random","sigma.sl_meanRMSE.csv"))
full.R2.sigma.sl=read.csv(file.path(outdir,filestr,"Random","sigma.sl_meanR2.csv"))
full.RMSE.disp=read.csv(file.path(outdir,filestr,"Random","disp_meanRMSE.csv"))
full.R2.disp=read.csv(file.path(outdir,filestr,"Random","disp_meanR2.csv"))
full.RMSE.sigma.disp=read.csv(file.path(outdir,filestr,"Random","sigma.disp_meanRMSE.csv"))
full.R2.sigma.disp=read.csv(file.path(outdir,filestr,"Random","sigma.disp_meanR2.csv"))

#Full model region kfold: get values from file
full.RMSE.sl.region=read.csv(file.path(outdir,filestr,"Region","sl_meanRMSE.csv"))
full.R2.sl.region=read.csv(file.path(outdir,filestr,"Region","sl_meanR2.csv"))
full.RMSE.sigma.sl.region=read.csv(file.path(outdir,filestr,"Region","sigma.sl_meanRMSE.csv"))
full.R2.sigma.sl.region=read.csv(file.path(outdir,filestr,"Region","sigma.sl_meanR2.csv"))
full.RMSE.disp.region=read.csv(file.path(outdir,filestr,"Region","disp_meanRMSE.csv"))
full.R2.disp.region=read.csv(file.path(outdir,filestr,"Region","disp_meanR2.csv"))
full.RMSE.sigma.disp.region=read.csv(file.path(outdir,filestr,"Region","sigma.disp_meanRMSE.csv"))
full.R2.sigma.disp.region=read.csv(file.path(outdir,filestr,"Region","sigma.disp_meanR2.csv"))

#Drop01 model: get vals
drop01.RMSE.sl=read.csv(file.path(outdir,filestr,"Random","GBM_01Drop","sl_meanRMSE.csv"))
drop01.R2.sl=read.csv(file.path(outdir,filestr,"Random","GBM_01Drop","sl_meanR2.csv"))
drop01.RMSE.sigma.sl=read.csv(file.path(outdir,filestr,"Random","GBM_01Drop","sigma.sl_meanRMSE.csv"))
drop01.R2.sigma.sl=read.csv(file.path(outdir,filestr,"Random","GBM_01Drop","sigma.sl_meanR2.csv"))
drop01.RMSE.disp=read.csv(file.path(outdir,filestr,"Random","GBM_01Drop","disp_meanRMSE.csv"))
drop01.R2.disp=read.csv(file.path(outdir,filestr,"Random","GBM_01Drop","disp_meanR2.csv"))
drop01.RMSE.sigma.disp=read.csv(file.path(outdir,filestr,"Random","GBM_01Drop","sigma.disp_meanRMSE.csv"))
drop01.R2.sigma.disp=read.csv(file.path(outdir,filestr,"Random","GBM_01Drop","sigma.disp_meanR2.csv"))

#Drop01 region kfold model: get vals
drop01.RMSE.sl.region=read.csv(file.path(outdir,filestr,"Region","GBM_01Drop","sl_meanRMSE.csv"))
drop01.R2.sl.region=read.csv(file.path(outdir,filestr,"Region","GBM_01Drop","sl_meanR2.csv"))
drop01.RMSE.sigma.sl.region=read.csv(file.path(outdir,filestr,"Region","GBM_01Drop","sigma.sl_meanRMSE.csv"))
drop01.R2.sigma.sl.region=read.csv(file.path(outdir,filestr,"Region","GBM_01Drop","sigma.sl_meanR2.csv"))
drop01.RMSE.disp.region=read.csv(file.path(outdir,filestr,"Region","GBM_01Drop","disp_meanRMSE.csv"))
drop01.R2.disp.region=read.csv(file.path(outdir,filestr,"Region","GBM_01Drop","disp_meanR2.csv"))
drop01.RMSE.sigma.disp.region=read.csv(file.path(outdir,filestr,"Region","GBM_01Drop","sigma.disp_meanRMSE.csv"))
drop01.R2.sigma.disp.region=read.csv(file.path(outdir,filestr,"Region","GBM_01Drop","sigma.disp_meanR2.csv"))

#null model: get vals
null.RMSE.sl=read.csv(file.path(outdir,filestr,"Random","GBM_Null","sl_meanRMSE.csv"))
null.R2.sl=read.csv(file.path(outdir,filestr,"Random","GBM_Null/sl_meanR2.csv"))
null.RMSE.sigma.sl=read.csv(file.path(outdir,filestr,"Random","GBM_Null","sigma.sl_meanRMSE.csv"))
null.R2.sigma.sl=read.csv(file.path(outdir,filestr,"Random","GBM_Null","sigma.sl_meanR2.csv"))
null.RMSE.disp=read.csv(file.path(outdir,filestr,"Random","GBM_Null","disp_meanRMSE.csv"))
null.R2.disp=read.csv(file.path(outdir,filestr,"Random","GBM_Null","disp_meanR2.csv"))
null.RMSE.sigma.disp=read.csv(file.path(outdir,filestr,"Random","GBM_Null","sigma.disp_meanRMSE.csv"))
null.R2.sigma.disp=read.csv(file.path(outdir,filestr,"Random","GBM_Null","sigma.disp_meanR2.csv"))

#null model region: get vals
null.RMSE.sl.region=read.csv(file.path(outdir,filestr,"Region","GBM_Null","sl_meanRMSE.csv"))
null.R2.sl.region=read.csv(file.path(outdir,filestr,"Region","GBM_Null","sl_meanR2.csv"))
null.RMSE.sigma.sl.region=read.csv(file.path(outdir,filestr,"Region","GBM_Null","sigma.sl_meanRMSE.csv"))
null.R2.sigma.sl.region=read.csv(file.path(outdir,filestr,"Region","GBM_Null","sigma.sl_meanR2.csv"))
null.RMSE.disp.region=read.csv(file.path(outdir,filestr,"Region","GBM_Null","disp_meanRMSE.csv"))
null.R2.disp.region=read.csv(file.path(outdir,filestr,"Region","GBM_Null","disp_meanR2.csv"))
null.RMSE.sigma.disp.region=read.csv(file.path(outdir,filestr,"Region","GBM_Null","sigma.disp_meanRMSE.csv"))
null.R2.sigma.disp.region=read.csv(file.path(outdir,filestr,"Region","GBM_Null","sigma.disp_meanR2.csv"))

#col 1,2: full model rmse, r2
#row 1: sl
model.sel.tbl[1,1]=round(full.RMSE.sl[1,2],3)
model.sel.tbl[1,2]=round(full.R2.sl[1,2],3)
#row 2: sigma sl
model.sel.tbl[2,1]=round(full.RMSE.sigma.sl[1,2],3)
model.sel.tbl[2,2]=round(full.R2.sigma.sl[1,2],3)
#row 3: disp
model.sel.tbl[3,1]=round(full.RMSE.disp[1,2],3)
model.sel.tbl[3,2]=round(full.R2.disp[1,2],3)
#row 4: sigma disp
model.sel.tbl[4,1]=round(full.RMSE.sigma.disp[1,2],3)
model.sel.tbl[4,2]=round(full.R2.sigma.disp[1,2],3)

#col 3,4: drop01 rmse, r2
#row 1: sl
model.sel.tbl[1,3]=round(drop01.RMSE.sl[1,2],3)
model.sel.tbl[1,4]=round(drop01.R2.sl[1,2],3)
#row 2: sigma sl
model.sel.tbl[2,3]=round(drop01.RMSE.sigma.sl[1,2],3)
model.sel.tbl[2,4]=round(drop01.R2.sigma.sl[1,2],3)
#row 3: disp
model.sel.tbl[3,3]=round(drop01.RMSE.disp[1,2],3)
model.sel.tbl[3,4]=round(full.R2.disp[1,2],3)
#row 4: sigma disp
model.sel.tbl[4,3]=round(drop01.RMSE.sigma.disp[1,2],3)
model.sel.tbl[4,4]=round(drop01.R2.sigma.disp[1,2],3)

#col 7,8: null model rmse, r2
#row 1: sl
model.sel.tbl[1,7]=round(null.RMSE.sl[1,2],3)
model.sel.tbl[1,8]=round(null.R2.sl[1,2],3)
#row 2: sigma sl
model.sel.tbl[2,7]=round(null.RMSE.sigma.sl[1,2],3)
model.sel.tbl[2,8]=round(null.R2.sigma.sl[1,2],3)
#row 3: disp
model.sel.tbl[3,7]=round(null.RMSE.disp[1,2],3)
model.sel.tbl[3,8]=round(null.R2.disp[1,2],3)
#row 4: sigma disp
model.sel.tbl[4,7]=round(null.RMSE.sigma.disp[1,2],3)
model.sel.tbl[4,8]=round(null.R2.sigma.disp[1,2],3)

#col 1,2: region full model rmse, r2
#row 1: sl
model.sel.tbl.region[1,1]=round(full.RMSE.sl.region[1,2],3)
model.sel.tbl.region[1,2]=round(full.R2.sl.region[1,2],3)
#row 2: sigma sl
model.sel.tbl.region[2,1]=round(full.RMSE.sigma.sl.region[1,2],3)
model.sel.tbl.region[2,2]=round(full.R2.sigma.sl.region[1,2],3)
#row 3: disp
model.sel.tbl.region[3,1]=round(full.RMSE.disp.region[1,2],3)
model.sel.tbl.region[3,2]=round(full.R2.disp.region[1,2],3)
#row 4: sigma disp
model.sel.tbl.region[4,1]=round(full.RMSE.sigma.disp.region[1,2],3)
model.sel.tbl.region[4,2]=round(full.R2.sigma.disp.region[1,2],3)

#col 3,4: region drop01 model rmse, r2
#row 1: sl
model.sel.tbl.region[1,3]=round(drop01.RMSE.sl.region[1,2],3)
model.sel.tbl.region[1,4]=round(drop01.R2.sl.region[1,2],3)
#row 2: sigma sl
model.sel.tbl.region[2,3]=round(drop01.RMSE.sigma.sl.region[1,2],3)
model.sel.tbl.region[2,4]=round(drop01.R2.sigma.sl.region[1,2],3)
#row 3: disp
model.sel.tbl.region[3,3]=round(drop01.RMSE.disp.region[1,2],3)
model.sel.tbl.region[3,4]=round(drop01.R2.disp.region[1,2],3)
#row 4: sigma disp
model.sel.tbl.region[4,3]=round(drop01.RMSE.sigma.disp.region[1,2],3)
model.sel.tbl.region[4,4]=round(drop01.R2.sigma.disp.region[1,2],3)

#col 7,8: region null model rmse, r2
#row 1: sl
model.sel.tbl.region[1,5]=round(null.RMSE.sl.region[1,2],3)
model.sel.tbl.region[1,6]=round(null.R2.sl.region[1,2],3)
#row 2: sigma sl
model.sel.tbl.region[2,5]=round(null.RMSE.sigma.sl.region[1,2],3)
model.sel.tbl.region[2,6]=round(null.R2.sigma.sl.region[1,2],3)
#row 3: disp
model.sel.tbl.region[3,5]=round(null.RMSE.disp.region[1,2],3)
model.sel.tbl.region[3,6]=round(null.R2.disp.region[1,2],3)
#row 4: sigma disp
model.sel.tbl.region[4,5]=round(null.RMSE.sigma.disp.region[1,2],3)
model.sel.tbl.region[4,6]=round(null.R2.sigma.disp.region[1,2],3)

model.sel.tbl.total=as.data.frame(matrix(nrow=4,ncol=12))
rownames(model.sel.tbl.total)=c("step length", "step length sigma", "displacement", "disp. sigma")
colnames(model.sel.tbl.total)=c("full RMSE random", "full r2 random", "drop 01 RMSE random", "drop 01 r2 random","null RMSE random", "null R2 random",
                                "full RMSE region", "full r2 region", "drop 01 RMSE region", "drop 01 r2 region", "null RMSE region", "null R2region")

model.sel.tbl.total[1:4,1:6]=model.sel.tbl
model.sel.tbl.total[1:4,9:12]=model.sel.tbl.region

write.csv(model.sel.tbl.total,file.path(outdir,filestr,"FigTab","ModelSelTblTotal.csv"))

# Make xsel/xvec -------

getXvec<-function(model.sel.tbl.total,out.opt){
  X_vec_list=vector(mode="list",length=4)
  names(X_vec_list)<-c("Xsl","sigmasl","Xdisp","sigmadisp")
  
  X_sel=data.frame(matrix(nrow=4,ncol=4))
  colnames(X_sel)<-c("response","X_vec","reg_ran","vars")
  
  for(i in 1:nrow(model.sel.tbl.total)){
    X_sel[i,1]=rownames(model.sel.tbl.total)[i]
    c_ind=which(model.sel.tbl.total[i,]==min(model.sel.tbl.total[i,]))
    X_sel[i,2]=colnames(model.sel.tbl.total)[c_ind]
  }
  
  #parse text in X sel: region or random
  X_sel[grep("region",X_sel$X_vec),3]<-"Region"
  X_sel[grep("random",X_sel$X_vec),3]<-"Random"
  
  #parse text in X sel: variable sel
  X_sel[grep("full",X_sel$X_vec),4]<-"Full"
  X_sel[grep("drop 01",X_sel$X_vec),4]<-"X_vec_01Drop"
  X_sel[grep("lasso",X_sel$X_vec),4]<-"X_vec_postlasso"
  X_sel[grep("null",X_sel$X_vec),4]<-"Null"
  
  if(out.opt=="X_sel"){
    return(X_sel)
  } else{
    
    for(i in 1:length(X_vec_list)){
      if(X_sel[i,4]=="Full"){
        X_vec_list[[i]]<-c(
          which(colnames(pigsums)=="sex"),
          which(colnames(pigsums)=="season"),
          which(colnames(pigsums)=="period"),
          which(colnames(pigsums)=="mean_tc"):
            which(colnames(pigsums)=="var_lc_24")
        )
      }
      if(X_sel[i,4]=="X_vec_01Drop"|X_sel[i,4]=="X_vec_postlasso"){
        X_vec_folder=file.path(outdir,filestr,X_sel[i,3],X_sel[i,4])
        Xfiles=list.files(X_vec_folder,full.names=TRUE)
        Xfile=Xfiles[grep(names(X_vec_list)[i],Xfiles)]
        X_vec_list[[i]]=read.csv(Xfile)[,2]
      }
      
      if(X_sel[i,4]=="Null"){
        X_vec_list[[i]]=rep(1,nrow(pigsums))
      }
    }
    
    return(X_vec_list)
  }
}


X_sel=getXvec(model.sel.tbl.total,"X_sel")

sl.opt.params.kfold=read.csv(file.path(outdir,filestr,X_sel[1,3],"sl_bestmodelparams.csv"))
sigma.sl.opt.params.kfold=read.csv(file.path(outdir,filestr,X_sel[2,3],"sigma.sl_bestmodelparams.csv"))
disp.opt.params.kfold=read.csv(file.path(outdir,filestr,X_sel[3,3],"disp_bestmodelparams.csv"))
sigma.disp.opt.params.kfold=read.csv(file.path(outdir,filestr,X_sel[4,3],"sigma.disp_bestmodelparams.csv"))

X_vec_list=getXvec(model.sel.tbl.total,"Xlist")

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
                       pigsums_sigmadisp,wash,1000,
                       sigma.disp.opt.params.kfold,
                       X_vec_list$Xdisp,
                       "displ_disp","gaussian",
                       path.out,"sigmadisp_pmat")

AverageQuarters(path.out,"sigmadisp_pmat")



