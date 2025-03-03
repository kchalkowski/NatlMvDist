
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
filestr<-"3MAR25_Runs"
outdir<-file.path(home,"4_Outputs")
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
sl.opt.params.kfold=read.csv(file.path(outdir,filestr,X_sel[1,3],"GBM_01Drop","sl_bestmodelparams.csv"))
sigma.sl.opt.params.kfold=read.csv(file.path(outdir,filestr,X_sel[2,3],"sigma.sl_bestmodelparams.csv"))
disp.opt.params.kfold=read.csv(file.path(outdir,filestr,X_sel[3,3],"GBM_01Drop","disp_bestmodelparams.csv"))
sigma.disp.opt.params.kfold=read.csv(file.path(outdir,filestr,X_sel[4,3],"sigma.disp_bestmodelparams.csv"))

#get helper shapefiles
usplot <- st_read(dsn = file.path(objdir,"usmapplot_best.shp"), layer = "usmapplot_best")

## format data -------

#format pigsums data
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

#format helper shapefiles for plotting
usplot<-st_as_sf(usplot)
usplot <- st_cast(usplot, "MULTIPOLYGON")
usplot=st_transform(usplot,crs=st_crs(wash))

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

# Predictions -----------------------

## Function for predictions ------------------
# get predictions for each response, each quarter, male and female
# average male/female predictions

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
    print(paste0("quartr ",q))

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
    print(paste0("period ",p))
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

## Get predictions -----------------------
pred_sl=AvgWashPreds(gbm.sl,wash,sl.opt.params.kfold,"sl")
pred_sigmasl=AvgWashPreds(gbm.sigma.sl,wash,sigma.sl.opt.params.kfold,"sigma_sl")
pred_disp=AvgWashPreds(gbm.disp,wash,disp.opt.params.kfold,"disp")
pred_sigmadisp=AvgWashPreds(gbm.sigma.disp,wash,sigma.disp.opt.params.kfold,"sigma_disp")
  
# Add mean predictions to watersheds --------------

washp=wash
washp=cbind(washp,pred_sl)
strings=paste0("sl_q",1:4)
colnames(washp)[(ncol(washp)-4):(ncol(washp)-1)]=strings

washp=cbind(washp,pred_sigmasl)
strings=paste0("sigmasl_q",1:4)
colnames(washp)[(ncol(washp)-4):(ncol(washp)-1)]=strings

washp=cbind(washp,pred_disp)
strings=paste0("disp_q",1:4)
colnames(washp)[(ncol(washp)-4):(ncol(washp)-1)]=strings

washp=cbind(washp,pred_sigmadisp)
strings=paste0("sigmadisp_q",1:4)
colnames(washp)[(ncol(washp)-4):(ncol(washp)-1)]=strings

# Save wash with predictions --------------

saveRDS(washp,file.path(objdir,"wash_preds.rds"))

# Plot prediction maps ------------------------

#crop with usplot 
#this can take a while
washp2=st_intersection(washp,usplot)  

#test version,crop to AL
test=washp2[washp2$State_Name=="ALABAMA",]
string="test"
colname="sl_q1"
response="sl"
map=test
WashMap<-function(map,response,colname,incleg){
  #parse(text=colname)
  rci=grep(paste0("^",response,"_q"),colnames(map))
  rcv=c(as.matrix(st_drop_geometry(map[,rci])))
  maxp=max(rcv)
  minp=min(rcv)
  if(incleg){
  washmap=ggplot() + 
    geom_sf(data=usplot, fill="black")+
    geom_sf(data = map, aes(fill = eval(parse(text=colname))),lwd=0)+
    scale_fill_viridis_c(name=string,limits=c(minp,maxp))+
    geom_sf(data=usplot, fill=NA, color="#EBEBEB")+
    theme_map()+
    theme(plot.margin = unit(c(0.1,0.1,0.1,0.1), "in"))+
    theme(legend.position="bottom",
          legend.text = element_text(angle = 90))
  } else{
    washmap=ggplot() + 
      geom_sf(data=usplot, fill="black")+
      geom_sf(data = map, aes(fill = eval(parse(text=colname))),lwd=0)+
      scale_fill_viridis_c(name="",limits=c(minp,maxp))+
      geom_sf(data=usplot, fill=NA, color="#EBEBEB")+
      theme_map()+
      theme(plot.margin = unit(c(0.1,0.1,0.1,0.1), "in"))+
      theme(legend.position="none")
  }
  #washmap    
    

  return(washmap)
}

## SL -------------------------------------
#map,response,colname,incleg
washq1.sl=WashMap(washp2,"^sl","sl_q1",FALSE)
washq2.sl=WashMap(washp2,"^sl","sl_q2",FALSE)
washq3.sl=WashMap(washp2,"^sl","sl_q3",FALSE)
washq4.sl=WashMap(washp2,"^sl","sl_q4",FALSE)
washleg=WashMap(washp2,"^sl","sl_q1",TRUE)

sl_pg=plot_grid(washq1.sl, washq2.sl, washq3.sl, washq4.sl, nrow = 2, labels=c("a","b","c","d"))
ggsave(file.path(outdir,filestr,"FigTab","Maps","sl_pmap.png"),plot=sl_pg,height=6.5,width=9,units="in")
ggsave(file.path(outdir,filestr,"FigTab","Maps","sl_legend.png"),plot=washleg,height=6.5,width=9,units="in")

## Sigma SL -------------------------------------
#map,response,colname,incleg
washq1.sigma.sl=WashMap(washp2,"^sigmasl","sigmasl_q1",FALSE)
washq2.sigma.sl=WashMap(washp2,"^sigmasl","sigmasl_q2",FALSE)
washq3.sigma.sl=WashMap(washp2,"^sigmasl","sigmasl_q3",FALSE)
washq4.sigma.sl=WashMap(washp2,"^sigmasl","sigmasl_q4",FALSE)
washlegsigma.sl=WashMap(washp2,"^sigmasl","sigmasl_q1",TRUE)

sigmasl_pg=plot_grid(washq1.sigma.sl, washq2.sigma.sl, washq3.sigma.sl, washq4.sigma.sl, nrow = 2, labels=c("a","b","c","d"))
ggsave(file.path(outdir,filestr,"FigTab","Maps","sigma.sl_pmap.png"),plot=sigmasl_pg,height=6.5,width=9,units="in")
ggsave(file.path(outdir,filestr,"FigTab","Maps","sigma.sl_legend.png"),plot=washlegsigma.sl,height=6.5,width=9,units="in")


## Disp -------------------------------------
#map,response,colname,incleg
washq1.disp=WashMap(washp2,"^disp","disp_q1",FALSE)
washq2.disp=WashMap(washp2,"^disp","disp_q2",FALSE)
washq3.disp=WashMap(washp2,"^disp","disp_q3",FALSE)
washq4.disp=WashMap(washp2,"^disp","disp_q4",FALSE)
washleg.disp=WashMap(washp2,"^disp","disp_q1",TRUE)

disp_pg=plot_grid(washq1.disp, washq2.disp, washq3.disp, washq4.disp, nrow = 2, labels=c("a","b","c","d"))
ggsave(file.path(outdir,filestr,"FigTab","Maps","disp_pmap.png"),plot=disp_pg,height=6.5,width=9,units="in")
ggsave(file.path(outdir,filestr,"FigTab","Maps","disp_legend.png"),plot=washleg.disp,height=6.5,width=9,units="in")

## Sigma Disp -------------------------------------
#map,response,colname,incleg
washq1.sigma.disp=WashMap(washp2,"^sigmadisp","sigmadisp_q1",FALSE)
washq2.sigma.disp=WashMap(washp2,"^sigmadisp","sigmadisp_q2",FALSE)
washq3.sigma.disp=WashMap(washp2,"^sigmadisp","sigmadisp_q3",FALSE)
washq4.sigma.disp=WashMap(washp2,"^sigmadisp","sigmadisp_q4",FALSE)
washleg.sigma.disp=WashMap(washp2,"^sigmadisp","sigmadisp_q1",TRUE)

sigma.disp_pg=plot_grid(washq1.sigma.disp, washq2.sigma.disp, washq3.sigma.disp, washq4.sigma.disp, nrow = 2, labels=c("a","b","c","d"))
ggsave(file.path(outdir,filestr,"FigTab","Maps","sigma.disp_pmap.png"),plot=sigma.disp_pg,height=6.5,width=9,units="in")
ggsave(file.path(outdir,filestr,"FigTab","Maps","sigma.disp_legend.png"),plot=washleg.sigma.disp,height=6.5,width=9,units="in")

## Make map grids -------------------

sl_pg=plot_grid(washq1.sl, washq2.sl, washq3.sl, washq4.sl, nrow = 2)
sigmasl_pg=plot_grid(washq1.sigma.sl, washq2.sigma.sl, washq3.sigma.sl, washq4.sigma.sl, nrow = 2)
disp_pg=plot_grid(washq1.disp, washq2.disp, washq3.disp, washq4.disp, nrow = 2)
sigmadisp_pg=plot_grid(washq1.sigma.disp, washq2.sigma.disp, washq3.sigma.disp, washq4.sigma.disp, nrow = 2)

