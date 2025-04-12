
# Purpose ---------------------

#The purpose of this script is to make figures using the output from the GBM models

# Process outline ---------------------

#1-model selection table
  #R2 and RMSE from random/region k-fold CV methods
  #comparison of feature selection methods

#2-preds vs obs scatter/histograms from each top model
  #preds vs obs from top random k-fold
  #preds vs obs from top region model

# Setup ---------------------

library(ggforce)
library(ggplot2)
library(tidyverse)
library(dismo)
library(gbm)
library(cowplot)
library(pdp)
library(ggpmisc)
library(ggpubr)
library(grid)
library(gridExtra)
library(viridis)
library(gridGraphics)

#set directories
home="/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/2_Projects/StatPigMvmt/Pipeline_R2"
outdir=file.path(home,"4_Outputs")
funcdir<-file.path(home,"1_Scripts","MakeFigures","Functions")
filestr<-"04APR25_Runs"
objdir<-file.path(home,"2_Data","Objects")
if(!dir.exists(file.path(outdir,filestr,"FigTab"))){dir.create(file.path(outdir,filestr,"FigTab"))}
#read in pigsums dataset

pigsums<-readRDS(file.path(home,"2_Data","Objects","dailyPigSums.rds"))
pigswsite<-readRDS(file.path(home,"2_Data","Objects","geolocsnatl_wDispl.rds"))
pigs<-readRDS(file.path(home,"2_Data","Objects","geolocsnatl_wDispl.rds"))

## Format pigsums df ------------------
colnames(pigsums)
#region got duplicated, fix this
colnames(pigsums)[3]<-"region"

#make sure all correct classes
num.cols=c(9:14,17:64)
cat.cols=c(1,2,3,4,5,6,7,8)
  
pigsums[,num.cols] <- lapply(pigsums[,num.cols],as.numeric)
pigsums[,cat.cols] <- lapply(pigsums[,cat.cols],as.factor)

#Get vector of study/region names
studydf=unique(pigsums$region)

#make cutoff set for sigma models
pigsums_sigmasl=pigsums[pigsums$not_na_sl>30,]
pigsums_sigmadisp=pigsums[pigsums$not_na_displ>30,]

#Remove NAs to run GBMs later for diagnostic stats
pigsums_sl<-pigsums[!is.na(pigsums$sl_mean),]
pigsums_displ<-pigsums[!is.na(pigsums$displ_mean),]

#Set responses to integer for Poisson models
pigsums_sl$sl_mean<-as.integer(pigsums_sl$sl_mean)
pigsums_displ$displ_mean<-as.integer(pigsums_displ$displ_mean)


#Source functions
source(file.path(funcdir,"NeatenVarNames.R"))
source(file.path(funcdir,"MakePlotGrid.R"))
source(file.path(funcdir,"CombinePredObsFunction.R"))
source(file.path(funcdir,"GetStudyCVTables.R"))
source(file.path(funcdir,"multiplotfunction.R"))
source(file.path(funcdir,"MakePdpGridFunction.R"))
source(file.path(funcdir,"MakeVarSizeDotPlots.R"))
source(file.path(home,"1_Scripts","Analysis","GBM_nestedCV","GBM_Functions","K_Split.R"))


# Make model selection table ---------------------------

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
model.sel.tbl[1,5]=round(null.RMSE.sl[1,2],3)
model.sel.tbl[1,6]=round(null.R2.sl[1,2],3)
#row 2: sigma sl
model.sel.tbl[2,5]=round(null.RMSE.sigma.sl[1,2],3)
model.sel.tbl[2,6]=round(null.R2.sigma.sl[1,2],3)
#row 3: disp
model.sel.tbl[3,5]=round(null.RMSE.disp[1,2],3)
model.sel.tbl[3,6]=round(null.R2.disp[1,2],3)
#row 4: sigma disp
model.sel.tbl[4,5]=round(null.RMSE.sigma.disp[1,2],3)
model.sel.tbl[4,6]=round(null.R2.sigma.disp[1,2],3)

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
colnames(model.sel.tbl.total)=c("full RMSE random", "full r2 random", "drop 01 RMSE random", "drop 01 r2 random", "null RMSE random", "null R2 random","full RMSE region", "full r2 region", "drop 01 RMSE region", "drop 01 r2 region", "null RMSE region", "null R2region")

model.sel.tbl.total[1:4,1:6]=model.sel.tbl
model.sel.tbl.total[1:4,7:12]=model.sel.tbl.region

#Removing R2
model.sel.tbl.total=model.sel.tbl.total[,-c(2,4,6,8,10,12,14,16)]


write.csv(model.sel.tbl.total,file.path(outdir,filestr,"FigTab","ModelSelTblTotal.csv"))

# Run GBMs -----------------

getXvec<-function(model.sel.tbl.total,out.opt){
  X_vec_list=vector(mode="list",length=4)
  names(X_vec_list)<-c("Xsl","sigma.sl","Xdisp","sigma.disp")
  
  X_sel=data.frame(matrix(nrow=4,ncol=4))
  colnames(X_sel)<-c("response","X_vec","reg_ran","vars")

    for(i in 1:nrow(model.sel.tbl.total)){
      X_sel[i,1]=rownames(model.sel.tbl.total)[i]
      c_ind=which(model.sel.tbl.total[i,]==min(model.sel.tbl.total[i,]))
      X_sel[i,2]=colnames(model.sel.tbl.total)[c_ind[length(c_ind)]]
    }
    
    #parse text in X sel: region or random
    X_sel[grep("region",X_sel$X_vec),3]<-"Region"
    X_sel[grep("random",X_sel$X_vec),3]<-"Random"
    
    #parse text in X sel: variable sel
    X_sel[grep("full",X_sel$X_vec),4]<-"Full"
    X_sel[grep("drop 01",X_sel$X_vec),4]<-"X_vec_01Drop"
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
          which(colnames(pigsums)=="mean_rd3"),
          which(colnames(pigsums)=="mean_drt"):
            which(colnames(pigsums)=="var_rd3"),
          which(colnames(pigsums)=="var_drt"):
          which(colnames(pigsums)=="var_lc_24")
        )
      }
      if(X_sel[i,4]=="X_vec_01Drop"){
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


## get opt params for each top model --------

X_sel=getXvec(model.sel.tbl.total,"X_sel")

#hard coding for 01 drop for now
sl.opt.params.kfold=read.csv(file.path(outdir,filestr,X_sel[1,3],"GBM_01Drop","sl_bestmodelparams.csv"))
sigma.sl.opt.params.kfold=read.csv(file.path(outdir,filestr,X_sel[2,3],"sigma.sl_bestmodelparams.csv"))
disp.opt.params.kfold=read.csv(file.path(outdir,filestr,X_sel[3,3],"disp_bestmodelparams.csv"))
sigma.disp.opt.params.kfold=read.csv(file.path(outdir,filestr,X_sel[4,3],"GBM_01Drop","sigma.disp_bestmodelparams.csv"))

## determine X vec for each response --------
X_vec_list=getXvec(model.sel.tbl.total,"Xlist")
colnames(pigsums_displ)[X_vec_list$Xdisp]

## Run GBMs --------
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
gbm.sigma.sl=gbm.fixed(data=pigsums_sigmasl, gbm.x=X_vec_list$sigma.sl, gbm.y=which(colnames(pigsums)=="sl_disp"),
                       learning.rate=sigma.sl.opt.params.kfold$learning.rate, 
                       tree.complexity=sigma.sl.opt.params.kfold$tree.complexity, 
                       n.trees=sigma.sl.opt.params.kfold$n.trees,
                       bag.fraction=sigma.sl.opt.params.kfold$bag.fraction,
                       family="gaussian") 
gbm.sigma.disp=gbm.fixed(data=pigsums_sigmadisp, gbm.x=X_vec_list$sigma.disp, gbm.y=which(colnames(pigsums)=="displ_disp"),
                         learning.rate=sigma.disp.opt.params.kfold$learning.rate, 
                         tree.complexity=sigma.disp.opt.params.kfold$tree.complexity, 
                         n.trees=sigma.disp.opt.params.kfold$n.trees,
                         bag.fraction=sigma.disp.opt.params.kfold$bag.fraction,
                         family="gaussian") 

# Relative influence plot grids -----------------

#export importance tables/figures
rel.inf.sl=summary(gbm.sl)
rel.inf.sigma.sl=summary(gbm.sigma.sl)
rel.inf.disp=summary(gbm.disp)
rel.inf.sigma.disp=summary(gbm.sigma.disp)

#for top models, transform RI into percentages
rel.inf.sl$rel.prop=rel.inf.sl$rel.inf/rel.inf.sl$rel.inf[1]
rel.inf.sigma.sl$rel.prop=rel.inf.sigma.sl$rel.inf/rel.inf.sigma.sl$rel.inf[1]
rel.inf.disp$rel.prop=rel.inf.disp$rel.inf/rel.inf.disp$rel.inf[1]
rel.inf.sigma.disp$rel.prop=rel.inf.sigma.disp$rel.inf/rel.inf.sigma.disp$rel.inf[1]

rel.inf.sl=rel.inf.sl[order(rel.inf.sl[,1]),]
rel.inf.sigma.sl=rel.inf.sigma.sl[order(rel.inf.sigma.sl[,1]),]
rel.inf.disp=rel.inf.disp[order(rel.inf.disp[,1]),]
rel.inf.sigma.disp=rel.inf.sigma.disp[order(rel.inf.sigma.disp[,1]),]

colnames(rel.inf.sl)[3]="S.L."
colnames(rel.inf.sigma.sl)[3]="S.L. var"
colnames(rel.inf.disp)[3]="Displ."
colnames(rel.inf.sigma.disp)[3]="Displ. var"

#join all rel inf dfs together
join1=full_join(rel.inf.sl,rel.inf.sigma.sl,by="var")
join2=full_join(join1,rel.inf.disp,by="var")
join3=full_join(join2,rel.inf.sigma.disp,by="var")

#select needed cols
join4=join3[,c(1,3,5,7,9)]

#get long format
join5=as.data.frame(pivot_longer(join4,2:5))

#Neaten names
#all good, just want to replace lc nums with something more informative
#grep(join5[,1])
lc_vals=c(11,12,21,22,23,24,31,41,42,43,51,52,71,72,73,74,81,82,90,95)
lc_abbrev=c("water","ice","dev_open","dev_low","dev_med","dev_hi","barren","decid_for","ever_for","mix_for","dw_scrub","shrub","grass_herb","sedge_herb","lichen","moss","pasture","crop","wood_wetl","emherb_wetl")
lc_names=c("Water","Ice","Dev. Open","Dev. Low","Dev. Med.","Dev. High","Barren","Decid. For.","Evergr. For.","Mixed For.","Dw. Scrub","Shrub","Grassland","Sedge","Lichen","Moss","Pasture","Crop","Woody Wetl.","Em. Herb. Wetl.")
lc_df=data.frame("vals"=lc_vals,"abbrev"=lc_abbrev,"names"=lc_names)

Rename_LC=function(df,col,lc_df){
  for(i in 1:nrow(df)){
    if(length(grep("lc",df[i,col]))>0){
      lcnam=lc_df[which(lc_df$vals==stringr::str_sub(df[i,col],-2L)),3]
      if(length(grep("mean",df[i,col]))>0){
        df[i,col]=paste("mean",lcnam,sep=" ")
      }
      if(length(grep("var",df[i,col]))>0){
        df[i,col]=paste("var",lcnam,sep=" ")
      }
      
      }
  }
  return(df)
}

join5=Rename_LC(join5,1,lc_df)

#Make heat map with proportional rel infl.
hm_ri=ggplot(join5, aes(name, var, fill= value)) + 
  geom_tile() +
  scale_fill_viridis(discrete=FALSE) +
  theme_minimal()+
  labs(x ="Response", y = "Predictor", fill="Rel. Inf. Prop.")+
  theme(axis.title.y = element_text("Arial", "bold", size=12),
        axis.title.x = element_text("Arial", "bold", size=12),
        legend.title = element_text("Arial", "bold", size=12),
        axis.text=element_text(size=12),
        legend.text=element_text(size=12))

#save heatmap
ggsave(file.path(outdir,filestr,"FigTab","heatmap.png"),plot=hm_ri,height=9,width=6.5,units="in")

#remove variables with 0 influence from tables (full models)
RemoveZeroInfluence=function(relinftbl){
  relinftbl=relinftbl[(relinftbl$rel.inf!=0),]
  return(relinftbl)
}

rel.inf.sl=RemoveZeroInfluence(rel.inf.sl)
rel.inf.sigma.sl=RemoveZeroInfluence(rel.inf.sigma.sl)
rel.inf.disp=RemoveZeroInfluence(rel.inf.disp)
rel.inf.sigma.disp=RemoveZeroInfluence(rel.inf.sigma.disp)

#make formatted rel influence plots
ri.sl=ggplot(rel.inf.sl, aes(x = reorder(var,rel.inf), y = rel.inf, fill=rel.inf))+
  geom_bar(stat="identity",show.legend=FALSE)+
  coord_flip()+
  labs(x="Independent variables", y="Relative Influence", fill="Rel. Inf.")+
  scale_fill_viridis_c(begin=0.2,end=0.4)+
  theme_minimal()+
  theme(text = element_text(size = 28))

ri.sigmasl=ggplot(rel.inf.sigma.sl, aes(x = reorder(var,rel.inf), y = rel.inf, fill=rel.inf))+
  geom_bar(stat="identity",show.legend=FALSE)+
  coord_flip()+
  labs(x="Independent variables", y="Relative Influence", fill="Rel. Inf.")+
  scale_fill_viridis_c(begin=0.2,end=0.4)+
  theme_minimal()+
  theme(text = element_text(size = 28))

ri.disp=ggplot(rel.inf.disp, aes(x = reorder(var,rel.inf), y = rel.inf, fill=rel.inf))+
  geom_bar(stat="identity",show.legend=FALSE)+
  coord_flip()+
  labs(x="Independent variables", y="Relative Influence", fill="Rel. Inf.")+
  scale_fill_viridis_c(begin=0.3,end=0.4)+
  theme_minimal()+
  theme(text = element_text(size = 28))

ri.sigmadisp=ggplot(rel.inf.sigma.disp, aes(x = reorder(var,rel.inf), y = rel.inf, fill=rel.inf))+
  geom_bar(stat="identity",show.legend=FALSE)+
  coord_flip()+
  labs(x="Independent variables", y="Relative Influence", fill="Rel. Inf.")+
  scale_fill_viridis_c(begin=0.2,end=0.4)+
  theme_minimal()+
  theme(text = element_text(size = 28))

#make plotgrids
relinf.pg=plot_grid(
  ri.sl,
  ri.sigmasl,
  ri.disp,
  ri.sigmadisp,
  ncol = 2,
  labels=c("a","b","c","d"),
  label_size=50
)

ggsave(file.path(outdir,filestr,"FigTab","relinf_pg.png"),plot=relinf.pg,height=9,width=6.5,units="in")


# FIXPartial dependence plots -----------------

numplots=12

#Make partial dependence plots for each response top model, top 12
#sl
MakePdpGrid(numplots,rel.inf.sl,gbm.sl,pigsums_sl,sl.opt.params.kfold,c(-0.35,0.35))
#sigmasl
MakePdpGrid(numplots,rel.inf.sigma.sl,gbm.sigma.sl,pigsums_sigmasl,sigma.sl.opt.params.kfold,c(-0.35,0.35))
#displ
MakePdpGrid(numplots,rel.inf.disp,gbm.disp,pigsums_displ,disp.opt.params.kfold,c(-0.9,0.9))
#sigmadispl
MakePdpGrid(numplots,rel.inf.sigma.disp,gbm.sigma.disp,pigsums_sigmadisp,sigma.disp.opt.params.kfold,c(-0.15,0.15))
#can't save as plot objects, need to save manually

# Cross-validation stats, region vs. random -----------------

#Get RMSE and R2 values for out of site/sample predictions
#studydf<-data.frame("region"=studydf)
#studydf=studydf$region

## Get basic CV tables -----------------

CVstats_sl.random=GetCVStats_Table(pigsums_sl,X_vec_list$Xsl,"sl_mean",sl.opt.params.kfold,"poisson","random",studydf,"CVtbl_only")
CVstats_sigma.sl.random=GetCVStats_Table(pigsums_sigmasl,X_vec_list$sigma.sl,"sl_disp",sigma.sl.opt.params.kfold,"gaussian","random",studydf,"CVtbl_only")
CVstats_disp.random=GetCVStats_Table(pigsums_displ,X_vec_list$Xdisp,"displ_mean",disp.opt.params.kfold,"poisson","random",studydf,"CVtbl_only")
CVstats_sigma.disp.random=GetCVStats_Table(pigsums_sigmadisp,X_vec_list$sigma.disp,"displ_disp",sigma.disp.opt.params.kfold,"gaussian","random",studydf,"CVtbl_only")

#save CVstats tables-- take a while to run
saveRDS(CVstats_sl.random,file.path(objdir,"CVstats_sl.random.rds"))
saveRDS(CVstats_sigma.sl.random,file.path(objdir,"CVstats_sigma.sl.random.rds"))
saveRDS(CVstats_disp.random,file.path(objdir,"CVstats_disp.random.rds"))
saveRDS(CVstats_sigma.disp.random,file.path(objdir,"CVstats_sigma.disp.random.rds"))

## order rows by state -----------------
#Make study counts to var dot plot sizes
studycounts=as.data.frame(pigsums_sl %>% 
                            group_by(region) %>% 
                            dplyr::summarise(numpigdays=sum(not_na_sl), 
                                             numpigs=n_distinct(animalid),
                                             State=first(state)))

## Get study IDs instead of names -----------------
colnames(CVstats_sl.random)[1]<-"region"
colnames(CVstats_sigma.sl.random)[1]<-"region"
colnames(CVstats_disp.random)[1]<-"region"
colnames(CVstats_sigma.disp.random)[1]<-"region"

snk=readRDS(file.path(objdir,"studynumkey.rds"))
colnames(snk)[1]<-"region"

CVstats_sl.random=left_join(CVstats_sl.random,snk)
CVstats_sigma.sl.random=left_join(CVstats_sigma.sl.random,snk)
CVstats_disp.random=left_join(CVstats_disp.random,snk)
CVstats_sigma.disp.random=left_join(CVstats_sigma.disp.random,snk)

## Make CV dot plots -----------------

#Draw dotplots
sl.rmse.dot=ggdraw(make.varsize.dotplots(CVstats_sl.random,studycounts,"RMSE","sl"))
sigmasl.rmse.dot=ggdraw(make.varsize.dotplots(CVstats_sigma.sl.random,studycounts, "RMSE","sigmasl"))
disp.rmse.dot=ggdraw(make.varsize.dotplots(CVstats_disp.random,studycounts, "RMSE","disp"))
sigmadisp.rmse.dot=ggdraw(make.varsize.dotplots(CVstats_sigma.disp.random,studycounts, "RMSE","sigmadisp"))

sl.r2.dot=ggdraw(make.varsize.dotplots(CVstats_sl.random,studycounts,"R2","sl"))
sigmasl.r2.dot=ggdraw(make.varsize.dotplots(CVstats_sigma.sl.random,studycounts, "R2","sigmasl"))
disp.r2.dot=ggdraw(make.varsize.dotplots(CVstats_disp.random,studycounts, "R2","disp"))
sigmadisp.r2.dot=ggdraw(make.varsize.dotplots(CVstats_sigma.disp.random,studycounts, "R2","sigmadisp"))

#save dotplots
path=file.path(outdir,filestr,"FigTab","dotplots_RMSE_R2")
if(!dir.exists(path)){dir.create(path)}

ggsave(file.path(path,"sl.rmse.dot.png"),plot=sl.rmse.dot,height=18,width=12,units="in")
ggsave(file.path(path,"sigmasl.rmse.dot.png"),plot=sigmasl.rmse.dot,height=18,width=12,units="in")
ggsave(file.path(path,"disp.rmse.dot.png"),plot=disp.rmse.dot,height=18,width=12,units="in")
ggsave(file.path(path,"sigmadisp.rmse.dot.png"),plot=sigmadisp.rmse.dot,height=18,width=12,units="in")

ggsave(file.path(path,"sl.r2.dot.png"),plot=sl.r2.dot,height=18,width=12,units="in")
ggsave(file.path(path,"sigmasl.r2.dot.png"),plot=sigmasl.r2.dot,height=18,width=12,units="in")
ggsave(file.path(path,"disp.r2.dot.png"),plot=disp.r2.dot,height=18,width=12,units="in")
ggsave(file.path(path,"sigmadisp.r2.dot.png"),plot=sigmadisp.r2.dot,height=18,width=12,units="in")

# Pred vs obs plots -----------------
CVstats_sl.random=GetCVStats_Table(pigsums_sl,X_vec_list$Xsl,"sl_mean",sl.opt.params.kfold,"poisson","random",studydf,"all")
CVstats_sigma.sl.random=GetCVStats_Table(pigsums_sigmasl,X_vec_list$sigma.sl,"sl_disp",sigma.sl.opt.params.kfold,"gaussian","random",studydf,"all")
CVstats_disp.random=GetCVStats_Table(pigsums_displ,X_vec_list$Xdisp,"displ_mean",disp.opt.params.kfold,"poisson","random",studydf,"all")
CVstats_sigma.disp.random=GetCVStats_Table(pigsums_sigmadisp,X_vec_list$sigma.disp,"displ_disp",sigma.disp.opt.params.kfold,"gaussian","random",studydf,"all")

saveRDS(CVstats_sl.random,file.path(objdir,"CVstats_sl.random_list.rds"))
saveRDS(CVstats_sigma.sl.random,file.path(objdir,"CVstats_sigma.sl.random_list.rds"))
saveRDS(CVstats_disp.random,file.path(objdir,"CVstats_disp.random_list.rds"))
saveRDS(CVstats_sigma.disp.random,file.path(objdir,"CVstats_sigma.disp.random_list.rds"))

#format prediction/test sets into dataframe for plotting
sl_predobs.df=CombinePredObs(CVstats_sl.random)
sigma.sl_predobs.df=CombinePredObs(CVstats_sigma.sl.random)
disp_predobs.df=CombinePredObs(CVstats_disp.random)
sigma.disp_predobs.df=CombinePredObs(CVstats_sigma.disp.random)

#format for plotting histograms
sl_predobs.df2=as.data.frame(pivot_longer(sl_predobs.df,cols=c("test","preds")))
sigma.sl_predobs.df2=as.data.frame(pivot_longer(sigma.sl_predobs.df,cols=c("test","preds")))
disp_predobs.df2=as.data.frame(pivot_longer(disp_predobs.df,cols=c("test","preds")))
sigma.disp_predobs.df2=as.data.frame(pivot_longer(sigma.disp_predobs.df,cols=c("test","preds")))

#get scatter plots
sl.po.scat=ggplot(sl_predobs.df, aes(test, preds, color = factor(CVmethod))) + 
  geom_point(aes(color=factor(CVmethod)), size=0.5, alpha=1)+
  geom_abline(intercept = 0, slope = 1,alpha=0.2)+
  xlim(0,3000)+
  ylim(0,3000)+
  theme_minimal()+
  facet_wrap(vars(region))+
  geom_mark_ellipse(aes(color = factor(CVmethod)))+
  scale_color_manual(values=c("#039e8c", "#832b9e"), name="CV method")+
  labs(x="observed", y="predicted")+
  theme(text = element_text(size = 12))+
  theme(axis.text.x = element_text(angle = 90, size=10),
        axis.text.y = element_text(size=9))

sigmasl.po.scat=ggplot(sigma.sl_predobs.df, aes(test, preds, color = factor(CVmethod))) + 
  geom_point(aes(color=factor(CVmethod)), size=0.5, alpha=1)+
  geom_abline(intercept = 0, slope = 1,alpha=0.2)+
  xlim(0,3)+
  ylim(0,3)+
  theme_minimal()+
  facet_wrap(vars(region))+
  geom_mark_ellipse(aes(color = factor(CVmethod)))+
  scale_color_manual(values=c("#039e8c", "#832b9e"), name="CV method")+
  labs(x="observed", y="predicted")+
  theme(text = element_text(size = 12))+
  theme(axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=9))

disp.po.scat=ggplot(disp_predobs.df, aes(test, preds, color = factor(region))) + 
  geom_point(aes(color=factor(CVmethod)), size=0.5, alpha=1)+
  geom_abline(intercept = 0, slope = 1)+
  xlim(0,12500)+
  ylim(0,12500)+
  theme_minimal()+
  facet_wrap(vars(region))+
  geom_mark_ellipse(aes(color = factor(CVmethod)))+
  scale_color_manual(values=c("#039e8c", "#832b9e"), name="CV method")+
  labs(x="observed", y="predicted")+
  theme(text = element_text(size = 12))+
  theme(axis.text.x = element_text(angle = 90, size=10),
        axis.text.y = element_text(size=9))

sigmadisp.po.scat=ggplot(sigma.disp_predobs.df, aes(test, preds, color = factor(region))) + 
  geom_point(aes(color=factor(CVmethod)), size=0.5, alpha=1)+
  geom_abline(intercept = 0, slope = 1)+
  xlim(0,3)+
  ylim(0,3)+
  theme_minimal()+
  facet_wrap(vars(region))+
  geom_mark_ellipse(aes(color = factor(CVmethod)))+
  scale_color_manual(values=c("#039e8c", "#832b9e"), name="CV method")+
  labs(x="observed", y="predicted")+
  theme(text = element_text(size = 12))+
  theme(axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=9))

#save po scats
path=file.path(outdir,filestr,"FigTab","pred_obs_scatters")
if(!dir.exists(path)){dir.create(path)}

ggsave(file.path(path,"sl_po_scatter.png"),plot=sl.po.scat,height=18,width=12,units="in")
ggsave(file.path(path,"sigmasl_po_scatter.png"),plot=sigmasl.po.scat,height=18,width=12,units="in")
ggsave(file.path(path,"disp_po_scatter.png"),plot=disp.po.scat,height=18,width=12,units="in")
ggsave(file.path(path,"sigmadisp_po_scatter.png"),plot=sigmadisp.po.scat,height=18,width=12,units="in")

# Density plots -----------------

#Formatting
sl_predobs.df2[sl_predobs.df2$name=="test",3]<-"observed"
sl_predobs.df2[sl_predobs.df2$name=="preds",3]<-"predicted"
sigma.sl_predobs.df2[sigma.sl_predobs.df2$name=="test",3]<-"observed"
sigma.sl_predobs.df2[sigma.sl_predobs.df2$name=="preds",3]<-"predicted"
disp_predobs.df2[disp_predobs.df2$name=="test",3]<-"observed"
disp_predobs.df2[disp_predobs.df2$name=="preds",3]<-"predicted"
sigma.disp_predobs.df2[sigma.disp_predobs.df2$name=="test",3]<-"observed"
sigma.disp_predobs.df2[sigma.disp_predobs.df2$name=="preds",3]<-"predicted"

sl.dens<-sl_predobs.df2%>%
  ggplot(aes(x=value, fill=name)) +
  geom_density(alpha=0.5)+ 
  labs(x= "Daily step length mean")+
  theme_minimal()+
  scale_fill_manual(values=c("#25858e", "#2ab07f"), name="CV method")+
  facet_wrap(vars(CVmethod))+
  theme(text = element_text(size = 40))

sigmasl.dens<-sigma.sl_predobs.df2%>%
  ggplot(aes(x=value, fill=name)) +
  geom_density(alpha=0.5)+ 
  labs(x= "Daily step length dispersion")+
  theme_minimal()+
  scale_fill_manual(values=c("#25858e", "#2ab07f"), name="CV method")+
  facet_wrap(vars(CVmethod))+
  theme(text = element_text(size = 40))

disp.dens<-disp_predobs.df2%>%
  ggplot(aes(x=value, fill=name)) +
  geom_density(alpha=0.5)+ 
  labs(x= "Daily displacement mean")+
  theme_minimal()+
  scale_fill_manual(values=c("#25858e", "#2ab07f"), name="CV method")+
  facet_wrap(vars(CVmethod))+
  theme(text = element_text(size = 40))

sigmadisp.dens<-sigma.disp_predobs.df2%>%
  ggplot(aes(x=value, fill=name)) +
  geom_density(alpha=0.5)+ 
  labs(x= "Daily displacement dispersion")+
  theme_minimal()+
  scale_fill_manual(values=c("#25858e", "#2ab07f"), name="CV method")+
  facet_wrap(vars(CVmethod))+
  theme(text = element_text(size = 40))

#save histos
path=file.path(outdir,filestr,"FigTab","histograms")
if(!dir.exists(path)){dir.create(path)}

ggsave(file.path(path,"sl_dens.png"),plot=sl.dens,height=18,width=12,units="in")
ggsave(file.path(path,"sigsl_dens.png"),plot=sigmasl.dens,height=18,width=12,units="in")
ggsave(file.path(path,"disp_dens.png"),plot=disp.dens,height=18,width=12,units="in")
ggsave(file.path(path,"sigdisp_dens.png"),plot=sigmadisp.dens,height=18,width=12,units="in")


# Violin plots -----------------------------

p1=ggplot(pigsums_sl, aes(x=fct_inorder(region), y=sl_mean)) + 
  geom_violin() + 
  coord_flip() + 
  geom_jitter(height=0,width=0.2,alpha=0.5) + theme(text = element_text(size = 40))+
  labs(x="study", y="average daily step length (m)")

p2=ggplot(pigsums_sigmasl, aes(x=fct_inorder(region), y=sl_disp)) + 
  geom_violin() + 
  coord_flip() + 
  geom_jitter(height=0,width=0.2,alpha=0.5) + theme(text = element_text(size = 40))+
  labs(x="study", y="average daily step length dispersion")

p3=ggplot(pigsums_displ, aes(x=fct_inorder(region), y=displ_mean)) + 
  geom_violin() + 
  coord_flip() + 
  geom_jitter(height=0,width=0.2,alpha=0.5) + theme(text = element_text(size = 40))+
  labs(x="study", y="average daily displacement (m)")

p4=ggplot(pigsums_sigmadisp, aes(x=fct_inorder(region), y=displ_disp)) + 
  geom_violin() + 
  coord_flip() + 
  geom_jitter(height=0,width=0.2,alpha=0.5) + theme(text = element_text(size = 40))+
  labs(x="study", y="average daily displacement dispersion")

path=file.path(outdir,filestr,"FigTab","violin_plots")
if(!dir.exists(path)){dir.create(path)}

png(file=file.path(path,"violin_grid.png"),
    width=1200, height=4000)
ggarrange(p1,p2,p3,p4,nrow=4,ncol=1,labels="auto",font.label=list(size=50,face="bold"),hjust=-0.2, vjust=1)
dev.off()



#Write any final output needed as objects ---------
saveRDS(X_sel,file.path(objdir,"X_sel.rds"))
saveRDS(X_vec_list,file.path(objdir,"X_vec_list.rds"))


