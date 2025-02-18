
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
home="/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/StatPigMvmt/Pipeline_R2"
outdir=file.path(home,"4_Outputs")

#read in pigsums dataset
pigsums<-readRDS(file.path(home,"2_Data","Objects","dailyPigSums.rds"))
pigswsite<-readRDS(file.path(home,"2_Data","Objects","geolocsnatl_wDispl.rds"))
pigs<-readRDS(file.path(home,"2_Data","Objects","geolocsnatl_wDispl.rds"))

#region got duplicated, fix this
pigsums2<-pigsums2[,-1]
pigsums2<-pigsums2[,-c(2)]
colnames(pigsums2)[ncol(pigsums2)]<-"region"

#make sure all correct classes
pigsums2[,c(2,3,6:48)] <- lapply(pigsums2[,c(2,3,6:48)],as.numeric)
pigsums2[,c(1,4,5,49)] <- lapply(pigsums2[,c(1,4,5,49)],as.factor)

pigsums=pigsums2
studydf=unique(pigsites[,c(1,2)])

#make cutoff set for sigma models
cutoff_sl=pigsums2[pigsums2$count>150,]

#Source functions
source(paste0(home,"/Scripts_Polished/Making_Plots_Tables/PlotFunctions/NeatenVarNames.R"))
source(paste0(home,"/Scripts_Polished/Making_Plots_Tables/PlotFunctions/MakePlotGrid.R"))
source(paste0(home,"/Scripts_Polished/Making_Plots_Tables/PlotFunctions/CombinePredObsFunction.R"))
source(paste0(home,"/Scripts_Polished/Making_Plots_Tables/PlotFunctions/GetStudyCVTables.R"))
source(paste0(home,"/Scripts_Polished/Making_Plots_Tables/PlotFunctions/multiplotfunction.R"))
source(paste0(home,"/Scripts_Polished/Making_Plots_Tables/PlotFunctions/MakePdpGridFunction.R"))
source(paste0(home,"/Scripts_Polished/Making_Plots_Tables/PlotFunctions/MakeVarSizeDotPlots.R"))
source(paste0(home,"/Scripts_Polished/GBM_nestedCV/K_Split.R"))

#########################
##Model selection table##
#########################

#make table with all mean R2 and RMSEs
#row for each response (5)
#col for each 2 metrics (RMSE and R2) and 4 models (full, drop01, lassodrop, null)

#initiate tables, region and random cv methods
model.sel.tbl=as.data.frame(matrix(nrow=5,ncol=8))
rownames(model.sel.tbl)=c("step length", "step length sigma", "displacement", "disp. sigma", "top average sl")
colnames(model.sel.tbl)=c("full RMSE", "full r2", "drop 01 RMSE", "drop 01 r2", "lasso RMSE", "lasso r2", "null RMSE", "null R2")
model.sel.tbl.region=as.data.frame(matrix(nrow=5,ncol=8))
rownames(model.sel.tbl.region)=c("step length", "step length sigma", "displacement", "disp. sigma", "top average sl")
colnames(model.sel.tbl.region)=c("full RMSE", "full r2", "drop 01 RMSE", "drop 01 r2", "lasso RMSE", "lasso r2", "null RMSE", "null R2")

#Full model: get values from file
full.RMSE.sl=read.csv(paste0(outdir,"05APR23_Runs/Random/sl_meanRMSE.csv"))
full.R2.sl=read.csv(paste0(outdir,"05APR23_Runs/Random/sl_meanR2.csv"))
full.RMSE.sigma.sl=read.csv(paste0(outdir,"05APR23_Runs/Random/sigma.sl_meanRMSE.csv"))
full.R2.sigma.sl=read.csv(paste0(outdir,"05APR23_Runs/Random/sigma.sl_meanR2.csv"))
full.RMSE.disp=read.csv(paste0(outdir,"05APR23_Runs/Random/disp_meanRMSE.csv"))
full.R2.disp=read.csv(paste0(outdir,"05APR23_Runs/Random/disp_meanR2.csv"))
full.RMSE.sigma.disp=read.csv(paste0(outdir,"05APR23_Runs/Random/sigma.disp_meanRMSE.csv"))
full.R2.sigma.disp=read.csv(paste0(outdir,"05APR23_Runs/Random/sigma.disp_meanR2.csv"))
full.RMSE.tenavg=read.csv(paste0(outdir,"05APR23_Runs/Random/tenavg_meanRMSE.csv"))
full.R2.tenavg=read.csv(paste0(outdir,"05APR23_Runs/Random/tenavg_meanR2.csv"))

#Full model region kfold: get values from file
full.RMSE.sl.region=read.csv(paste0(outdir,"05APR23_Runs/Region/sl_meanRMSE.csv"))
full.R2.sl.region=read.csv(paste0(outdir,"05APR23_Runs/Region/sl_meanR2.csv"))
full.RMSE.sigma.sl.region=read.csv(paste0(outdir,"05APR23_Runs/Region/sigma.sl_meanRMSE.csv"))
full.R2.sigma.sl.region=read.csv(paste0(outdir,"05APR23_Runs/Region/sigma.sl_meanR2.csv"))
full.RMSE.disp.region=read.csv(paste0(outdir,"05APR23_Runs/Region/disp_meanRMSE.csv"))
full.R2.disp.region=read.csv(paste0(outdir,"05APR23_Runs/Region/disp_meanR2.csv"))
full.RMSE.sigma.disp.region=read.csv(paste0(outdir,"05APR23_Runs/Region/sigma.disp_meanRMSE.csv"))
full.R2.sigma.disp.region=read.csv(paste0(outdir,"05APR23_Runs/Region/sigma.disp_meanR2.csv"))
full.RMSE.tenavg.region=read.csv(paste0(outdir,"05APR23_Runs/Region/tenavg_meanRMSE.csv"))
full.R2.tenavg.region=read.csv(paste0(outdir,"05APR23_Runs/Region/tenavg_meanR2.csv"))

#Drop01 model: get vals
drop01.RMSE.sl=read.csv(paste0(outdir,"05APR23_Runs/Random/GBM_01Drop/sl_meanRMSE.csv"))
drop01.R2.sl=read.csv(paste0(outdir,"05APR23_Runs/Random/GBM_01Drop/sl_meanR2.csv"))
drop01.RMSE.sigma.sl=read.csv(paste0(outdir,"05APR23_Runs/Random/GBM_01Drop/sigma.sl_meanRMSE.csv"))
drop01.R2.sigma.sl=read.csv(paste0(outdir,"05APR23_Runs/Random/GBM_01Drop/sigma.sl_meanR2.csv"))
drop01.RMSE.disp=read.csv(paste0(outdir,"05APR23_Runs/Random/GBM_01Drop/disp_meanRMSE.csv"))
drop01.R2.disp=read.csv(paste0(outdir,"05APR23_Runs/Random/GBM_01Drop/disp_meanR2.csv"))
drop01.RMSE.sigma.disp=read.csv(paste0(outdir,"05APR23_Runs/Random/GBM_01Drop/sigma.disp_meanRMSE.csv"))
drop01.R2.sigma.disp=read.csv(paste0(outdir,"05APR23_Runs/Random/GBM_01Drop/sigma.disp_meanR2.csv"))
drop01.RMSE.tenavg=read.csv(paste0(outdir,"05APR23_Runs/Random/GBM_01Drop/tenavg_meanRMSE.csv"))
drop01.R2.tenavg=read.csv(paste0(outdir,"05APR23_Runs/Random/GBM_01Drop/tenavg_meanR2.csv"))

#Drop01 region kfold model: get vals
drop01.RMSE.sl.region=read.csv(paste0(outdir,"05APR23_Runs/Region/GBM_01Drop/sl_meanRMSE.csv"))
drop01.R2.sl.region=read.csv(paste0(outdir,"05APR23_Runs/Region/GBM_01Drop/sl_meanR2.csv"))
drop01.RMSE.sigma.sl.region=read.csv(paste0(outdir,"05APR23_Runs/Region/GBM_01Drop/sigma.sl_meanRMSE.csv"))
drop01.R2.sigma.sl.region=read.csv(paste0(outdir,"05APR23_Runs/Region/GBM_01Drop/sigma.sl_meanR2.csv"))
drop01.RMSE.disp.region=read.csv(paste0(outdir,"05APR23_Runs/Region/GBM_01Drop/disp_meanRMSE.csv"))
drop01.R2.disp.region=read.csv(paste0(outdir,"05APR23_Runs/Region/GBM_01Drop/disp_meanR2.csv"))
drop01.RMSE.sigma.disp.region=read.csv(paste0(outdir,"05APR23_Runs/Region/GBM_01Drop/sigma.disp_meanRMSE.csv"))
drop01.R2.sigma.disp.region=read.csv(paste0(outdir,"05APR23_Runs/Region/GBM_01Drop/sigma.disp_meanR2.csv"))
drop01.RMSE.tenavg.region=read.csv(paste0(outdir,"05APR23_Runs/Region/GBM_01Drop/tenavg_meanRMSE.csv"))
drop01.R2.tenavg.region=read.csv(paste0(outdir,"05APR23_Runs/Region/GBM_01Drop/tenavg_meanR2.csv"))

#lasso model: get vals
lasso.RMSE.sl=read.csv(paste0(outdir,"05APR23_Runs/Random/GBM_Lasso/sl_meanRMSE.csv"))
lasso.R2.sl=read.csv(paste0(outdir,"05APR23_Runs/Random/GBM_Lasso/sl_meanR2.csv"))
lasso.RMSE.sigma.sl=read.csv(paste0(outdir,"05APR23_Runs/Random/GBM_Lasso/sigma.sl_meanRMSE.csv"))
lasso.R2.sigma.sl=read.csv(paste0(outdir,"05APR23_Runs/Random/GBM_Lasso/sigma.sl_meanR2.csv"))
lasso.RMSE.disp=read.csv(paste0(outdir,"05APR23_Runs/Random/GBM_Lasso/disp_meanRMSE.csv"))
lasso.R2.disp=read.csv(paste0(outdir,"05APR23_Runs/Random/GBM_Lasso/disp_meanR2.csv"))
lasso.RMSE.sigma.disp=read.csv(paste0(outdir,"05APR23_Runs/Random/GBM_Lasso/sigma.disp_meanRMSE.csv"))
lasso.R2.sigma.disp=read.csv(paste0(outdir,"05APR23_Runs/Random/GBM_Lasso/sigma.disp_meanR2.csv"))
lasso.RMSE.tenavg=read.csv(paste0(outdir,"05APR23_Runs/Random/GBM_Lasso/tenavg_meanRMSE.csv"))
lasso.R2.tenavg=read.csv(paste0(outdir,"05APR23_Runs/Random/GBM_Lasso/tenavg_meanR2.csv"))

#lasso model region kfold: get vals
lasso.RMSE.sl.region=read.csv(paste0(outdir,"05APR23_Runs/Region/GBM_Lasso/sl_meanRMSE.csv"))
lasso.R2.sl.region=read.csv(paste0(outdir,"05APR23_Runs/Region/GBM_Lasso/sl_meanR2.csv"))
lasso.RMSE.sigma.sl.region=read.csv(paste0(outdir,"05APR23_Runs/Region/GBM_Lasso/sigma.sl_meanRMSE.csv"))
lasso.R2.sigma.sl.region=read.csv(paste0(outdir,"05APR23_Runs/Region/GBM_Lasso/sigma.sl_meanR2.csv"))
lasso.RMSE.disp.region=read.csv(paste0(outdir,"05APR23_Runs/Region/GBM_Lasso/disp_meanRMSE.csv"))
lasso.R2.disp.region=read.csv(paste0(outdir,"05APR23_Runs/Region/GBM_Lasso/disp_meanR2.csv"))
lasso.RMSE.sigma.disp.region=read.csv(paste0(outdir,"05APR23_Runs/Region/GBM_Lasso/sigma.disp_meanRMSE.csv"))
lasso.R2.sigma.disp.region=read.csv(paste0(outdir,"05APR23_Runs/Region/GBM_Lasso/sigma.disp_meanR2.csv"))
lasso.RMSE.tenavg.region=read.csv(paste0(outdir,"05APR23_Runs/Region/GBM_Lasso/tenavg_meanRMSE.csv"))
lasso.R2.tenavg.region=read.csv(paste0(outdir,"05APR23_Runs/Region/GBM_Lasso/tenavg_meanR2.csv"))

#null model: get vals
null.RMSE.sl=read.csv(paste0(outdir,"05APR23_Runs/Random/GBM_Null/sl_meanRMSE.csv"))
null.R2.sl=read.csv(paste0(outdir,"05APR23_Runs/Random/GBM_Null/sl_meanR2.csv"))
null.RMSE.sigma.sl=read.csv(paste0(outdir,"05APR23_Runs/Random/GBM_Null/sigma.sl_meanRMSE.csv"))
null.R2.sigma.sl=read.csv(paste0(outdir,"05APR23_Runs/Random/GBM_Null/sigma.sl_meanR2.csv"))
null.RMSE.disp=read.csv(paste0(outdir,"05APR23_Runs/Random/GBM_Null/disp_meanRMSE.csv"))
null.R2.disp=read.csv(paste0(outdir,"05APR23_Runs/Random/GBM_Null/disp_meanR2.csv"))
null.RMSE.sigma.disp=read.csv(paste0(outdir,"05APR23_Runs/Random/GBM_Null/sigma.disp_meanRMSE.csv"))
null.R2.sigma.disp=read.csv(paste0(outdir,"05APR23_Runs/Random/GBM_Null/sigma.disp_meanR2.csv"))
null.RMSE.tenavg=read.csv(paste0(outdir,"05APR23_Runs/Random/GBM_Null/tenavg_meanRMSE.csv"))
null.R2.tenavg=read.csv(paste0(outdir,"05APR23_Runs/Random/GBM_Null/tenavg_meanR2.csv"))

#null model region: get vals
null.RMSE.sl.region=read.csv(paste0(outdir,"05APR23_Runs/Region/GBM_Null/sl_meanRMSE.csv"))
null.R2.sl.region=read.csv(paste0(outdir,"05APR23_Runs/Region/GBM_Null/sl_meanR2.csv"))
null.RMSE.sigma.sl.region=read.csv(paste0(outdir,"05APR23_Runs/Region/GBM_Null/sigma.sl_meanRMSE.csv"))
null.R2.sigma.sl.region=read.csv(paste0(outdir,"05APR23_Runs/Region/GBM_Null/sigma.sl_meanR2.csv"))
null.RMSE.disp.region=read.csv(paste0(outdir,"05APR23_Runs/Region/GBM_Null/disp_meanRMSE.csv"))
null.R2.disp.region=read.csv(paste0(outdir,"05APR23_Runs/Region/GBM_Null/disp_meanR2.csv"))
null.RMSE.sigma.disp.region=read.csv(paste0(outdir,"05APR23_Runs/Region/GBM_Null/sigma.disp_meanRMSE.csv"))
null.R2.sigma.disp.region=read.csv(paste0(outdir,"05APR23_Runs/Region/GBM_Null/sigma.disp_meanR2.csv"))
null.RMSE.tenavg.region=read.csv(paste0(outdir,"05APR23_Runs/Region/GBM_Null/tenavg_meanRMSE.csv"))
null.R2.tenavg.region=read.csv(paste0(outdir,"05APR23_Runs/Region/GBM_Null/tenavg_meanR2.csv"))

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
#row 5: sigma disp
model.sel.tbl[5,1]=round(full.RMSE.tenavg[1,2],3)
model.sel.tbl[5,2]=round(full.R2.tenavg[1,2],3)

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
#row 5: sigma disp
model.sel.tbl[5,3]=round(drop01.RMSE.tenavg[1,2],3)
model.sel.tbl[5,4]=round(drop01.R2.tenavg[1,2],3)

#col 5,6: lasso drop
#row 1: sl
model.sel.tbl[1,5]=round(lasso.RMSE.sl[1,2],3)
model.sel.tbl[1,6]=round(lasso.R2.sl[1,2],3)
#row 2: sigma sl
model.sel.tbl[2,5]=round(lasso.RMSE.sigma.sl[1,2],3)
model.sel.tbl[2,6]=round(lasso.R2.sigma.sl[1,2],3)
#row 3: disp
model.sel.tbl[3,5]=round(lasso.RMSE.disp[1,2],3)
model.sel.tbl[3,6]=round(lasso.R2.disp[1,2],3)
#row 4: sigma disp
model.sel.tbl[4,5]=round(lasso.RMSE.sigma.disp[1,2],3)
model.sel.tbl[4,6]=round(lasso.R2.sigma.disp[1,2],3)
#row 5: sigma disp
model.sel.tbl[5,5]=round(lasso.RMSE.tenavg[1,2],3)
model.sel.tbl[5,6]=round(lasso.R2.tenavg[1,2],3)

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
#row 5: sigma disp
model.sel.tbl[5,7]=round(null.RMSE.tenavg[1,2],3)
model.sel.tbl[5,8]=round(null.R2.tenavg[1,2],3)

#write.csv(model.sel.tbl,paste(outdir,"ModelSelTbl_14DEC22.csv"))

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
#row 5: sigma disp
model.sel.tbl.region[5,1]=round(full.RMSE.tenavg.region[1,2],3)
model.sel.tbl.region[5,2]=round(full.R2.tenavg.region[1,2],3)

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
#row 5: sigma disp
model.sel.tbl.region[5,3]=round(drop01.RMSE.tenavg.region[1,2],3)
model.sel.tbl.region[5,4]=round(drop01.R2.tenavg.region[1,2],3)

#col 5,6: region lasso model rmse, r2
#row 1: sl
model.sel.tbl.region[1,5]=round(lasso.RMSE.sl.region[1,2],3)
model.sel.tbl.region[1,6]=round(lasso.R2.sl.region[1,2],3)
#row 2: sigma sl
model.sel.tbl.region[2,5]=round(lasso.RMSE.sigma.sl.region[1,2],3)
model.sel.tbl.region[2,6]=round(lasso.R2.sigma.sl.region[1,2],3)
#row 3: disp
model.sel.tbl.region[3,5]=round(lasso.RMSE.disp.region[1,2],3)
model.sel.tbl.region[3,6]=round(lasso.R2.disp.region[1,2],3)
#row 4: sigma disp
model.sel.tbl.region[4,5]=round(lasso.RMSE.sigma.disp.region[1,2],3)
model.sel.tbl.region[4,6]=round(lasso.R2.sigma.disp.region[1,2],3)
#row 5: sigma disp
model.sel.tbl.region[5,5]=round(lasso.RMSE.tenavg.region[1,2],3)
model.sel.tbl.region[5,6]=round(lasso.R2.tenavg.region[1,2],3)

#col 7,8: region null model rmse, r2
#row 1: sl
model.sel.tbl.region[1,7]=round(null.RMSE.sl.region[1,2],3)
model.sel.tbl.region[1,8]=round(null.R2.sl.region[1,2],3)
#row 2: sigma sl
model.sel.tbl.region[2,7]=round(null.RMSE.sigma.sl.region[1,2],3)
model.sel.tbl.region[2,8]=round(null.R2.sigma.sl.region[1,2],3)
#row 3: disp
model.sel.tbl.region[3,7]=round(null.RMSE.disp.region[1,2],3)
model.sel.tbl.region[3,8]=round(null.R2.disp.region[1,2],3)
#row 4: sigma disp
model.sel.tbl.region[4,7]=round(null.RMSE.sigma.disp.region[1,2],3)
model.sel.tbl.region[4,8]=round(null.R2.sigma.disp.region[1,2],3)
#row 5: sigma disp
model.sel.tbl.region[5,7]=round(null.RMSE.tenavg.region[1,2],3)
model.sel.tbl.region[5,8]=round(null.R2.tenavg.region[1,2],3)

#write.csv(model.sel.tbl.region,paste(outdir,"ModelSelTbl_region_14DEC22.csv"))

model.sel.tbl.total=as.data.frame(matrix(nrow=5,ncol=16))
rownames(model.sel.tbl.total)=c("step length", "step length sigma", "displacement", "disp. sigma", "top average sl")
colnames(model.sel.tbl.total)=c("full RMSE random", "full r2 random", "drop 01 RMSE random", "drop 01 r2 random", "lasso RMSE random", "lasso r2 random", "null RMSE random", "null R2 random",
                          "full RMSE region", "full r2 region", "drop 01 RMSE region", "drop 01 r2 region", "lasso RMSE region", "lasso r2 region", "null RMSE region", "null R2region")

model.sel.tbl.total[1:5,1:8]=model.sel.tbl
model.sel.tbl.total[1:5,9:16]=model.sel.tbl.region

model.sel.tbl.total=model.sel.tbl.total[,-c(2,4,6,8,10,12,14,16)]

write.csv(model.sel.tbl.total,paste(outdir,"ModelSelTblTotal_26APR22.csv"))

#get opt params for each top model
sl.opt.params.kfold=read.csv(paste0(outdir,"05APR23_Runs/Random/sl_bestmodelparams.csv"))
sigma.sl.opt.params.kfold=read.csv(paste0(outdir,"05APR23_Runs/Random/sigma.sl_bestmodelparams.csv"))
disp.opt.params.kfold=read.csv(paste0(outdir,"05APR23_Runs/Random/disp_bestmodelparams.csv"))
sigma.disp.opt.params.kfold=read.csv(paste0(outdir,"05APR23_Runs/Random/sigma.disp_bestmodelparams.csv"))
tenavg.opt.params.kfold=read.csv(paste0(outdir,"05APR23_Runs/Random/tenavg_bestmodelparams.csv"))

sl.opt.params.kfold_region=read.csv(paste0(outdir,"05APR23_Runs/Region/sl_bestmodelparams.csv"))
sigma.sl.opt.params.kfold_region=read.csv(paste0(outdir,"05APR23_Runs/Region/sigma.sl_bestmodelparams.csv"))
disp.opt.params.kfold_region=read.csv(paste0(outdir,"05APR23_Runs/Region/disp_bestmodelparams.csv"))
sigma.disp.opt.params.kfold_region=read.csv(paste0(outdir,"05APR23_Runs/Region/sigma.disp_bestmodelparams.csv"))
tenavg.opt.params.kfold_region=read.csv(paste0(outdir,"05APR23_Runs/Region/tenavg_bestmodelparams.csv"))

#read in X vecs as needed
#RANDOM: all full model X_vec
#sl-full region
#sl.sigma-full random
#displ-random drop 01
#tenavg-full region

#Random
#sl-full
#sl sigma-full
#disp-drop 01
#disp sigma- full
#top avg-full

#Region
#sl-full
#sl sigma-full
#disp-lasso
#disp sigma- full
#top avg-full

#sl,sl.sigma,tenavg
X_vec.start=c(3,4,6:13,18,21:35,38:44,46:48)

#displ random
X.vec.disp.01drop=read.csv(paste0(outdir,"05APR23_Runs/Random/X_vec_01Drop/Xdisp.csv"))
X.vec.disp.01drop=X.vec.disp.01drop[,2]

#displ region
X.vec.disp.lasso.region=read.csv(paste0(outdir,"05APR23_Runs/Region/X_vec_postlasso/Xdisp_lasso.csv"))
X.vec.disp.lasso.region=X.vec.disp.lasso.region[,2]

#REGION: 
#sl region 01 drop
#X.vec.sl_01drop.region=read.csv(paste0(outdir,"05APR23_Runs/Region/X_vec_01Drop/Xsl.csv"))
#X.vec.sl_01drop.region=X.vec.sl_01drop.region[,2]
#sigma sl lasso region - full model, X_vec.start
#X.vec.sigma.sl.lasso.region=read.csv(paste0(outdir,"X_vec_postlasso_13DEC22_region/Xsigmasl_lasso.csv"))
#X.vec.sigma.sl.lasso.region=X.vec.sigma.sl.lasso.region[,2]
#disp region 01 drop region
#X.vec.disp.01drop.region=read.csv(paste0(outdir,"05APR23_Runs/Region/X_vec_01Drop/Xdisp.csv"))
#X.vec.disp.01drop.region=X.vec.disp.01drop.region[,2]
#tenavg lasso region
#X.vec.tenavg.01drop.region=read.csv(paste0(outdir,"05APR23_Runs/Region/X_vec_01Drop/Xtenavg.csv"))
#X.vec.tenavg.01drop.region=X.vec.tenavg.01drop.region[,2]


gbm.sl=gbm.fixed(data=pigsums2, gbm.x=X_vec.start, gbm.y=which(colnames(pigsums)=="sl_"),
                 learning.rate=sl.opt.params.kfold$learning.rate, 
                 tree.complexity=sl.opt.params.kfold$tree.complexity, 
                 n.trees=sl.opt.params.kfold$n.trees,
                 bag.fraction=sl.opt.params.kfold$bag.fraction,
                 family="poisson") 
gbm.disp=gbm.fixed(data=pigsums2, gbm.x=X.vec.disp.01drop, gbm.y=which(colnames(pigsums)=="displacement"),
                   learning.rate=disp.opt.params.kfold$learning.rate, 
                   tree.complexity=disp.opt.params.kfold$tree.complexity, 
                   n.trees=disp.opt.params.kfold$n.trees,
                   bag.fraction=disp.opt.params.kfold$bag.fraction,
                   family="poisson") 
gbm.sigma.sl=gbm.fixed(data=cutoff_sl, gbm.x=X_vec.start, gbm.y=which(colnames(pigsums)=="sigma_sl"),
                       learning.rate=sigma.sl.opt.params.kfold$learning.rate, 
                       tree.complexity=sigma.sl.opt.params.kfold$tree.complexity, 
                       n.trees=sigma.sl.opt.params.kfold$n.trees,
                       bag.fraction=sigma.sl.opt.params.kfold$bag.fraction,
                       family="gaussian") 
gbm.sigma.disp=gbm.fixed(data=cutoff_sl, gbm.x=X_vec.start, gbm.y=which(colnames(pigsums)=="sigma_disp"),
                         learning.rate=sigma.disp.opt.params.kfold$learning.rate, 
                         tree.complexity=sigma.disp.opt.params.kfold$tree.complexity, 
                         n.trees=sigma.disp.opt.params.kfold$n.trees,
                         bag.fraction=sigma.disp.opt.params.kfold$bag.fraction,
                         family="gaussian") 
gbm.tenavg=gbm.fixed(data=pigsums, gbm.x=X_vec.start, gbm.y=which(colnames(pigsums)=="tenavg"),
                     learning.rate=tenavg.opt.params.kfold$learning.rate, 
                     tree.complexity=tenavg.opt.params.kfold$tree.complexity, 
                     n.trees=tenavg.opt.params.kfold$n.trees,
                     bag.fraction=tenavg.opt.params.kfold$bag.fraction,
                     family="poisson") 
#Region
#sl-full
#sl sigma-full
#disp-lasso
#disp sigma- full
#top avg-full
#region kfold models
gbm.sl_region=gbm.fixed(data=pigsums2, gbm.x=X_vec.start, gbm.y=which(colnames(pigsums)=="sl_"),
                 learning.rate=sl.opt.params.kfold_region$learning.rate, 
                 tree.complexity=sl.opt.params.kfold_region$tree.complexity, 
                 n.trees=sl.opt.params.kfold_region$n.trees,
                 bag.fraction=sl.opt.params.kfold_region$bag.fraction,
                 family="poisson")

gbm.sigma.sl_region=gbm.fixed(data=cutoff_sl, gbm.x=X_vec.start, gbm.y=which(colnames(pigsums)=="sigma_sl"),
                        learning.rate=sigma.sl.opt.params.kfold_region$learning.rate, 
                        tree.complexity=sigma.sl.opt.params.kfold_region$tree.complexity, 
                        n.trees=sigma.sl.opt.params.kfold_region$n.trees,
                        bag.fraction=sigma.sl.opt.params.kfold_region$bag.fraction,
                        family="gaussian") 

gbm.disp_region=gbm.fixed(data=pigsums2, gbm.x=X.vec.disp.lasso.region, gbm.y=which(colnames(pigsums)=="displacement"),
                        learning.rate=disp.opt.params.kfold_region$learning.rate, 
                        tree.complexity=disp.opt.params.kfold_region$tree.complexity, 
                        n.trees=disp.opt.params.kfold_region$n.trees,
                        bag.fraction=disp.opt.params.kfold_region$bag.fraction,
                        family="poisson") 

gbm.sigma.disp_region=gbm.fixed(data=cutoff_sl, gbm.x=X_vec.start, gbm.y=which(colnames(pigsums)=="sigma_disp"),
                          learning.rate=sigma.disp.opt.params.kfold_region$learning.rate, 
                          tree.complexity=sigma.disp.opt.params.kfold_region$tree.complexity, 
                          n.trees=sigma.disp.opt.params.kfold_region$n.trees,
                          bag.fraction=sigma.disp.opt.params.kfold_region$bag.fraction,
                          family="gaussian") 

gbm.tenavg_region=gbm.fixed(data=pigsums2, gbm.x=X_vec.start, gbm.y=which(colnames(pigsums)=="tenavg"),
                     learning.rate=tenavg.opt.params.kfold$learning.rate, 
                     tree.complexity=tenavg.opt.params.kfold$tree.complexity, 
                     n.trees=tenavg.opt.params.kfold$n.trees,
                     bag.fraction=tenavg.opt.params.kfold$bag.fraction,
                     family="poisson") 

#export importance tables/figures
rel.inf.sl=summary(gbm.sl)
rel.inf.sigma.sl=summary(gbm.sigma.sl)
rel.inf.disp=summary(gbm.disp)
rel.inf.sigma.disp=summary(gbm.sigma.disp)
rel.inf.tenavg=summary(gbm.tenavg)

rel.inf.sl_region=summary(gbm.sl_region)
rel.inf.sigma.sl_region=summary(gbm.sigma.sl_region)
rel.inf.disp_region=summary(gbm.disp_region)
rel.inf.sigma.disp_region=summary(gbm.sigma.disp_region)
rel.inf.tenavg_region=summary(gbm.tenavg_region)

#for top models, transform RI into percentages
rel.inf.sl_region$rel.prop=rel.inf.sl_region$rel.inf/rel.inf.sl_region$rel.inf[1]
rel.inf.sigma.sl$rel.prop=rel.inf.sigma.sl$rel.inf/rel.inf.sigma.sl$rel.inf[1]
rel.inf.disp$rel.prop=rel.inf.disp$rel.inf/rel.inf.disp$rel.inf[1]
rel.inf.sigma.disp_region$rel.prop=rel.inf.sigma.disp_region$rel.inf/rel.inf.sigma.disp_region$rel.inf[1]

rel.inf.sl_region=rel.inf.sl_region[order(rel.inf.sl_region[,1]),]
rel.inf.sigma.sl=rel.inf.sigma.sl[order(rel.inf.sigma.sl[,1]),]
rel.inf.disp=rel.inf.disp[order(rel.inf.disp[,1]),]
rel.inf.sigma.disp_region=rel.inf.sigma.disp_region[order(rel.inf.sigma.disp_region[,1]),]

colnames(rel.inf.sl_region)[3]="rel.prop.sl"
colnames(rel.inf.sigma.sl)[3]="rel.prop.sigmasl"
colnames(rel.inf.disp)[3]="rel.prop.disp"
colnames(rel.inf.sigma.disp_region)[3]="rel.prop.sigmadisp"

join1=full_join(rel.inf.sl_region,rel.inf.sigma.sl,by="var")
join2=full_join(join1,rel.inf.disp,by="var")
join3=full_join(join2,rel.inf.sigma.disp_region,by="var")

join4=join3[,c(1,3,5,7,9)]

join5=as.data.frame(pivot_longer(join4,2:5))

#Fix names, neaten
join5[join5$name=="rel.prop.disp",2]<-"Displ."
join5[join5$name=="rel.prop.sigmadisp",2]<-"Displ. var"
join5[join5$name=="rel.prop.sigmasl",2]<-"S.L. var"
join5[join5$name=="rel.prop.sl",2]<-"S.L."

join5[join5$var=="dayl..s.",1]<-"Daylight"
join5[join5$var=="daylrange",1]<-"Dayl. range"
join5[join5$var=="Droughtrange",1]<-"Drought range"
join5[join5$var=="LCsums",1]<-"LC sums"
join5[join5$var=="LCT_1",1]<-"LC 1"
join5[join5$var=="LCT_10",1]<-"LC 10"
join5[join5$var=="LCT_14",1]<-"LC 14"
join5[join5$var=="LCT_15",1]<-"LC 15"
join5[join5$var=="LCT_16",1]<-"LC 16"
join5[join5$var=="LCT_17",1]<-"LC 17"
join5[join5$var=="LCT_18",1]<-"LC 18"
join5[join5$var=="LCT_3",1]<-"LC 3"
join5[join5$var=="LCT_4",1]<-"LC 4"
join5[join5$var=="LCT_5",1]<-"LC 5"
join5[join5$var=="LCT_6",1]<-"LC 6"
join5[join5$var=="LCT_7",1]<-"LC 7"
join5[join5$var=="LCT_8",1]<-"LC 8"
join5[join5$var=="LCT_9",1]<-"LC 9"
join5[join5$var=="Mastrange",1]<-"Mast range"
join5[join5$var=="prcp..mm.day.",1]<-"Prcp"
join5[join5$var=="prcprange",1]<-"Prcp rng."
join5[join5$var=="Prd2",1]<-"Prox rd."
join5[join5$var=="Prd2range",1]<-"Prox rd. range"
join5[join5$var=="Prd2range",1]<-"Prox rd. range"
join5[join5$var=="Rug2",1]<-"Rugd."
join5[join5$var=="Rug2range",1]<-"Rugd. range"
join5[join5$var=="TempRange",1]<-"Temp. Range"
join5[join5$var=="tmax..deg.c.",1]<-"Tmax"
join5[join5$var=="tmaxrange",1]<-"Tmax range"
join5[join5$var=="tmin..deg.c.",1]<-"Tmin"
join5[join5$var=="tminrange",1]<-"Tmin range"

#Make heat map with proportional rel infl.
ggplot(join5, aes(name, var, fill= value)) + 
  geom_tile() +
  scale_fill_viridis(discrete=FALSE) +
  theme_minimal()+
  labs(x ="Response", y = "Predictor", fill="Rel. Inf. Prop.")+
  theme(axis.title.y = element_text("Arial", "bold", size=12),
        axis.title.x = element_text("Arial", "bold", size=12),
        legend.title = element_text("Arial", "bold", size=12),
        axis.text=element_text(size=12),
        legend.text=element_text(size=12))

#Neaten var names for plots
rel.inf.sl=neaten.var.names.full(rel.inf.sl)
rel.inf.sigma.sl=neaten.var.names.full(rel.inf.sigma.sl)
rel.inf.disp=neaten.var.names.full(rel.inf.disp)
rel.inf.sigma.disp=neaten.var.names.full(rel.inf.sigma.disp)
rel.inf.tenavg=neaten.var.names.full(rel.inf.tenavg)

rel.inf.sl_region=neaten.var.names.full(rel.inf.sl_region)
rel.inf.sigma.sl_region=neaten.var.names.full(rel.inf.sigma.sl_region)
rel.inf.disp_region=neaten.var.names.full(rel.inf.disp_region)
rel.inf.sigma.disp_region=neaten.var.names.full(rel.inf.sigma.disp_region)
rel.inf.tenavg_region=neaten.var.names.full(rel.inf.tenavg_region)

#remove variables with 0 influence from tables (full models)
RemoveZeroInfluence=function(relinftbl){
  relinftbl=relinftbl[(relinftbl$rel.inf!=0),]
  return(relinftbl)
}

rel.inf.sl=RemoveZeroInfluence(rel.inf.sl)
rel.inf.sl_region=RemoveZeroInfluence(rel.inf.sl_region)
rel.inf.sigma.sl=RemoveZeroInfluence(rel.inf.sigma.sl)
rel.inf.sigma.sl_region=RemoveZeroInfluence(rel.inf.sigma.sl_region)
rel.inf.disp=RemoveZeroInfluence(rel.inf.disp)
rel.inf.disp_region=RemoveZeroInfluence(rel.inf.disp_region)
rel.inf.sigma.disp=RemoveZeroInfluence(rel.inf.sigma.disp)
rel.inf.sigma.disp_region=RemoveZeroInfluence(rel.inf.sigma.disp_region)
rel.inf.tenavg=RemoveZeroInfluence(rel.inf.tenavg)
rel.inf.tenavg_region=RemoveZeroInfluence(rel.inf.tenavg_region)

#make formatted rel influence plots
ri.sl.random=ggplot(rel.inf.sl, aes(x = reorder(var,rel.inf), y = rel.inf, fill=rel.inf))+
  geom_bar(stat="identity",show.legend=FALSE)+
  coord_flip()+
  labs(x="Independent variables", y="Relative Influence", fill="Rel. Inf.")+
  scale_fill_viridis_c(begin=0.2,end=0.4)+
  theme_minimal()+
  theme(text = element_text(size = 28))

ri.sl.region=ggplot(rel.inf.sl_region, aes(x = reorder(var,rel.inf), y = rel.inf, fill=rel.inf))+
  geom_bar(stat="identity",show.legend=FALSE)+
  coord_flip()+
  labs(x="Independent variables", y="Relative Influence", fill="Rel. Inf.")+
  scale_fill_viridis_c(begin=0.2,end=0.4,option="magma")+
  theme_minimal()+
  theme(text = element_text(size = 28))


ri.sigmasl.random=ggplot(rel.inf.sigma.sl, aes(x = reorder(var,rel.inf), y = rel.inf, fill=rel.inf))+
  geom_bar(stat="identity",show.legend=FALSE)+
  coord_flip()+
  labs(x="Independent variables", y="Relative Influence", fill="Rel. Inf.")+
  scale_fill_viridis_c(begin=0.2,end=0.4)+
  theme_minimal()+
  theme(text = element_text(size = 28))

ri.sigmasl.region=ggplot(rel.inf.sigma.sl_region, aes(x = reorder(var,rel.inf), y = rel.inf, fill=rel.inf))+
  geom_bar(stat="identity",show.legend=FALSE)+
  coord_flip()+
  labs(x="Independent variables", y="Relative Influence", fill="Rel. Inf.")+
  scale_fill_viridis_c(begin=0.2,end=0.4,option="magma")+
  theme_minimal()+
  theme(text = element_text(size = 28))

ri.disp.random=ggplot(rel.inf.disp, aes(x = reorder(var,rel.inf), y = rel.inf, fill=rel.inf))+
  geom_bar(stat="identity",show.legend=FALSE)+
  coord_flip()+
  labs(x="Independent variables", y="Relative Influence", fill="Rel. Inf.")+
  scale_fill_viridis_c(begin=0.3,end=0.4)+
  theme_minimal()+
  theme(text = element_text(size = 28))

ri.disp.region=ggplot(rel.inf.disp_region, aes(x = reorder(var,rel.inf), y = rel.inf, fill=rel.inf))+
  geom_bar(stat="identity",show.legend=FALSE)+
  coord_flip()+
  labs(x="Independent variables", y="Relative Influence", fill="Rel. Inf.")+
  scale_fill_viridis_c(begin=0.2,end=0.4,option="magma")+
  theme_minimal()+
  theme(text = element_text(size = 28))

ri.sigmadisp.random=ggplot(rel.inf.sigma.disp, aes(x = reorder(var,rel.inf), y = rel.inf, fill=rel.inf))+
  geom_bar(stat="identity",show.legend=FALSE)+
  coord_flip()+
  labs(x="Independent variables", y="Relative Influence", fill="Rel. Inf.")+
  scale_fill_viridis_c(begin=0.2,end=0.4)+
  theme_minimal()+
  theme(text = element_text(size = 28))

ri.sigmadisp.region=ggplot(rel.inf.sigma.disp_region, aes(x = reorder(var,rel.inf), y = rel.inf, fill=rel.inf))+
  geom_bar(stat="identity",show.legend=FALSE)+
  coord_flip()+
  labs(x="Independent variables", y="Relative Influence", fill="Rel. Inf.")+
  scale_fill_viridis_c(begin=0.2,end=0.4,option="magma")+
  theme_minimal()+
  theme(text = element_text(size = 28))

ri.tenavg.random=ggplot(rel.inf.tenavg, aes(x = reorder(var,rel.inf), y = rel.inf, fill=rel.inf))+
  geom_bar(stat="identity")+
  coord_flip()+
  labs(x="Independent variables", y="Relative Influence")+
  theme_minimal()+
  theme(text = element_text(size = 28))

ri.tenavg.region=ggplot(rel.inf.tenavg_region, aes(x = reorder(var,rel.inf), y = rel.inf, fill=rel.inf))+
  geom_bar(stat="identity",show.legend=FALSE)+
  coord_flip()+
  labs(x="Independent variables", y="Relative Influence")+
  scale_fill_viridis_c(begin=0.3,end=0.4,option="magma")+
  theme_minimal()+
  theme(text = element_text(size = 28))

#make plotgrids
relinf.pg=plot_grid(
  ri.sl.random,
  ri.sl.region,
  ri.sigmasl.random,
  ri.sigmasl.region,
  ri.disp.random,
  ri.disp.region,
  ri.sigmadisp.random,
  ri.sigmadisp.region,
  ncol = 2,
  labels=c("a","b","c","d","e","f","g","h"),
  label_size=50
)

relinf.pg.sl=plot_grid(
  ri.sl.random,
  ri.sl.region,
  ri.sigmasl.random,
  ri.sigmasl.region,
  ncol = 2,
  labels=c("a","b","c","d"),
  label_size=50
)

relinf.pg.disp=plot_grid(
  ri.disp.random,
  ri.disp.region,
  ri.sigmadisp.random,
  ri.sigmadisp.region,
  ncol = 2,
  labels=c("a","b","c","d"),
  label_size=50
)

png(file=paste(home,"Outputs/YFigureOutputs/26APR23_FigureSet/rel_infl_charts/relinfplots_08MAY_sl.png",sep="/"),
    width=2550, height=2550)
relinf.pg.sl
dev.off()

png(file=paste(home,"Outputs/YFigureOutputs/26APR23_FigureSet/rel_infl_charts/relinfplots_08MAY_disp.png",sep="/"),
    width=2550, height=2550)
relinf.pg.disp
dev.off()







numplots=12
#Make partial dependence plots for each response top model, top 12
MakePdpGrid(numplots,rel.inf.sl_region,gbm.sl_region,pigsums2,sl.opt.params.kfold_region,c(-0.35,0.35))

MakePdpGrid(numplots,rel.inf.sigma.sl,gbm.sigma.sl,cutoff_sl,sigma.sl.opt.params.kfold,c(-0.35,0.35))

MakePdpGrid(numplots,rel.inf.disp,gbm.disp,pigsums2,disp.opt.params.kfold,c(-0.9,0.9))

MakePdpGrid(numplots,rel.inf.sigma.disp_region,gbm.sigma.disp_region,pigsums2,sigma.disp.opt.params.kfold,c(-0.05,0.05))

#MakePdpGrid(numplots,rel.inf.tenavg_region,gbm.tenavg_region,pigsums2,tenavg.opt.params.kfold_region,c(-0.5,0.5))

#######################################################
  
#gbm.interactions(gbm.disp_region)


names(gbm.sl$gbm.call)[1] <- "dataframe"
names(gbm.sigma.sl$gbm.call)[1] <- "dataframe"
names(gbm.disp$gbm.call)[1] <- "dataframe"

#Get RMSE and R2 values for out of site/sample predictions

#Random
#sl-full
#sl sigma-full
#disp-drop 01
#disp sigma- full
#top avg-full

#Region
#sl-full
#sl sigma-full
#disp-lasso
#disp sigma- full
#top avg-full

CVstats_sl.random=GetCVStats_Table(pigsums2,X_vec.start,"sl_",sl.opt.params.kfold,"poisson","random",studydf)
CVstats_sigma.sl.random=GetCVStats_Table(cutoff_sl,X_vec.start,"sigma_sl",sigma.sl.opt.params.kfold,"gaussian","random",studydf)
CVstats_disp.random=GetCVStats_Table(pigsums2,X.vec.disp.01drop,"displacement",disp.opt.params.kfold,"poisson","random",studydf)
CVstats_sigma.disp.random=GetCVStats_Table(cutoff_sl,X_vec.start,"sigma_disp",sigma.disp.opt.params.kfold,"gaussian","random",studydf)
#CVstats_tenavg.random=GetCVStats_Table(pigsums2,X_vec.start,"tenavg",tenavg.opt.params.kfold,"poisson","random",studydf)

CVstats_sl.region=GetCVStats_Table(pigsums2,X_vec.start,"sl_",sl.opt.params.kfold_region,"poisson","region",studydf)
CVstats_sigma.sl.region=GetCVStats_Table(cutoff_sl,X_vec.start,"sigma_sl",sigma.sl.opt.params.kfold_region,"gaussian","region",studydf)
CVstats_disp.region=GetCVStats_Table(pigsums2,X.vec.disp.lasso.region,"displacement",disp.opt.params.kfold_region,"poisson","region",studydf)
CVstats_sigma.disp.region=GetCVStats_Table(cutoff_sl,X_vec.start,"sigma_disp",sigma.disp.opt.params.kfold_region,"gaussian","region",studydf)
#CVstats_tenavg.region=GetCVStats_Table(pigsums2,X_vec.start,"tenavg",tenavg.opt.params.kfold_region,"poisson","region",studydf)

#combine CV method sets
CVstats_sl.total=rbind(CVstats_sl.random[[1]],CVstats_sl.region[[1]])
CVstats_sigma.sl.total=rbind(CVstats_sigma.sl.random[[1]],CVstats_sigma.sl.region[[1]])
CVstats_disp.total=rbind(CVstats_disp.random[[1]],CVstats_disp.region[[1]])
CVstats_sigma.disp.total=rbind(CVstats_sigma.disp.random[[1]],CVstats_sigma.disp.region[[1]])
#CVstats_tenavg.total=rbind(CVstats_tenavg.random[[1]],CVstats_tenavg.region[[1]])



#Remove Mitchell from these descriptive plots because only two points-- 
#difficult to get accurate estimate of prediction strength
#CVstats_sl.total=CVstats_sl.total[!(CVstats_sl.total$region=="Mitchell"),]
#CVstats_sigma.sl.total=CVstats_sigma.sl.total[!(CVstats_sigma.sl.total$region=="Mitchell"),]
#CVstats_disp.total=CVstats_disp.total[!(CVstats_disp.total$region=="Mitchell"),]
#CVstats_sigma.disp.total=CVstats_sigma.disp.total[!(CVstats_sigma.disp.total$region=="Mitchell"),]
#CVstats_tenavg.total=CVstats_tenavg.total[!(CVstats_tenavg.total$region=="Mitchell"),]

#Redo the means for both random/region
#CVstats_sl.total[CVstats_sl.total$type=="random",][25,3]=mean(CVstats_sl.total[CVstats_sl.total$type=="random",3][1:24])
#CVstats_sl.total[CVstats_sl.total$type=="region",][25,3]=mean(CVstats_sl.total[CVstats_sl.total$type=="region",3][1:24])

#CVstats_sigma.sl.total[CVstats_sigma.sl.total$type=="random",][23,3]=mean(CVstats_sigma.sl.total[CVstats_sigma.sl.total$type=="random",3][1:22])
#CVstats_sigma.sl.total[CVstats_sigma.sl.total$type=="region",][23,3]=mean(CVstats_sigma.sl.total[CVstats_sigma.sl.total$type=="region",3][1:22])

#CVstats_disp.total[CVstats_disp.total$type=="random",][25,3]=mean(CVstats_disp.total[CVstats_disp.total$type=="random",3][1:24])
#CVstats_disp.total[CVstats_disp.total$type=="region",][25,3]=mean(CVstats_disp.total[CVstats_disp.total$type=="region",3][1:24])

#CVstats_sigma.disp.total[CVstats_sigma.disp.total$type=="random",][23,3]=mean(CVstats_sigma.disp.total[CVstats_sigma.disp.total$type=="random",3][1:22])
#CVstats_sigma.disp.total[CVstats_sigma.disp.total$type=="region",][23,3]=mean(CVstats_sigma.disp.total[CVstats_sigma.disp.total$type=="region",3][1:22])
#CVstats_tenavg.total[CVstats_tenavg.total$type=="region",][25,3]=mean(CVstats_tenavg.total[CVstats_tenavg.total$type=="region",3][1:24])

#left join to get state added for ordering
regstats=unique(pigswsite[,c(3,7)])
CVstats_sl.total=left_join(CVstats_sl.total,regstats,by="region")
CVstats_sigma.sl.total=left_join(CVstats_sigma.sl.total,regstats,by="region")
CVstats_disp.total=left_join(CVstats_disp.total,regstats,by="region")
CVstats_sigma.disp.total=left_join(CVstats_sigma.disp.total,regstats,by="region")
#CVstats_tenavg.total=left_join(CVstats_tenavg.total,regstats,by="region")

CVstats_sl.total=unique(CVstats_sl.total[,-(5)])
CVstats_sigma.sl.total=unique(CVstats_sigma.sl.total[,-(5)])
CVstats_disp.total=unique(CVstats_disp.total[,-(5)])
CVstats_sigma.disp.total=unique(CVstats_sigma.disp.total[,-(5)])
#CVstats_tenavg.total=unique(CVstats_tenavg.total[,-(5)])

#order levels for plotting
CVstats_sl.total=CVstats_sl.total[order(CVstats_sl.total$State,CVstats_sl.total$region),]
CVstats_sl.total$region <- factor(CVstats_sl.total$region, levels=unique(CVstats_sl.total$region))
CVstats_sigma.sl.total=CVstats_sigma.sl.total[order(CVstats_sigma.sl.total$State,CVstats_sigma.sl.total$region),]
CVstats_sigma.sl.total$region <- factor(CVstats_sigma.sl.total$region, levels=unique(CVstats_sigma.sl.total$region))
CVstats_disp.total=CVstats_disp.total[order(CVstats_disp.total$State,CVstats_disp.total$region),]
CVstats_disp.total$region <- factor(CVstats_disp.total$region, levels=unique(CVstats_disp.total$region))
CVstats_sigma.disp.total=CVstats_sigma.disp.total[order(CVstats_sigma.disp.total$State,CVstats_sigma.disp.total$region),]
CVstats_sigma.disp.total$region <- factor(CVstats_sigma.disp.total$region, levels=unique(CVstats_sigma.disp.total$region))
#CVstats_tenavg.total=CVstats_tenavg.total[order(CVstats_tenavg.total$State,CVstats_tenavg.total$region),]
#CVstats_tenavg.total$region <- factor(CVstats_tenavg.total$region, levels=unique(CVstats_tenavg.total$region))

##Want to add some descriptive stats to change dot sizes in R2/RMSE plots
#need counts summed per study (ie total pig days)
#and total num pigs (count num unique pig id's per study)
#then left join to each CVstats table
#studycounts=as.data.frame(pigsums2 %>% 
#                            group_by(region) %>% 
#                            dplyr::summarise(numpigdays=sum(count), 
#                                             numpigs=n_distinct(id)))

#get avg.days.region


##remove means, won't plot nicely
CVstats_sl.total=CVstats_sl.total[CVstats_sl.total$region!="MEANS",]
CVstats_sigma.sl.total=CVstats_sigma.sl.total[CVstats_sigma.sl.total$region!="MEANS",]
CVstats_disp.total=CVstats_disp.total[CVstats_disp.total$region!="MEANS",]
CVstats_sigma.disp.total=CVstats_sigma.disp.total[CVstats_sigma.disp.total$region!="MEANS",]
#CVstats_tenavg.total=CVstats_tenavg.total[CVstats_tenavg.total$region!="MEANS",]

#replace table, do left join and then order by new study ID
region.names=c("LindsayNC",
               "LindsayCC",
               "Mitchell",
               "LindsaySC",
               "Tejon",
               "Noble",
               "Nate Snow",
               "Susan",
               "Camp Bullis",
               "Tyler",
               "Hartley",
               "Potts",
               "Noxubee",
               "Kurt",
               "Steve",
               "Bill",
               "FL2",
               "Raoul",
               "contact",
               "SC3",
               "Jim",
               "Jim2",
               "SREL_Vacuum",
               "SREL_contact",
               "PhD")
study.ID=as.character(c(1,
           2,
           3,
           4,
           5,
           6,
           7,
           8,
           9,
           10,
           11,
           12,
           13,
           14,
           15,
           16,
           17,
           18,
           19,
           20,
           21,
           22,
           23,
           24,
           25))
newID.tbl=data.frame(region=region.names,ID=study.ID)

#force 9,8,7 to go after 10 when character
newID.tbl[newID.tbl$ID==9,2]<-"09"
newID.tbl[newID.tbl$ID==8,2]<-"08"
newID.tbl[newID.tbl$ID==7,2]<-"07"
newID.tbl[newID.tbl$ID==6,2]<-"06"
newID.tbl[newID.tbl$ID==5,2]<-"05"
newID.tbl[newID.tbl$ID==4,2]<-"04"
newID.tbl[newID.tbl$ID==3,2]<-"03"
newID.tbl[newID.tbl$ID==2,2]<-"02"
newID.tbl[newID.tbl$ID==1,2]<-"01"

#do join with new ID.tbl to get new IDs replaced
CVstats_sl.total=left_join(CVstats_sl.total,newID.tbl,by="region")
CVstats_sl.total$region=CVstats_sl.total$ID
CVstats_sl.total=CVstats_sl.total[order(as.numeric(CVstats_sl.total$region)),]

CVstats_sigma.sl.total=left_join(CVstats_sigma.sl.total,newID.tbl,by="region")
CVstats_sigma.sl.total$region=CVstats_sigma.sl.total$ID
CVstats_sigma.sl.total=CVstats_sigma.sl.total[order(as.numeric(CVstats_sigma.sl.total$region)),]

CVstats_disp.total=left_join(CVstats_disp.total,newID.tbl,by="region")
CVstats_disp.total$region=CVstats_disp.total$ID
CVstats_disp.total=CVstats_disp.total[order(as.numeric(CVstats_disp.total$region)),]

CVstats_sigma.disp.total=left_join(CVstats_sigma.disp.total,newID.tbl,by="region")
CVstats_sigma.disp.total$region=CVstats_sigma.disp.total$ID
CVstats_sigma.disp.total=CVstats_sigma.disp.total[order(as.numeric(CVstats_sigma.disp.total$region)),]

studycounts=left_join(avg.days.region,newID.tbl,by="region")
studycounts$region=studycounts$ID
studycounts=studycounts[order(as.numeric(studycounts$region)),]
studycounts[studycounts$ID==1,1]<-"01"
studycounts[studycounts$ID==2,1]<-"02"
studycounts[studycounts$ID==3,1]<-"03"
studycounts[studycounts$ID==4,1]<-"04"
studycounts[studycounts$ID==5,1]<-"05"
studycounts[studycounts$ID==6,1]<-"06"
studycounts[studycounts$ID==7,1]<-"07"
studycounts[studycounts$ID==8,1]<-"08"
studycounts[studycounts$ID==9,1]<-"09"

#studycounts$region<-studycounts$ID
#3500, 9106, 13764, 17128, 23423, 27516

#test
#CVstats_sl.total$region<-as.character(CVstats_sl.total$region)
#studycounts$region<-as.character(studycounts$region)

#studycounts2<-studycounts
#studycounts=studycounts[,-c(2,3,6)]

colnames(studycounts)[which(colnames(studycounts)=="avgdays")]<-"numpigdays"
class(studycounts$region)
class(CVstats$region)

sl.rmse.dot=ggdraw(make.varsize.dotplots(CVstats_sl.total,studycounts,"RMSE","sl"))
sigmasl.rmse.dot=ggdraw(make.varsize.dotplots(CVstats_sigma.sl.total,studycounts, "RMSE","sigmasl"))
disp.rmse.dot=ggdraw(make.varsize.dotplots(CVstats_disp.total,studycounts, "RMSE","disp"))
sigmadisp.rmse.dot=ggdraw(make.varsize.dotplots(CVstats_sigma.disp.total,studycounts, "RMSE","sigmadisp"))
#tenavg.rmse.dot=ggdraw(make.varsize.dotplots(CVstats_tenavg.total,studycounts, "RMSE","tenavg"))

sl.r2.dot=ggdraw(make.varsize.dotplots(CVstats_sl.total,studycounts,"R2","sl"))
sigmasl.r2.dot=ggdraw(make.varsize.dotplots(CVstats_sigma.sl.total,studycounts, "R2","sigmasl"))
disp.r2.dot=ggdraw(make.varsize.dotplots(CVstats_disp.total,studycounts, "R2","disp"))
sigmadisp.r2.dot=ggdraw(make.varsize.dotplots(CVstats_sigma.disp.total,studycounts, "R2","sigmadisp"))
#tenavg.r2.dot=ggdraw(make.varsize.dotplots(CVstats_tenavg.total,studycounts, "R2","tenavg"))

ggarrange(sl.rmse.dot,
          sigmasl.rmse.dot,
          disp.rmse.dot,
          sigmadisp.rmse.dot,
          #tenavg.rmse.dot,
          common.legend=TRUE, 
          legend="right",
          ncol=2,nrow=2,labels=c("a","b","c","d"), 
          font.label=list(size=45,face="bold"))

ggarrange(sl.r2.dot,
          sigmasl.r2.dot,
          disp.r2.dot,
          sigmadisp.r2.dot,
          #tenavg.r2.dot,
          common.legend=TRUE, 
          legend="right",
          ncol=2,nrow=2,labels=c("a","b","c","d"), 
          font.label=list(size=45,face="bold"))

#####Pred vs obs plots!

#format prediction/test sets into dataframe for plotting
sl_predobs.df=CombinePredObs(CVstats_sl.random,CVstats_sl.region)
sigma.sl_predobs.df=CombinePredObs(CVstats_sigma.sl.random,CVstats_sigma.sl.region)
disp_predobs.df=CombinePredObs(CVstats_disp.random,CVstats_disp.region)
sigma.disp_predobs.df=CombinePredObs(CVstats_sigma.disp.random,CVstats_sigma.disp.region)
#tenavg_predobs.df=CombinePredObs(CVstats_tenavg.random,CVstats_tenavg.region)

#force single digits to go after 10 when character
newID.tbl[newID.tbl$ID==1,2]<-"01"
newID.tbl[newID.tbl$ID==2,2]<-"02"
newID.tbl[newID.tbl$ID==3,2]<-"03"
newID.tbl[newID.tbl$ID==4,2]<-"04"
newID.tbl[newID.tbl$ID==5,2]<-"05"
newID.tbl[newID.tbl$ID==6,2]<-"06"
newID.tbl[newID.tbl$ID==7,2]<-"07"
newID.tbl[newID.tbl$ID==8,2]<-"08"
newID.tbl[newID.tbl$ID==9,2]<-"09"


#do join with new ID.tbl to get new IDs replaced
sl_predobs.df=left_join(sl_predobs.df,newID.tbl,by="region")
sl_predobs.df$region=sl_predobs.df$ID
sl_predobs.df=sl_predobs.df[order(as.numeric(sl_predobs.df$region)),]

sigma.sl_predobs.df=left_join(sigma.sl_predobs.df,newID.tbl,by="region")
sigma.sl_predobs.df$region=sigma.sl_predobs.df$ID
sigma.sl_predobs.df=sigma.sl_predobs.df[order(as.numeric(sigma.sl_predobs.df$region)),]

disp_predobs.df=left_join(disp_predobs.df,newID.tbl,by="region")
disp_predobs.df$region=disp_predobs.df$ID
disp_predobs.df=disp_predobs.df[order(as.numeric(disp_predobs.df$region)),]

sigma.disp_predobs.df=left_join(sigma.disp_predobs.df,newID.tbl,by="region")
sigma.disp_predobs.df$region=sigma.disp_predobs.df$ID
sigma.disp_predobs.df=sigma.disp_predobs.df[order(as.numeric(sigma.disp_predobs.df$region)),]

#format for plotting histograms
sl_predobs.df2=as.data.frame(pivot_longer(sl_predobs.df,cols=c("test","preds")))
sigma.sl_predobs.df2=as.data.frame(pivot_longer(sigma.sl_predobs.df,cols=c("test","preds")))
disp_predobs.df2=as.data.frame(pivot_longer(disp_predobs.df,cols=c("test","preds")))
sigma.disp_predobs.df2=as.data.frame(pivot_longer(sigma.disp_predobs.df,cols=c("test","preds")))
#tenavg_predobs.df2=as.data.frame(pivot_longer(tenavg_predobs.df,cols=c("test","preds")))


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


#tenavg.po.scat=ggplot(tenavg_predobs.df, aes(test, preds, color = factor(region))) + 
#  geom_point(aes(color=factor(CVmethod)), size=0.5, alpha=1)+
#  geom_abline(intercept = 0, slope = 1)+
#  xlim(0,10000)+
#  ylim(0,10000)+
#  theme_minimal()+
#  facet_wrap(vars(region))+
#  geom_mark_ellipse(aes(color = factor(CVmethod)))+
#  scale_color_manual(values=c("#039e8c", "#832b9e"), name="CV method")+
#  labs(x="observed", y="predicted")+
#  theme(text = element_text(size = 12))+
#  theme(axis.text.x = element_text(angle = 90, size=10),
#        axis.text.y = element_text(size=9))

#tenavg.po.scat

png(file=paste(home,"Outputs/YFigureOutputs/26APR23_FigureSet/pred_obs_scatters/sl_po_scatter_26APR23.png",sep="/"),
    width=770, height=730)
sl.po.scat 
dev.off()

png(file=paste(home,"Outputs/YFigureOutputs/26APR23_FigureSet/pred_obs_scatters/sigmasl_po_scatter_26APR23.png",sep="/"),
    width=770, height=730)
sigmasl.po.scat
dev.off()

png(file=paste(home,"Outputs/YFigureOutputs/26APR23_FigureSet/pred_obs_scatters/disp_po_scatter_26APR23.png",sep="/"),
    width=770, height=730)
disp.po.scat
dev.off()

png(file=paste(home,"Outputs/YFigureOutputs/26APR23_FigureSet/pred_obs_scatters/sigmadisp_po_scatter_26APR23.png",sep="/"),
    width=770, height=730)
sigmadisp.po.scat
dev.off()

#########

sl_predobs.df2[sl_predobs.df2$name=="test",3]<-"observed"
sl_predobs.df2[sl_predobs.df2$name=="preds",3]<-"predicted"

sigma.sl_predobs.df2[sigma.sl_predobs.df2$name=="test",3]<-"observed"
sigma.sl_predobs.df2[sigma.sl_predobs.df2$name=="preds",3]<-"predicted"

disp_predobs.df2[disp_predobs.df2$name=="test",3]<-"observed"
disp_predobs.df2[disp_predobs.df2$name=="preds",3]<-"predicted"

sigma.disp_predobs.df2[sigma.disp_predobs.df2$name=="test",3]<-"observed"
sigma.disp_predobs.df2[sigma.disp_predobs.df2$name=="preds",3]<-"predicted"

#tenavg_predobs.df2[tenavg_predobs.df2$name=="test",3]<-"observed"
#tenavg_predobs.df2[tenavg_predobs.df2$name=="preds",3]<-"predicted"


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

#tenavg.dens<-tenavg_predobs.df2%>%
#  ggplot(aes(x=value, fill=name)) +
#  geom_density(alpha=0.5)+ 
#  labs(x= "Top 10% of daily step length means")+
#  theme_minimal()+
#  scale_fill_manual(values=c("#25858e", "#2ab07f"), name="CV method")+
#  facet_wrap(vars(CVmethod))+
#  theme(text = element_text(size = 40))

png(file=paste(home,"Outputs/YFigureOutputs/26APR23_FigureSet/pred_obs_densitys/pr_dens_grid_26APR23.png",sep="/"),
    width=2550, height=3300)
ggarrange(sl.dens,
          sigmasl.dens,
          disp.dens,
          sigmadisp.dens,
          #tenavg.dens,
          common.legend=TRUE, 
          legend="right",
          ncol=1,nrow=4,labels=c("a","b","c","d"), 
          font.label=list(size=45,face="bold"))
dev.off()

####Violin plots

#do join with new ID.tbl to get new IDs replaced
pigsums3=left_join(pigsums2,newID.tbl,by="region")
pigsums3$region=pigsums3$ID
pigsums3=pigsums3[order(as.numeric(pigsums3$region)),]

p1=ggplot(pigsums3, aes(x=fct_inorder(region), y=sl_)) + 
  geom_violin() + 
  coord_flip() + 
  geom_jitter(height=0,width=0.2,alpha=0.5) + theme(text = element_text(size = 40))+
  labs(x="study", y="average daily step length (m)")

p2=ggplot(pigsums3, aes(x=fct_inorder(region), y=sigma_sl)) + 
  geom_violin() + 
  coord_flip() + 
  geom_jitter(height=0,width=0.2,alpha=0.5) + theme(text = element_text(size = 40))+
  labs(x="study", y="average daily step length dispersion")

p3=ggplot(pigsums3, aes(x=fct_inorder(region), y=displacement)) + 
  geom_violin() + 
  coord_flip() + 
  geom_jitter(height=0,width=0.2,alpha=0.5) + theme(text = element_text(size = 40))+
  labs(x="study", y="average daily displacement (m)")

p4=ggplot(pigsums3, aes(x=fct_inorder(region), y=sigma_disp)) + 
  geom_violin() + 
  coord_flip() + 
  geom_jitter(height=0,width=0.2,alpha=0.5) + theme(text = element_text(size = 40))+
  labs(x="study", y="average daily displacement dispersion")

#p5=ggplot(pigsums2, aes(x=fct_inorder(num), y=tenavg)) + 
#  geom_violin() + 
#  coord_flip() + 
#  geom_jitter(height=0,width=0.2,alpha=0.5) + theme(text = element_text(size = 40))+
#  labs(x="study", y="average of top ten percent of daily average step lengths")

ggarrange(p1,p2,p3,p4,nrow=2,ncol=2,labels="auto",font.label=list(size=50,face="bold"),hjust=-0.2, vjust=1)




#scratch
######################################################3


density_grid.sl=MakePlotGrid(sl_predobs.df2,region.counts)
density_grid.disp=MakePlotGrid(disp_predobs.df2,region.counts)
density_grid.sl

#step length sigma
ggplot() + 
  geom_point(aes(y=df.sigma.sl$preds, x=df.sigma.sl$obs), color="black",size=0.5,alpha=1)+
  geom_abline(intercept = 0, slope = 1)+
  xlim(0,3)+
  ylim(0,3)

sigma.sl_predobs.df2%>%
  ggplot(aes(x=value, fill=name)) +
  geom_density(alpha=0.3)+ 
  labs(x= "Step Length")

ggplot() + 
  geom_point(aes(y=disp_predobs.df$preds, x=disp_predobs.df$test), color="black",size=0.5,alpha=1)+
  geom_abline(intercept = 0, slope = 1)+
  xlim(0,20000)+
  ylim(0,20000)+
  facet_wrap(vars(disp_predobs.df$region))

disp_predobs.df2%>%
  ggplot(aes(x=value, fill=name)) +
  geom_density(alpha=0.3)+ 
  scale_x_log10()+
  labs(x= "Step Length")

ggplot() + 
  geom_point(aes(y=df.sigma.disp$preds, x=df.sigma.disp$obs), color="black",size=0.5,alpha=1)+
  geom_abline(intercept = 0, slope = 1)+
  xlim(0,3)+
  ylim(0,3)

df2.sigma.disp%>%
  ggplot(aes(x=value, fill=name)) +
  geom_density(alpha=0.3)+ 
  scale_x_log10()+
  labs(x= "Step Length")

ggplot() + 
  geom_point(aes(y=df.tenavg$preds, x=df.tenavg$obs), color="black",size=0.5,alpha=1)+
  geom_abline(intercept = 0, slope = 1)+
  xlim(0,22000)+
  ylim(0,22000)

df2.tenavg%>%
  ggplot(aes(x=value, fill=name)) +
  geom_density(alpha=0.3)+ 
  scale_x_log10()+
  labs(x= "Step Length")

#coord with lexi on resource selection stuff




gbm.sl.int=gbm.interactions(gbm.sl)
gbm.sl.int$rank.list

gbm.sigmasl.int=gbm.interactions(gbm.sigma.sl)
gbm.sigmasl.int$rank.list

gbm.disp.int=gbm.interactions(gbm.disp)
gbm.disp.int$rank.list

pd <- partial(gbm.sl_region, pred.var = c("LCsums", "TCrange"),train=pigsums2,n.trees=disp.opt.params.kfold$n.trees)
pdp::partial

pd <- partial(gbm.tenavg_region, pred.var = c("LCsums", "TCrange"),train=pigsums2,n.trees=disp.opt.params.kfold$n.trees)
plotPartial(pd,train=pigsums2,rug=TRUE)





