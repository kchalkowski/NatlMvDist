
# Purpose --------------------------------

# This script does predictions from top models for 
#sl/disp responses for all watershed/seasons

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
library(ggplot2)
library(ggthemes)

## set working directories -------
#local:
home<-"/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/2_Projects/StatPigMvmt/Pipeline_R2"
filestr<-"04APR25_Runs"
outdir<-file.path(home,"4_Outputs",filestr)
objdir<-file.path(home,"2_Data","Objects")
pmat_path=file.path(outdir,"UncPredMats")
wash=readRDS(file.path(objdir,"wash_envcov_final.rds"))

#get helper shapefiles
usplot <- st_read(dsn = file.path(objdir,"usmapplot_best.shp"), layer = "usmapplot_best")
#format helper shapefiles for plotting
usplot<-st_as_sf(usplot)
usplot <- st_cast(usplot, "MULTIPOLYGON")
usplot=st_transform(usplot,crs=st_crs(wash))

## grab uncertainty matrices -------

CalcUncertainty<-function(pr){
prm=rowMeans(pr)
for(i in 1:nrow(pr)){
  x=pr[i,]
  varn=sqrt((sum((x-prm[i])^2)/length(pr)))/prm[i]
  if(i==1){
    varn_all=varn
  } else{
    varn_all=c(varn_all,varn)
  }
}
return(varn_all)
}

AvgSeasonUnc<-function(lf_subset){
  for(i in 1:length(lf_subset)){
    lfi=readRDS(lf_subset[i])
    varn=CalcUncertainty(lfi)
    if(i==1){
      varn_all=varn
    } else{
      varn_all=cbind(varn_all,varn)
    }
  } 
  varn_all<-rowMeans(varn_all)
  return(varn_all)
}

colname="sl_unc"
map=wash
string="sl_uncertainty"
WashMap<-function(map,colname,string){
    washmap=ggplot() + 
      geom_sf(data=usplot, fill="black")+
      geom_sf(data = map, aes(fill = eval(parse(text=colname))),lwd=0)+
      scale_fill_viridis_c(name=string)+
      geom_sf(data=usplot, fill=NA, color="#EBEBEB")+
      theme_map()+
      theme(plot.margin = unit(c(0.1,0.1,0.1,0.1), "in"))+
      theme(legend.position="bottom",
            legend.text = element_text(angle = 90,size=20),
            legend.title = element_text(size=20))
  return(washmap)
}


#pr=pr_sl
#map=wash
#string="sl_uncertainty"
UncPredMaps<-function(pr,map,string){
  varn_all=CalcUncertainty(pr)
  map$unc=varn_all
  washmap=WashMap(map,"unc",string)
  ggsave(file.path(outdir,"FigTab","UncMaps",paste0(string,".png")),plot=washmap,height=6.5,width=9,units="in")
}

lf=list.files(path=pmat_path,full.names=TRUE)
pr_sl=readRDS(lf[21])
pr_sigsl=readRDS(lf[16])
pr_sigdisp=readRDS(lf[11])
pr_disp=readRDS(lf[6])

UncPredMaps(pr_sl,wash,"S.L. Mean Uncertainty")
UncPredMaps(pr_sigsl,wash,"S.L. Var Uncertainty")
UncPredMaps(pr_sigdisp,wash,"Displ. Var Uncertainty")
UncPredMaps(pr_disp,wash,"Displ. Mean Uncertainty")

#Get values for manuscript
sl_u=CalcUncertainty(pr_sl)
sigsl_u=CalcUncertainty(pr_sigsl)
sigdisp_u=CalcUncertainty(pr_sigdisp)
disp_u=CalcUncertainty(pr_disp)

minmax<-function(unc){
  c(min(unc),max(unc))
}

minmax(sl_u)
minmax(sigsl_u)
minmax(sigdisp_u)
minmax(disp_u)



