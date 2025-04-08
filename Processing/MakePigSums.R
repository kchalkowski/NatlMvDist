##Averaging response variables step length and displacement for each pig in DisplacementResults.csv
##and, calculating sigma (dispersion) of step length and displacement distribution for each pig

###################
##Process Outline##
###################

#1-create dummy columns for categorical variables
#2-for each pig in a loop....
  #a-average numeric variables, and sum dummy columns for categorical variables
  #b-get proportions of total use of categorical variables for each pig
  #c-calculate sigma of displacement and step length distributions for each pig
#3-additional tidying/formatting
  #a-get step lengths and displacements as integers
  #b-remove any NAs
#4-write out the dataset for boosted regression analysis

#input: allsteps_daily.rds
#output: dailyPigSums.rds

################
##Script Setup##
################

#Load libraries
require(tidyverse)
require(dplyr)
require(job)
require(fastDummies)
require(fitdistrplus)

#set working directories
#local:
home<-"/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/StatPigMvmt"
outdir<-"/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/StatPigMvmt/Pipeline_R2/2_Data/Objects/"

#read data
dsteps=readRDS(paste0(outdir,"allsteps_daily.rds"))

#get vector of unique pig ids
dpigidvec=unique(dsteps$animalid)

#make id col for pig id, period, season
dsteps$animalid_per_season=paste(dsteps$animalid,dsteps$period,dsteps$season,sep="_")

##################
##Make functions##
##################
#GABi_FI9_p1_q3
#stepvec=steps[steps$animalid_per_season=="21958_d_Y048_p1_q4",]$sl_
#stepvec=steps[steps$animalid_per_season=="01_10_RaoulFL2_p1_q2",]$sl_
getMean<-function(stepvec){
  stepvec<-stepvec[!is.na(stepvec)]
  if(length(stepvec)>6){
    fit.gamma <- fitdist(stepvec+0.0001, distr = "gamma", method = "mle",lower=c(0,0), start=list(shape=1,rate=1))
    sigma=1/fit.gamma$estimate[[1]]
    meang=fit.gamma$estimate[[1]]/fit.gamma$estimate[[2]]
  } else{
    sigma=NA
    meang=NA
  }
  
  return(meang)
}

getDispersion<-function(stepvec){
  stepvec<-stepvec[!is.na(stepvec)]
  if(length(stepvec)>6){
    fit.gamma <- fitdist(stepvec+0.0001, distr = "gamma", method = "mle",lower=c(0,0), start=list(shape=1,rate=1))
    sigma=1/fit.gamma$estimate[[1]]
    meang=fit.gamma$estimate[[1]]/fit.gamma$estimate[[2]]
  } else{
    sigma=NA
    meang=NA
  }
  
  return(sigma)
}

#steps=dsteps
#steps=dsteps[dsteps$animalid_per_season=="GABi_FI9_p1_q3",]
steps=dsteps
#sl_mean=getMean(steps$sl_)
#MakePigSums(steps,idcol,keep_cols,env_cols,step_cols)
MakePigSums<-function(steps,idcol,keep_cols,env_cols,step_cols){
steps2=steps %>% 
    dplyr::group_by(animalid_per_season) %>%
    dplyr::summarise(
      across(keep_cols,first),
      across(env_cols,mean,.names = "mean_{.col}"),
      across(env_cols,var,.names = "var_{.col}"),
      not_na_sl = sum(!is.na(sl_)),
      not_na_displ = sum(!is.na(displ)))
      #sl_mean=getMean(sl_),
      #sl_disp=getDispersion(sl_),
      #displ_mean=getMean(displ),
      #displ_disp=getDispersion(displ))

ids=unique(steps$animalid_per_season)
df_r=data.frame("animalid_per_season"=ids,"sl_mean"=NA,"sl_disp"=NA,"displ_mean"=NA,"displ_disp"=NA)
for(i in 1:length(ids)){
sid=steps[steps$animalid_per_season==ids[i],]
#df_r$animalid_per_season[i]=ids[i]
df_r$sl_mean[i]=getMean(sid$sl_)
df_r$sl_disp[i]=getDispersion(sid$sl_)
df_r$displ_mean[i]=getMean(sid$displ)
df_r$displ_disp[i]=getDispersion(sid$displ)
}

steps2=left_join(steps2,df_r,by="animalid_per_season")

return(steps2)
}



###########################
##Get Summaries by Animal##
###########################

#daily summaries per pig
id_col="animalid_per_season"
keep_cols=c(
  which(colnames(dsteps)=="animalid"),
  which(colnames(dsteps)=="study"),
  which(colnames(dsteps)=="sex"),
  which(colnames(dsteps)=="age"),
  which(colnames(dsteps)=="state"),
  which(colnames(dsteps)=="season"),
  which(colnames(dsteps)=="period")
)
env_cols=c(18:30,42:56) #tc:tmin,lc98:lc
step_cols=c(which(colnames(dsteps)=="sl_"),
            which(colnames(dsteps)=="displ"))

#recode drought as numeric, make 'None' -1
dsteps$drt[dsteps$drt=="None"]<-(-1)
dsteps$drt[dsteps$drt=="D0"]<-0
dsteps$drt[dsteps$drt=="D1"]<-1
dsteps$drt[dsteps$drt=="D2"]<-2
dsteps$drt[dsteps$drt=="D3"]<-3
dsteps$drt<-as.integer(dsteps$drt)

dstep.sums=MakePigSums(dsteps,idcol,keep_cols,env_cols,step_cols)

#Note:
#NAs are generated for sl or displ when there are fewer than 6 steps for a given animal/period/season

#lots of NAs in var cols correspond to zero variance
dstep.sums[,grep("var_",colnames(dstep.sums))][is.na(dstep.sums[,grep("var_",colnames(dstep.sums))])]<-0

#cols that shouldn't have NAs, verify that there are no NAs
#all columns except age and sex (some indiv didn't have this info)
test=dstep.sums[,c(1:3,6:64)]
test[!complete.cases(test),]

#save out
saveRDS(as.data.frame(dstep.sums),paste0(outdir,"dailyPigSums.rds"))




