
# Purpose ---------------------

#The purpose of this script is to make a summary table of all the studies/pigs in dataset


# Setup ---------------------

#libraries
library(dplyr)
library(tidyr)

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
studyref=read.csv(file.path(home,"2_Data","Input","Pigs","studydf.csv"))

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

studydf=pigswsite %>%
  group_by(state,study,sex) %>%
  dplyr::summarise(
    Npigs=n_distinct(animalid)
  )
sdf2=pivot_wider(studydf,id_cols=1:2,names_from=c(3),values_from=c(4))

studydf2=pigswsite %>%
  group_by(state,study) %>%
  dplyr::summarise(
    Ngeo=n(),
    strtdt=as.Date(min(datetime)),
    enddt=as.Date(max(datetime))
  )

sdf3=left_join(studydf2,sdf2,by=c("state","study"))
colnames(sdf3)[c(6,7,8)]=c("NF","NM","NU")
sdf3=as.data.frame(sdf3)
sdf3[is.na(sdf3)]<-0

sdf3$num=c(1,2,3,0,0,0,15,16,11,14,13,24,12,5,0,0,0,0,7,9,6,4,8,10)
sdf3[sdf3$study=="Contact_FL",]$num<-22
sdf3[sdf3$study=="Raoul_FL1",]$num<-20
sdf3[sdf3$study=="Raoul_FL2",]$num<-21
sdf3[sdf3$study=="Jim_SC",]$num<-18
sdf3[sdf3$study=="Kilgo_USFS_SC3",]$num<-17
sdf3[sdf3$study=="SREL_Contact_SC",]$num<-20
sdf3[sdf3$study=="SREL_Vacuum_SC",]$num<-19

#Save key for GBM figure outputs
saveRDS(sdf3[,c(2,9)],file.path(objdir,"studynumkey.rds"))

#save outputs
saveRDS(sdf3,file.path(objdir,"studydf.rds"))
write.csv(sdf3,file.path(outdir,filestr,"FigTab","study_table.csv"))
