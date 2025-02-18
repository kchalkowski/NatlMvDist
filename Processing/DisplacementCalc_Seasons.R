##Calculating displacements for U.S. Pig Movement Modeling and Prediction

#Input: geolocs_ntl.csv
#Output: geolocsnatl_wDispl.rds

###################
##Process Outline##
###################

#1. set key identifier for each season
#2. find first day of tracking within each season, get average location over first 7 days
#3. by season, calc displacement from starting location for each pig

################
##Script Setup##
################

require(sf)
require(lubridate)
require(dplyr)

home<-"/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/StatPigMvmt/"
indir<-"/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/StatPigMvmt/Pipeline_R2/2_Data/Input/Pigs/"
outdir<-"/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/StatPigMvmt/Pipeline_R2/2_Data/Objects/"

pigs=read.csv(paste0(indir,"geolocsnatl.csv"))
pigs=pigs[,-1]

pigs$datetime<-Neat.Dates.POSIXct(pigs$datetime,tz="UTC")
pigs$year=lubridate::year(pigs$datetime)

####################
## Set season key ##
####################

#q1-jan-mar
#q2-apr-jun
#q3-jul-sep
#q4-oct-dec

months=lubridate::month(pigs$datetime)

pigs$season=NA
pigs$season[which(months==1|months==2|months==3)]<-"q1"
pigs$season[which(months==4|months==5|months==6)]<-"q2"
pigs$season[which(months==7|months==8|months==9)]<-"q3"
pigs$season[which(months==10|months==11|months==12)]<-"q4"

##################################
## Summarise pigs in each season ##
##################################
#summarise pigs in each season
#if there is <X duration of tracking in any one season, remove points from dataset

#cutoff: avg first 7 days each season
#dont calc displacement for pigseasons with < 14 days total tracking

#Group by animal,year,season
#order by date
#get start date within each pig/season
#calculate diff time between that point and all other points
pigs2=pigs %>% 
  dplyr::group_by(study,animalid, year, season) %>%
  dplyr::arrange(datetime,.by_group=TRUE) %>%
  mutate(strtdt=min(datetime)) %>%
  mutate(difftime=as.integer(difftime(datetime,strtdt,units="days")))

#set column with 1/0 if is in start week
#these will be averaged to get mean start point location
#NA in difftime means is first point, set to 0
pigs2$difftime[is.na(pigs2$difftime)]<-0
pigs2$is.startwk <- 0
pigs2[pigs2$difftime<=7,]$is.startwk <- 1

##############################################
## Avg first 7 days locs within each season ##
##############################################

#convert lat long to albers equal area conic (epsg 5070)
pigsf<-st_as_sf(pigs2,coords=c("longitude","latitude"),crs=st_crs(4326))
pigsf<-st_transform(pigsf,crs=st_crs(5070))
pigs2$X=st_coordinates(pigsf)[,1]
pigs2$Y=st_coordinates(pigsf)[,2]

#get start location for each pig/season
  start.locs=
    pigs2[pigs2$is.startwk==1,] %>% 
  group_by(study,animalid, season, year) %>%
  dplyr::arrange(datetime,.by_group=TRUE) %>%
  dplyr::summarize(startX=mean(X),startY=mean(Y))

#rejoin back to pigs with key
  start.locs$key=paste(start.locs$animalid,start.locs$season,start.locs$year,sep="_")
  pigs3=left_join(pigs2,start.locs)

##############################################
#### Calculate displacement by pig/season ####
##############################################

#group by pig season
  #get displ, use pythagorean method
pigs4=pigs3 %>% group_by(animalid, year, season) %>%
    dplyr::arrange(datetime,.by_group=TRUE) %>%
    mutate(displ=sqrt((X-startX)^2+(Y-startY)^2)) %>%
    as.data.frame()
  
################################################
#### Set displ as NA for certain conditions ####
################################################

  #displ should be NA for: all starting weeks (these are average start points)
  #any pig/seasons with ndays < 14
  
alsum=
    pigs4 %>% 
  group_by(study,animalid) %>%
  dplyr::summarise(strtdt=min(datetime),enddt=max(datetime))
alsum$dur=difftime(alsum$enddt,alsum$strtdt,units="days")

sesum=
  pigs4 %>% 
  group_by(study,animalid,season,year) %>%
  dplyr::summarise(strtdt=min(datetime),enddt=max(datetime))

sesum$dur=difftime(sesum$enddt,sesum$strtdt,units="days")
sesum$key=paste(sesum$animalid,sesum$season,sesum$year,sep="_")
dodispkey=sesum[sesum$dur>=14,]$key

#det which ones to do displacement on
pigs4$do.disp<-0
pigs4$key=paste(pigs4$animalid,pigs4$season,pigs4$year,sep="_")
pigs4[pigs4$key%in%dodispkey,]$do.disp<-1
pigs4[pigs4$difftime<7,]$do.disp<-0

#set some disp to na that meet conditions
pigs4[pigs4$do.disp==0,]$displ<-NA

#############################
#### Clean up formatting ####
#############################

#remove unneeded key columns
#pigs4=pigs4[,c(1:21,23,28,29:33)]
pigs4=
  pigs4[,-c(which(colnames(pigs4)=="key"),
  which(colnames(pigs4)=="do.disp"),
  which(colnames(pigs4)=="is.startwk"),
  which(colnames(pigs4)=="difftime"),
  which(colnames(pigs4)=="strtdt")
)]

#Write out as rds
outdir<-"/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/StatPigMvmt/Pipeline_R2/2_Data/Objects/"
saveRDS(pigs4,paste0(outdir,"geolocsnatl_wDispl.rds"))


