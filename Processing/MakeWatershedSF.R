#Set path of pipeline
home<-"/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/StatPigMvmt/Pipeline_R2"

# Purpose -------------------------------

#the purpose of this script is to make an sf object of US watersheds
#within the wild pig distribution, joined with env covs, 
#to predict pig sl/displacement mean and dispersion for 
#national pig movement prediction

# Input/Output -----------------------------------------------------------------

#inputs: multiple env cov files, watershed file, pig distribution file
#output: sf watershed file with env covs, wash_envcov_final.rds

# Process ----------------------------------------------------------------------

#1. gather inputs
    # a. watersheds shapefile
    # b. u.s. wp county distribution
    # c. all env cov files

#2. intersect watersheds with wp county distribution

#3. intersect/overlay cropped wash with static env covs
    
#4. intersect/overlay cropped wash with temporally varying env covs

# Setup ------------------------------------------------------------------------

#set file paths 
indir=file.path(home,"2_Data","Input",fsep = .Platform$file.se)
objdir=file.path(home,"2_Data","Objects",fsep = .Platform$file.se)

#load libraries
if (!require("pacman")) install.packages("pacman")
pacman::p_load(raster,sf,mactheknife,exactextractr)
library(daymetr)

#read env cov filenames
envls=list.files(file.path(indir,"Env_Cov",fsep = .Platform$file.se),full.names=TRUE)

#pig distribution shapefile, merged
pigd<-st_read(mactheknife::resolve_alias(envls[grep("merged_conus_alias",envls)]))

#HUC 12 wash data
#wash<-st_read(mactheknife::resolve_alias(envls[grep("WBD_HUC12_GCPO10k_alias",envls)]))
#wash<-st_read("~/Downloads/WBD_National_GPKG/WBD_National_GPKG.gpkg")
wash=st_read("~/Downloads/WBD_National_GPKG/WBD_National_GPKG.gpkg",layer="WBDHU12")

#Env covs for joining spat covs
#NLCD 
lcd<-raster(mactheknife::resolve_alias(envls[grep("NLCD_2021_alias",envls)]))
#TC
tc<-raster(mactheknife::resolve_alias(envls[grep("nlcd_tc_2021_alias",envls)]))
#Roads
rd1<-raster("/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Env_Data/Roads/2018/road_density_rasters/final_roads1.tif")
rd2<-raster("/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Env_Data/Roads/2018/road_density_rasters/final_roads2.tif")
rd3<-raster("/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Env_Data/Roads/2018/road_density_rasters/final_roads3.tif")

#Rgd
rgd<-raster(mactheknife::resolve_alias(envls[grep("rgd_5070_alias",envls)]))

#Drought files (temporally varying, multiple)
#drt<-st_read(mactheknife::resolve_alias(envls[grep("USDM_20211228_alias",envls)]))
read_alias<-function(x){
  st_read(mactheknife::resolve_alias(x))
}
drts=as.list(envls[grep("USDM_2021",envls)])
drt_list=lapply(drts,read_alias)

#transform/compareCRS for all drt 
#x=drt_list[[1]]
drt_transform=function(x){
  st_transform(x,crs=CRS("epsg:5070"))
}
drt_list=lapply(drt_list,drt_transform)

#Get formatted filenames into list
env_names=list.files(file.path(indir,"Env_Cov",fsep = .Platform$file.se),full.names=FALSE)
drt_names=env_names[grep("USDM_2021",env_names)]
drt_names=stringr::str_sub(drt_names,10L,-5L)
names(drt_list)=drt_names

#get indices for each quarter
#q1-> months 01-03
q1i=grep("^01|^02|^03",drt_names)
#q2-> months 04-06
q2i=grep("^04|^05|^06",drt_names)
#q3-> months 07-09
q3i=grep("^07|^08|^09",drt_names)
#q4-> months 10-12
q4i=grep("^10|^11|^12",drt_names)

#split drt_list into quadrants
drt_listq1=drt_list[q1i]
drt_listq2=drt_list[q2i]
drt_listq3=drt_list[q3i]
drt_listq4=drt_list[q4i]

#Mast
mast<-raster(mactheknife::resolve_alias(envls[grep("mast_5070_alias",envls)]))

#LCD and tc need reset crs to aea5070
crs(lcd)=CRS("epsg:5070")
crs(tc)=CRS("epsg:5070")
crs(rgd)=CRS("epsg:5070")
rds=st_transform(rds,crs=st_crs(lcd))
#drt=st_transform(drt,crs=st_crs(lcd))
wash=st_transform(wash,crs=st_crs(lcd))
#compare CRS for all spatials
compareCRS(lcd,wash)
#compareCRS(lcd,rds)
compareCRS(lcd,rgd)
compareCRS(lcd,mast)
compareCRS(pigd,wash)
pigd=st_transform(pigd,crs=st_crs(wash))
drt_compare=function(x){
  compareCRS(x,lcd)
}
lapply(drt_list,drt_compare)
# Create watershed basemap -------------------------------

#crop by pig distribution polygon
washsf=st_intersection(wash,pigd)

#check geometries
for(shape in 1:nrow(washsf)){
  sf1 <- washsf[shape,,drop=FALSE]
  if(!st_is_valid(sf1)){
    print(paste0(shape," is not valid"))
  }
}

#Remove any watersheds that are islands, water, or closed basin
washsf=washsf[washsf$hutype_description!="Water"&
         washsf$hutype_description!="Island"&
         washsf$hutype_description!="Closed Basin",]

#Get centroids of each watershed-- use this for joining temporal data
for(shape in 1:nrow(washsf)){
  sf1 <- washsf[shape,,drop=FALSE]
  sf1_ctr=st_centroid(sf1)
  if(shape==1){
    wctrs=sf1_ctr
  } else{
    wctrs=rbind(wctrs,sf1_ctr)
  }
}

#summarize daymetr output vals in list by quarter
dm_mv<-function(dat){
  
  dayl_mn=mean(dat[,"dayl..s."])
  dayl_var=var(dat[,"dayl..s."])
  
  prcp_mn=mean(dat[,"prcp..mm.day."])
  prcp_var=var(dat[,"prcp..mm.day."])
  
  tmin_mn=mean(dat[,"tmin..deg.c."])
  tmin_var=var(dat[,"tmin..deg.c."])
  
  tmax_mn=mean(dat[,"tmax..deg.c."])
  tmax_var=var(dat[,"tmax..deg.c."])
  
  out=data.frame("dayl_mn"=dayl_mn,
                 "dayl_var"=dayl_var,
                 "prcp_mn"=prcp_mn,
                 "prcp_var"=prcp_var,
                 "tmin_mn"=tmin_mn,
                 "tmin_var"=tmin_var,
                 "tmax_mn"=tmax_mn,
                 "tmax_var"=tmax_var)
  
  return(out)
}

# Savepoint 1 ------------------------------------
saveRDS(washsf,file.path(objdir,"washsf_tmp.rds"))
saveRDS(wctrs,file.path(objdir,"wctrs_tmp.rds"))
saveRDS(drt_list,file.path(objdir,"drt_list_tmp.rds"))

washsf=readRDS(file.path(objdir,"washsf_tmp.rds"))
wctrs=readRDS(file.path(objdir,"wctrs_tmp.rds"))
drt_list=readRDS(file.path(objdir,"drt_list_tmp.rds"))

# Set functions -------------------------------

#for aggregating lc values
#x=sf1
meanvar_cover<-function(x){
  x$coverage_fraction=round(x$coverage_fraction)
  if(any(x$coverage_fraction==1)){
  x=x[x$coverage_fraction==1,]
  } #if no full squares, just accumulate total, small polygon
  x=as.data.frame(table(x$value))
  x$total=sum(x$Freq)
  x$zero=x$total-x$Freq
  out=data.frame("lc"=x$Var1,"lc_mn"=rep(0,length(x$Var1)),"lc_var"=rep(0,length(x$Var1)))
  for(l in 1:length(x$Var1)){
    tv=c(rep(1,times=x$Freq[l]),rep(0,times=x$zero[l]))
    out[l,2]=mean(tv)
    out[l,3]=var(tv)
  }
 
  return(out)
}

#for all static covars except LC
meanvar_tc<-function(x){
  x$coverage_fraction=round(x$coverage_fraction)
  x=x[x$coverage_fraction==1,]
  mn_val=mean(x$value)
  var_val=var(x$value)
  
  return(data.frame("tc_mn"=mn_val,"tc_var"=var_val))
}


#for drt
intersect_pull_drt=function(x){
  drt=st_intersection(x,ctri)
  if(nrow(drt)>0){
    score=drt$DM
    if(length(score)>1){
      score=mean(score)
    }
  } else{
    score=0
    #drt_var=0
  }
  return(score)
}

# Overlay static env covs -------------------------------
#env covs that don't change by season
#Values needed
#all LC, mean and sd
#tc, mean and sd
#mast- mean and sd
#rgd- mean and sd
#rds- mean and sd

#loop through
#nrow(washsf) #total rows: 40052
#print(i)
for(i in 1:nrow(washsf)){
if(i%%100==0){print(i)}
#i=913
#subset
sf1 <- washsf[i,,drop=FALSE]

#get land class
#x=exact_extract(lcd, sf1, coverage_area = TRUE, summarize_df = TRUE, fun = sum_cover)
x=exact_extract(lcd, sf1, coverage_area = FALSE, summarize_df = TRUE, fun=meanvar_cover)
#x=exact_extract(lcd, sf1, coverage_area = FALSE)

x=tidyr::pivot_wider(x, names_from="lc",values_from=c(lc_mn,lc_var))

#rename to match washsf columns
#x[,1]=paste0("lc_mn_",x[,1])
#pivot wider to match format of washsf
#x=tidyr::pivot_wider(x, names_from="value",values_from="proportion")
#sf2=dplyr::bind_cols(sf1,x)

sf2=dplyr::bind_cols(sf1,x)

#tc, mean and sd
x=exact_extract(tc, sf2, coverage_area = FALSE, summarize_df = TRUE, fun=meanvar_tc)
#colnames(x)=c("tc","tc_mn","tc_var")
#x=tidyr::pivot_wider(x, names_from="tc",values_from=c(tc_mn,tc_var))

sf3=dplyr::bind_cols(sf2,x)

#rd1- mean and sd
rd1x=exact_extract(rd1, sf3, coverage_area = FALSE, summarize_df = TRUE, fun=meanvar_tc)
colnames(rd1x)=c("rd1_mn","rd1_var")

#rd2- mean and sd
rd2x=exact_extract(rd2, sf3, coverage_area = FALSE, summarize_df = TRUE, fun=meanvar_tc)
colnames(rd1x)=c("rd2_mn","rd2_var")

#rd3- mean and sd
rd3x=exact_extract(rd3, sf3, coverage_area = FALSE, summarize_df = TRUE, fun=meanvar_tc)
colnames(rd1x)=c("rd3_mn","rd3_var")

#mast- mean and sd
mastx=exact_extract(mast, sf3, coverage_area = FALSE, summarize_df = TRUE, fun=meanvar_tc)
colnames(mastx)=c("mast_mn","mast_var")

#rgd- mean and sd
rgdx=exact_extract(rgd, sf3, coverage_area = FALSE, summarize_df = TRUE, fun=meanvar_tc)
colnames(rgdx)=c("rgd_mn","rgd_var")

sf3=dplyr::bind_cols(sf3,rd1x)
sf3=dplyr::bind_cols(sf3,rd2x)
sf3=dplyr::bind_cols(sf3,rd3x)
sf3=dplyr::bind_cols(sf3,mastx)
sf3=dplyr::bind_cols(sf3,rgdx)
sf3=dplyr::bind_cols(sf3,rdsx)

if(i==1){
washsf_out=sf3
} else{
washsf_out=dplyr::bind_rows(washsf_out,sf3)
}

} #end for loop for static env covs here

# Savepoint 2 -------------------------------
#save watersheds with temporally static env covs
#saveRDS(washsf_out,file.path(objdir,"washsfout_tmp.rds"))
washsf_out=readRDS(file.path(objdir,"washsfout_tmp.rds"))

# Overlay temporal env covs -------------------------------
#env covs that do change by season
#nrow(wctrs) #40052
#for(i in 1:10){
for(i in 1:nrow(wctrs)){
if(i%%100==0){print(i)}

ctri=wctrs[i,]

#use centroids to get vals
#get temporal mean and var for diff seasons
#drt-- pull vals for all four, get mean and variance
#1-intersect ctri with each drt in list
  #make function to use in lapply for each drt.list
drt_q1_vals=lapply(drt_listq1,intersect_pull_drt)
drt_q1_mn=mean(unlist(drt_q1_vals))
drt_q1_var=var(unlist(drt_q1_vals))

drt_q2_vals=lapply(drt_listq2,intersect_pull_drt)
drt_q2_mn=mean(unlist(drt_q2_vals))
drt_q2_var=var(unlist(drt_q2_vals))

drt_q3_vals=lapply(drt_listq3,intersect_pull_drt)
drt_q3_mn=mean(unlist(drt_q3_vals))
drt_q3_var=var(unlist(drt_q3_vals))

drt_q4_vals=lapply(drt_listq4,intersect_pull_drt)
drt_q4_mn=mean(unlist(drt_q4_vals))
drt_q4_var=var(unlist(drt_q4_vals))

drt_vals=data.frame("drt_q1_mn"=drt_q1_mn,
           "drt_q1_var"=drt_q1_var,
           "drt_q2_mn"=drt_q2_mn,
           "drt_q2_var"=drt_q2_var,
           "drt_q3_mn"=drt_q3_mn,
           "drt_q3_var"=drt_q3_var,
           "drt_q4_mn"=drt_q4_mn,
           "drt_q4_var"=drt_q4_var)

#convert ctri to lat/long
ctri_4326=st_transform(ctri,crs=st_crs(4326))
ct.lon=st_coordinates(ctri_4326)[,1]
ct.lat=st_coordinates(ctri_4326)[,2]

#Need trycatch bc some watersheds missing coverage in daymetr?
#just skip these ones to get for loop to continue
tryCatch({
data=download_daymet(site="ctri",lat=ct.lat,lon=ct.lon,start=2021,end=2021,silent=TRUE)


#?tryCatch

dat=data$data
colnames(dat)
dat=dat[,c(2,3,4,7,8)]

#need to for each quarter....
#1-q1, dates 01/01-3/31
  #q2,04/01-06/30
  #q3,07/01-09/30
  #q4,10/01-12/31
q1yd=seq(from=as.integer(lubridate::yday(as.Date("2021/01/01"))),to=as.integer(lubridate::yday(as.Date("2021/03/31"))),by=1)
q2yd=seq(from=as.integer(lubridate::yday(as.Date("2021/04/01"))),to=as.integer(lubridate::yday(as.Date("2021/06/30"))),by=1)
q3yd=seq(from=as.integer(lubridate::yday(as.Date("2021/07/01"))),to=as.integer(lubridate::yday(as.Date("2021/09/30"))),by=1)
q4yd=seq(from=as.integer(lubridate::yday(as.Date("2021/10/01"))),to=as.integer(lubridate::yday(as.Date("2021/12/31"))),by=1)

#make into list to use lapply
dm_qlist=list(dat[q1yd,],dat[q2yd,],dat[q3yd,],dat[q4yd,])

#do operations using lapply and custom function
dm_vals=lapply(dm_qlist,dm_mv)

#adjust colnames to get q naming
colnames(dm_vals[[1]])=paste0(colnames(dm_vals[[1]]),"_q1")
colnames(dm_vals[[2]])=paste0(colnames(dm_vals[[2]]),"_q2")
colnames(dm_vals[[3]])=paste0(colnames(dm_vals[[3]]),"_q3")
colnames(dm_vals[[4]])=paste0(colnames(dm_vals[[4]]),"_q4")

ctri=dplyr::bind_cols(ctri,dm_vals[[1]])
ctri=dplyr::bind_cols(ctri,dm_vals[[2]])
ctri=dplyr::bind_cols(ctri,dm_vals[[3]])
ctri=dplyr::bind_cols(ctri,dm_vals[[4]])
}, error=function(e) {print(paste0("Skipping daymetr vals for watershed ", i))})

#merge all 
ctri=dplyr::bind_cols(ctri,drt_vals)

if(i==1){
  sf_tempvar=ctri
} else{
  sf_tempvar=dplyr::bind_rows(sf_tempvar,ctri)
}

} #close for loop through all centers

# Savepoint 3 -----------------------------------------------------------------

#saveRDS(sf_tempvar,file.path(objdir,"washctr_tempvar.rds"))
#sf_tempvar=readRDS(file.path(objdir,"washctr_tempvar.rds"))

# Bind static and temporal covs together -----------------------------------------------------------------

#trim to columns needed from wash ctr object
sf_tempvar2=sf_tempvar[,23:62]
sf_tempvar2=st_drop_geometry(sf_tempvar2)

#bind cols
washsf_final=dplyr::bind_cols(washsf_out,sf_tempvar2)

# Save final output ------------------------------------------------------------

saveRDS(washsf_final,file.path(objdir,"wash_envcov_final_2.rds"))



