###Check national geolocs data to make sure formatted correctly before analysis

##########################
###### Script Setup ######
##########################
#set location for base folder to current dir
home<-"/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/StatPigMvmt/Pipeline_R2"
geo<-read.csv(paste0(home,"/2_Data/Input/Pigs/","geolocsnatl.csv"))
#geo<-read.csv("/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/StatPigMvmt/Data/Pigs/Original_Geolocs/Geolocs_Capdata/joined/geolocsnatl.csv")

#Formatting:
#remove junk cols
geo<-geo[,-1]

#set all tz to UTC
geo$tz<-"UTC"

#set classes
#study, animalid, tz, sex, age, state, lc
factors<-c(1,2,4,9,10,11,13)
#longitude, latitude, tc, mast, rgd, rds, drt, dayl, prcp, tmax, tmin
nums<-c(5,6,12,14,15,16,17,18,19,20,21)
geo[factors] <- lapply(geo[factors], as.factor)
geo[nums] <- lapply(geo[nums], as.numeric)

##############################
###### Source functions ######
##############################

source(paste0(home,"/1_Scripts/Processing/NeatDates.R"))

#date time check
dat_time_check_fn <- function(dat_time,study,dtformat) {
  require(lubridate)
  require(dplyr)
  
  #Neat dates function takes care of some issues that as.POSIXct introduces
    #otherwise weird stuff happens sometimes-- certain times get removed, and dates reformatted without timestamps
  #print("formatting datetimes with Neat.Dates.POSIXc")
  #dates1=NeatDates(dat_time)
    #sitll want as.POSIXct here to double check that all are formatted the same
  dates=as.POSIXct(dat_time,tz="UTC",format= dtformat)
  if (anyNA(dates)) stop("Error: dates must be in format '1999-12-31 23:59:59' ")
  if (any(year(dates)<1999)) stop("Error: incorrectly formatted year (year < 1999)")
  if (any(year(dates)<1999)) stop("Error: incorrectly formatted year (year < 1999)")
  stdt=data.frame("study"=study,dt=dates)
  
  #check if any studies don't have hour that goes above 12. 
  #would mean that dates with AM/PM dropped off at some pt
  stdt=stdt %>% dplyr::group_by(study) |> nest() 
  for(i in 1:nrow(stdt)){
  if(!any(hour(stdt$data[[i]]$dt)>12)) stop(paste0("Hour for study ",stdt[[1]][i]," doesn't go above 12, may be missing AM/PM in date format"))
  }
  
}

#all sanity checks (includes date time check within)
NtlGeolocSanityCheck<-function(geo.all,study,sex,age,animalid,datetime,lat,lon,tz,DOP,dtformat,epsg_utm_code){
  #no animalIDs duplicated between studies
  #convert cols to indiv vars
  require(spData)
  require(sf)
  require(amt)
  
  print("checking for duplicated animalids")
  try(if(any(duplicated(unique(geo.all[,c(study,animalid)])$animalid))) 
    stop("animalid duplicated between studies"))
  
  studies=unique(geo.all[,study])
  for(s in 1:length(studies)){
    geo=geo.all[geo.all[,which(colnames(geo.all)==study)]==studies[s],]
    print(paste0("Checking study ID: ",studies[s]))
    
  #check datetime formats
  print("checking datetime formatting")
  dat_time_check_fn(geo[,datetime],study,dtformat)
  
  #check that all lat longs fall within u.s.
  #make warning because could be a point that just falls outside of a border due to gps error
  print("checking lat long validity")
  geosf=st_as_sf(geo,coords=c(lon,lat),crs=st_crs(4326))
  us_states<-st_transform(us_states,crs=st_crs(4326))
  us_int=st_intersection(geosf,us_states)
  if(!(nrow(geosf)==nrow(us_int))) warning("Some lat longs fall outside of U.S. boundaries")
  
  #print levels of each factor incl sex, age
  print("Printing levels of select factors for review")
  print("Sex:")
  print(levels(geo[,sex]))
  print("Age:")
  print(levels(geo[,age]))
  
  #check DOPs
  if(any(geo[,DOP]>10|anyNA(geo[,DOP]))) warning("Some DOP are NA or >10")
 
  #check tz
  if(any(geo[,tz]!="UTC"|anyNA(geo[,tz]))) warning("Some tz are not UTC")
  
  #Make proofing summaries
  
   
  }
  
}

#summarize geolocs by study
summ_geo<-function(geo,dtformat,lon,lat,utm_epsg,animalid,study){
  require(tidyr)
  require(sf)
  require(plyr)
  require(dplyr)

#formatting
geo[,datetime]=as.POSIXct(geo[,datetime],tz="UTC",format= dtformat)
geosf<-st_as_sf(geo,coords=c(lon,lat),crs=st_crs(4326))
geosf<-st_transform(geosf,crs=st_crs(utm_epsg))
geo$x=st_coordinates(geosf)[,1]
geo$y=st_coordinates(geosf)[,2]
studies=unique(geo[,study])

#loop through studies
for(i in 1:length(studies)){
  print(paste0("starting study ",studies[i]))
sgeo=geo[geo[,study]==studies[i],]

#get unique IDs, loop through each animal
IDs=unique(sgeo[,animalid])

#subset to each animal
#j=1
for(j in 1:length(IDs)){
animal=sgeo[sgeo[,animalid]==IDs[j],]
animal=animal[order(animal[,datetime]),]

#make time cols for calc fixrate through each track
animal$time2=dplyr::lead(animal[,datetime])
animal$difftime=difftime(animal$time2,animal[,datetime],units="mins")

#x and y start for calc displacements
animal$x.start=animal$x[1]
animal$y.start=animal$y[2]
animal$displacement=sqrt((animal$x-animal$x.start)^2+(animal$y-animal$x.start)^2)

#lead locs for calc step lengths
animal$x2=dplyr::lead(animal$x)
animal$y2=dplyr::lead(animal$y)
animal$sl=sqrt((animal$x2-animal$x)^2+(animal$y2-animal$y)^2)

#make consolidated df with all animals in study
if(j==1&i==1){
  study.df=animal
} else{
  study.df=rbind(study.df,animal)
}

}
}

study.summary=study.df %>% 
            group_by(study,animalid) %>% 
            dplyr::summarise(numfix=n(),
                             strtdt=min(datetime),
                             enddt=max(datetime))

pig.summary=study.df %>% 
            group_by(study,animalid) %>% 
            dplyr::summarise(numfix=n(),
                             strtdt=min(datetime),
                             enddt=max(datetime))

pig.summary$dur=difftime(pig.summary$enddt,pig.summary$strtdt,units="days")

study.pig.summary=pig.summary %>% 
                  group_by(study) %>%
                  dplyr::summarise(trkdur_mean=mean(dur),
                                   trkdur_min=min(dur),
                                   trkdur_max=max(dur),
                                   numfix_mean=mean(numfix),
                                   numfix_min=min(numfix),
                                   numfix_max=max(numfix))

sex.unique=unique(study.df[,c(study,animalid,sex)])
sex.summary.l=sex.unique %>% group_by(study,sex) %>%
                dplyr::summarise(num=n())

sex.summary=tidyr::pivot_wider(sex.summary.l,names_from=sex,values_from=num)

study.summary <- study.df %>%
                 group_by(study) %>%
                 dplyr::summarise(strtdate=min(datetime),
                                  enddate=max(datetime),
                                  mlat=mean(latitude),
                                  mlon=mean(longitude))

sldisp.summary <- study.df[!is.na(study.df$sl),] %>%
                  group_by(study) %>%
                  dplyr::summarise(meanfr=mean(difftime),
                                   minfr=min(difftime),
                                   maxfr=max(difftime),
                                   meansl=mean(sl),
                                   minsl=min(sl),
                                   maxsl=max(sl),
                                   meandisp=mean(displacement),
                                   mindisp=min(displacement),
                                   maxdisp=max(displacement))

study.summary.out=join_all(list(sldisp.summary,
              study.summary,
              sex.summary,
              study.pig.summary), by='study', type='left')

return(study.summary.out)

}

#Plot geo U.S.
plot_geoUS<-function(geo,us_states,epsg){
  us_states<-st_transform(us_states,crs=st_crs(epsg))

  geo2=geo %>% group_by(study) %>% dplyr::summarise(mlat=mean(latitude),mlon=mean(longitude))
  geosf=st_as_sf(geo2,coords=c("mlon","mlat"),crs=st_crs(epsg))
  
  plot(st_geometry(us_states))
  plot(st_geometry(geosf),col="red", add=TRUE)
  
}

########################
###### Run Checks ######
########################

#arguments
animalid="animalid"
datetime="datetime"
lat="latitude"
lon="longitude"
sex="sex"
tz="tz"
age="age"
study="study"
DOP="DOP"
utm_epsg=32614
dtformat="%Y-%m-%d %H:%M:%S"
epsg=4326

NtlGeolocSanityCheck(geo,study,sex,age,animalid,datetime,lat,lon,tz,DOP,dtformat,epsg_utm_code)
geosums=summ_geo(geo,dtformat,lon,lat,utm_epsg,animalid,study)
plot_geoUS(geo,us_states)



