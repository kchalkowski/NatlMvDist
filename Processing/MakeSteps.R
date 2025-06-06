##Calculating step lengths for national geolocation dataset
##output both daily and weekly scales
  #daily needed for natl geolocation pub
  #weekly needed for ASF modeling parameters

#input: geolocsnatl_wDispl.rds
#output: 
  #allsteps_daily.rds
  #allsteps_weekly.rds

###############
#### Setup ####
###############

#load libraries
require(amt)
require(plyr)
require(dplyr)
require(sf)
require(NatlWPGeolocv1.0.0)
require(lubridate)
require(purrr)

#set dirs
home<-"/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/StatPigMvmt/"
outdir<-"/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/StatPigMvmt/Pipeline_R2/2_Data/Objects/"

#load geolocation data
pigs=readRDS(paste0(outdir,"geolocsnatl_wDispl.rds"))

#####################################################
#### Resampling to match fixrate for all studies ####
#####################################################

#get to amt track format
pigtk=amt::mk_track(pigs,X,Y,datetime,id=animalid,all_cols=TRUE)

#get hours to separate into periods
pigtk$hour=lubridate::hour(pigtk$t_)

#take look at sampling rates
ssr=pigtk |>
  amt::summarize_sampling_rate_many(c("animalid"),time_unit="hour")%>% as.data.frame()

#viz
hist(ssr$median)
hist(ssr$min)

#nest by IDs
ptk1=pigtk |> nest(data=-"animalid")

#Resample all to 6 hours, time interval median shared by all in dataset
ptk1=ptk1 |> mutate(resampled = map(data, function(x)
                    x |> track_resample(rate = hours(6), 
                                        tolerance = hours(1))))

#####################################
#### Daily sl/disp by resampling ####
#####################################
#Resample steps to get daily-scale step lengths and displacements

#remove any bursts with lt 8 locs included (48 hrs min)
ptk1=ptk1 %>%
  mutate(filtered=map(resampled,function(x)
    x |> group_by(burst_) %>% filter(n()>=8) %>% ungroup()))

#get summary of how many remaining tracks in filtered for each group
ptk1=ptk1 %>% mutate(n=map_dbl(filtered,nrow))

#remove those with zero
ptk.f=ptk1[which(ptk1$n > 0),]

as.data.frame(ptk1$filtered[[3]])

#subset by hours and get steps
ptk.f=ptk.f |> 
  mutate(p1 = map(filtered, function(x)
    x |> filter(hour<=23&hour>=18) |>
      steps_by_burst(keep_cols='start')
  ),
  p2 = map(filtered, function(x)
    x |> filter(hour<=17&hour>=12) |>
      steps_by_burst(keep_cols='start')
  ),
  p3 = map(filtered, function(x)
    x |> filter(hour<=11&hour>=6) |>
      steps_by_burst(keep_cols='start')
  ),
  p4 = map(filtered, function(x)
    x |> filter(hour<=5&hour>=0) |>
      steps_by_burst(keep_cols='start')
  )
  )

#View to check
#as.data.frame(ptk.f$p2[[3]])

#retain only steps with at least 21-27 hours time interval
ptk.f=ptk.f |>
  mutate(p1_24 = map(p1, function(x)
    x |> filter(dt_<=27&dt_>=21)
  ),
  p2_24 = map(p2, function(x)
    x |> filter(dt_<=27&dt_>=21)
  ),
  p3_24 = map(p3, function(x)
    x |> filter(dt_<=27&dt_>=21)
  ),
  p4_24 = map(p4, function(x)
    x |> filter(dt_<=27&dt_>=21)
  )
  )

#View to check
#as.data.frame(ptk.f$p1_24[[3]])

#now need to combine p1-p4
p1=ptk.f |> dplyr::select(animalid, p1_24) |> unnest(cols = p1_24) %>% as.data.frame()
p2=ptk.f |> dplyr::select(animalid, p2_24) |> unnest(cols = p2_24) %>% as.data.frame()
p3=ptk.f |> dplyr::select(animalid, p3_24) |> unnest(cols = p3_24) %>% as.data.frame()
p4=ptk.f |> dplyr::select(animalid, p4_24) |> unnest(cols = p4_24) %>% as.data.frame()

p1$period="p1"
p2$period="p2"
p3$period="p3"
p4$period="p4"

#get together into one dataframe
allsteps_daily=rbind(p1,p2,p3,p4)

#do dummy coding for lc, to match weekly
allsteps_daily=allsteps_daily %>%
  mutate(n = 1) %>%
  tidyr::pivot_wider(names_from = lc, values_from = n, 
                     names_prefix = 'lc_', values_fill = list(n = 0)) %>%
  as.data.frame()


saveRDS(allsteps_daily,paste0(outdir,"allsteps_daily.rds"))
