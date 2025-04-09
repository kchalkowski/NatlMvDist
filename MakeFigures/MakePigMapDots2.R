#The purpose of this script is to make a U.S. map with dots indicating study sites
#size the dots for n pigs and n days

#################
##Set up script##
#################
library("maps")
library(ggforce)
library(ggplot2)
library(tidyverse)
library(dismo)
library(gbm)
library(cowplot)
library(pdp)
library(gridExtra)
library(plyr)
library(dplyr)
library(usmap)
library(ggrepel)
library(sf)
#install.packages("usmap")
library("rnaturalearth")
library("rnaturalearthdata")


home="/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/StatPigMvmt"
outdir="/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/StatPigMvmt/Outputs/"

#read needed tables
studydf<-read.csv(paste0(outdir,"YProcessingObjects/studydf.csv"))
studydf<-sdf3
sf_use_s2(FALSE)

#aggregate locations
avglocs<-as.data.frame(pigswsite) %>% 
                         group_by(study) %>% 
                         dplyr::summarise(mean.lat.x=mean(latitude),
                                          mean.lon.x=mean(longitude))

#reproject to map onto us_map
geolocations<-data.frame(Longitude=avglocs$mean.lon.x,Latitude=avglocs$mean.lat.x)
geoloc=usmap_transform(geolocations,input_names=c("Longitude","Latitude"))
avglocs=cbind(avglocs,geoloc)
avglocs=avglocs[order(-avglocs$num.pigs),]
avglocs=left_join(avglocs,studydf,by="study")

#Make overall U.S. map for pig locations
pos <- position_jitter(width=0.5, seed=1)
plot_usmap(regions = "states",color="black") +
  theme_minimal()+theme(panel.grid.major=element_blank(),
                        panel.grid.minor=element_blank(),
                        panel.background=element_blank(),
                        axis.text=element_blank(),
                        axis.title=element_blank(),
                        legend.position="top")+
  #geom_point(data=avglocs,aes(x=x,y=y,shape=region,color=region))+
  geom_sf(data=geoloc)+
  scale_shape_manual(values=seq(1,25,1))+
  scale_color_manual(values=seq(1,25,1))


world<-ne_countries(scale="medium",returnclass="sf")
world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)

states <- st_as_sf(map("state", plot = FALSE, fill = TRUE))
head(states)
states <- cbind(states, st_coordinates(st_centroid(states)))

#Make three zoomed-in maps: Michigan, California, Southeast,TX
#make subsets
avglocs.Big=avglocs[avglocs$state!="SC"&avglocs$state!="FL",]
avglocs.SC=avglocs[avglocs$state=="SC",]
avglocs.CA=avglocs[avglocs$state=="CA",]
avglocs.FL=avglocs[avglocs$state=="FL",]


####Plots

  #Plot CA points
ggplot(data = world) +
  geom_sf(color="white",fill=NA,stroke=10) +
  geom_sf(data = states, color="#a9ebcd", fill = "black",linewidth=0.2) + 
  theme_minimal()+theme(panel.grid.major=element_blank(),
                        panel.grid.minor=element_blank(),
                        panel.background=element_blank(),
                        axis.title=element_blank(),
                        legend.position="left")+
  coord_sf(xlim = c(-124, -117), ylim = c(33.5, 39.5), expand = FALSE)+
  geom_point(data=avglocs.CA,aes(x=mean.lon.x,y=mean.lat.x),
             pch=21,size=7,stroke=1,colour="#000000",fill="#cafc03")+
  scale_size_continuous(range = c(3,30))+
  theme(text = element_text(size = 28))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0, hjust=0))


#Plot SC points
ggplot(data = world) +
  geom_sf(color="white",fill="black",stroke=10) +
  geom_sf(data = states, color="#a9ebcd", fill = "black",linewidth=0.2) + 
  theme_minimal()+theme(panel.grid.major=element_blank(),
                        panel.grid.minor=element_blank(),
                        panel.background=element_blank(),
                        axis.title=element_blank())+
  theme(text = element_text(size = 28,family="sans"))+
  #geom_sf(data = sites, size = 4, shape = 23, fill = "darkred") +
  coord_sf(xlim = c(-81.85, -81.5), ylim = c(33.35, 33.15), expand = FALSE)+
  geom_point(data=avglocs.SC,aes(x=mean.lon.x,y=mean.lat.x,fill=study),
             pch=21,size=7,stroke=1,colour="#000000")+
  scale_size_continuous(range = c(1,10))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0, hjust=0))
  

#Plot FL points
ggplot(data = world) +
  geom_sf(color="white",fill="black",stroke=10) +
  geom_sf(data = states, color="#a9ebcd", fill = "black",linewidth=0.2) + 
  theme_minimal()+theme(panel.grid.major=element_blank(),
                        panel.grid.minor=element_blank(),
                        panel.background=element_blank(),
                        #axis.text=element_blank(),
                        axis.title=element_blank())+
  theme(text = element_text(size = 28))+
coord_sf(xlim = c(-81.3, -81), ylim = c(27.05, 27.25), expand = FALSE)+
  geom_point(data=avglocs.FL,aes(x=mean.lon.x,y=mean.lat.x,fill=study),
             pch=21,size=7,stroke=1,colour="#000000")+
  scale_size_continuous(range = c(1,10))+
  theme(text = element_text(size = 28,family="sans"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0, hjust=0))


#plot US
ggplot() +
  geom_sf(color="white",fill="black",stroke=10) +
  geom_sf(data = states, color="#a9ebcd", fill = "black",linewidth=0.2) + 
  theme_minimal()+theme(panel.grid.major=element_blank(),
                        panel.grid.minor=element_blank(),
                        panel.background=element_blank(),
                        axis.title=element_blank(),
                        legend.position="left")+
  geom_point(data=avglocs.Big,aes(x=mean.lon.x,y=mean.lat.x,fill=study),
             pch=21,size=5,stroke=1,colour="#000000")+
  scale_size_continuous(range = c(3,30))+
  theme(text = element_text(size = 28,family="sans"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0, hjust=0))

sdf3$study
