# Original code written by JP Rippe (jpr6mg@gmail.com)
# Updated by Hannah Aichelman April 22, 2020

# This script creates the zoomed in and inset map for figure 1. These images were then put together using Adobe Illustrator.

#install and load necessary libraries
#install.packages("gpclib", type="source")
library(mapproj)
library(dplyr)
library(ggplot2)
library(ggmap)
detach("package:ggmap", unload=TRUE)
library(rgeos)
library(rgdal)
library(maps)
library(mapdata)
library(maptools)

#set working directory
setwd("~/Documents/BU/NEU_experiment/Site_Map")

#You will need to download, gunzip, and place the gshhs_f.b.gz file from NOAA in your working directory. 
#The file can be found here: https://www.ngdc.noaa.gov/mgg/shorelines/data/gshhs/oldversions/version1.2/
#clip gshhs_f.b to the broader area that you want to include in your map. gshhs_f.b is the high definition noaa coastline layer
if (!rgeosStatus()) gpclibPermit()
gshhs.f.b <- "gshhs_f.b"
#Crop global map layer to desired extent
sf1 <- getRgshhsMap(gshhs.f.b, xlim = c(-110, -55), ylim = c(0, 35)) %>%
  fortify()

#Read in coordinates of sampling sites
a=read.csv('GPS_Coordinates.csv')

#plot zoomed in map of Belize collection sites
colorz <- c("SCMR"="slategray4", "PHMR"="black")
map_zoom <- ggplot() + 
  geom_polygon(data=sf1, aes(x=long, y=lat, group = group), fill = "grey70", color='black', lwd = 0.4)+
  geom_point(data=a[c(1:3)], aes(x=site_long, y=site_lat, shape = Site), size=4.5, color=colorz) +
  scale_shape_manual(values=c(16,16))+
  xlab("Longitude")+
  ylab("Latitude")+
  theme_cowplot()+
  theme(legend.position = "none")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border=element_blank())+
  panel_border(size=1.8, colour="black",remove=FALSE)+
  coord_fixed(ratio=1, xlim = c(-89.25,-87.8), ylim = c(15.6,17))+
  scale_y_continuous(breaks=seq(15.6,17, 0.2))
map_zoom

# save ggplot figure
ggsave(file="map_zoom2.pdf", map_zoom, width = 8, height = 6, units = c("in"), useDingbats=FALSE)

#plot map further away to use for inset
map_distance <- ggplot() + 
  geom_polygon(data=sf1, aes(x=long, y=lat, group = group), fill = "grey70", color='black', lwd = 0.1)+
  #geom_point(data=a[c(1:3)], aes(x=site_long, y=site_lat, shape = Site), size=3, col = 'black') +
  #scale_shape_manual(values=c(15,17))+
  theme_bw()+
  theme(axis.title.x=element_blank(),
          #axis.text.x=element_blank(),
          #axis.ticks.x=element_blank(),
          axis.title.y=element_blank())+
          #axis.text.y=element_blank()) +
          #axis.ticks.y=element_blank())+
          #panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border=element_blank())+
  panel_border(size=1, colour="black",remove=FALSE)+
  coord_fixed(ratio=1, xlim = c(-100,-70), ylim = c(10,30))+
  scale_x_continuous(breaks=seq(-100,-70, 10))
map_distance

#save figure
ggsave(file="~/Documents/Manuscripts/AstrangiaPhysiology/Map_Figure/map_distance.pdf", map_distance, width = 5, height = 5, units = c("in"), useDingbats=FALSE)

