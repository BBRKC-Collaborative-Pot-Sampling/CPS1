##############################################################################################
### 2023 BBRKC Cooperative Pot Survey 

# Read ABSC and ADFG T/S logger data for plotting and database export

# Author: May 2023 - jdweems

### 1) Load and summarise ABSC logger data
### 2) Load and summerise ADFG logger data
### 3) Load Master Stations grid
### 4) Sampling and EDA Maps
### 5) Combine data sets, export data to .csv, and revisualize T/S

##############################################################################################

# set parent working directory to extracted zip file
setwd("V:/BSAI King Crab/2023BBRKC/DATA/Oceanography/2023_BBRKC_ALL_TEMPS_SAL")
#rm(list = ls())
getwd()

# libraries
library(tidyverse)
library(na.tools)
library(ggplot2)
library(fuzzyjoin)
library(PBSmapping) # `PBSmapping` has a better defined coastline than the `worldHires` database. 
data('nepacLLhigh') # Alaska based upon PBSmapping data.
library(viridis)



### 1) Load and summarise ABSC logger data

# import ABSC data
absc <- read.csv('./AMCC_FV_NRT_d113_5950_44d0.csv', skip = 2, header = FALSE)
header <- scan('./AMCC_FV_NRT_d113_5950_44d0.csv', nlines = 1, what = character())
header <- as.list(strsplit(header, ",")[[1]])
names(absc) <- header
head(absc); str(absc); # view(absc)

# sort and filter data for haul average bottom depth and temperature
absc.2 <- absc %>% 
  filter(segment_type == 'Fishing') %>% 
  group_by(vessel_id, tow_id) %>% 
  summarise(Latitude = mean(latitude),
            Longitude = mean(longitude),
            AveDepth = mean(depth),
            AveTemp = mean(temperature),
            SDTemp = sd(temperature)
)
head(absc.2); str(absc.2); # view(absc.2)



### 2) Load and summerise ADFG logger data

# import ADFG data
adfg <- read.csv('./CPS1_2023BB_RKC_LoggerData_ADFG.csv')
head(adfg); str(adfg); # view(adfg)

# sort and filter data for haul average bottom depth and temperature
adfg.2 <- adfg %>% 
  group_by(Vessel, SPN) %>% 
  summarise(Latitude = mean(Latitude),
            Longitude = mean(Longitude),
            Logger = mean(Logger),
            AveDepth = mean(Depth_m),
            AveTemp = mean(Temperature_C., na.rm = T),
            SDTemp = sd(Temperature_C., na.rm = T),
            AveSal = mean(Salinity_PSU, na.rm = T),
            SDSal = sd(Salinity_PSU, na.rm = T)
  )
head(adfg.2); str(adfg.2); # view(adfg.2)



### 3) Load Master Stations grid

# import master grid stations and locations from Leah's survey design
# jdweems - added lat lon decimal degree conversions to spreadsheet
stn <- read.csv('./Pot_IDs_and_Coordinates.csv') # saved as .csv UTF-8
head(stn); str(stn); # view(absc)
stn.2 <- stn[ ,c("Pot_ID", "Lat_DD", "Lon_DD")]
stn.2 <- rename(stn.2, c(Latitude = Lat_DD))
stn.2 <- rename(stn.2, c(Longitude = Lon_DD))
head(stn.2); str(stn.2); # view(stn.2)



### 3) Join stn level data streams by Pot_Id

# join stn to adfg
adfg.3 <- adfg.2 %>% fuzzyjoin::geo_inner_join(stn.2, max_dist = 1)
str(adfg.3); names(adfg.3)
adfg.3 <- adfg.3[,c(11,12,13,1,5,3,4,6:10)]

adfg.3 <- rename(adfg.3, c(Latitude_grid = Latitude.y))
adfg.3 <- rename(adfg.3, c(Longitude_grid = Longitude.y))
adfg.3 <- rename(adfg.3, c(Latitude_pot = Latitude.x))
adfg.3 <- rename(adfg.3, c(Longitude_pot = Longitude.x))

str(adfg.3); names(adfg.3); # view(adfg.3)

# join stn to absc
absc.3 <- absc.2 %>% fuzzyjoin::geo_inner_join(stn.2, max_dist = 1)
str(absc.3); names(absc.3)
absc.3 <- absc.3[,c(8,9,10,1:7)]

absc.3 <- rename(absc.3, c(Latitude_grid = Latitude.y))
absc.3 <- rename(absc.3, c(Longitude_grid = Longitude.y))
absc.3 <- rename(absc.3, c(Latitude_pot = Latitude.x))
absc.3 <- rename(absc.3, c(Longitude_pot = Longitude.x))

absc.3 <- rename(absc.3, c(Vessel = vessel_id))
absc.3 <- rename(absc.3, c(Logger = tow_id))
absc.3$Vessel[absc.3$Vessel == "120"] <-"Silver Spray"
absc.3$Vessel[absc.3$Vessel == "119"] <-"Summer Bay"
absc.3$AveSal <- NA
absc.3$SDSal <- NA
str(absc.3); names(absc.3);

adfg.3$Source <- "ADFG"
absc.3$Source <- "ABSC"

dim(adfg.3); dim(absc.3)
names(adfg.3); names(absc.3)
#view(adfg.3)
#view(absc.3)



### 4) Sampling and EDA Maps

# basemap
iso <- ggplot() + 
  geom_polygon(data=nepacLLhigh,aes(x=X,y=Y,group=PID), fill=8, color='black') +
  coord_map(xlim = c(-165,-159),ylim = c(54.75,58)) +
  # coord_quickmap(xlim = c(-170, -160),ylim = c(62.75,65.62)) + # doesn't work with Seth Bathymetry
  xlab(expression(paste(Longitude^o,~'W'))) +
  ylab(expression(paste(Latitude^o,~'N'))) +
  theme_bw() +
  theme(panel.background = element_rect(fill = 'white'), 
        panel.grid.major = element_line(color = 'white')) #+
iso

# red king crab savings area feature
rkcsa <- data.frame(x = c(-162, -162, -164, -164),
                    y = c(56, 57, 57, 56))
rkcssa <- data.frame(x = c(-162, -162, -164, -164),
                     y = c(56, 56.1666, 56.1666, 56))
m2.mooring <- data.frame(x = c(-164.06), y = c(56.87), z = c("M2"))

# adfg.3 map
iso.adfg.3 <- iso + 
  # annotations
  geom_polygon(rkcsa, mapping=aes(x, y), color = "green", linewidth=1, fill = NA) + 
  geom_polygon(rkcssa, mapping=aes(x, y), color = "green", linewidth=1, fill = NA, linetype = 2) +
  geom_text(data = m2.mooring, mapping=aes(x, y, label = z), colour="blue", size=1.75) +
  annotate(geom="text", x = -165, y = 58, label = "2023 BBRKC Winter/Spring Pot Survey", size=3, hjust = 0) +
  # pot locations
  geom_text(data = adfg.3, mapping=aes(Longitude_pot, Latitude_pot, 
                                        label = Pot_ID, color = AveTemp), size=1.75) +
  #geom_text(data = absc.3, mapping=aes(Longitude_pot, Latitude_pot, 
  #                                      label = Pot_ID, color = AveTemp), size=1.75) +
  scale_colour_viridis_c(begin = 0, end = 0.85, direction = 1, option = "A")
iso.adfg.3

# absc.3 map
iso.absc.3 <- iso + 
  # annotations
  geom_polygon(rkcsa, mapping=aes(x, y), color = "green", linewidth=1, fill = NA) + 
  geom_polygon(rkcssa, mapping=aes(x, y), color = "green", linewidth=1, fill = NA, linetype = 2) +
  geom_text(data = m2.mooring, mapping=aes(x, y, label = z), colour="blue", size=1.75) +
  annotate(geom="text", x = -165, y = 58, label = "2023 BBRKC Winter/Spring Pot Survey", size=3, hjust = 0) +
  # pot locations
  #geom_text(data = adfg.3, mapping=aes(Longitude_pot, Latitude_pot, 
  #                                     label = Pot_ID, color = AveTemp), size=1.75) +
  geom_text(data = absc.3, mapping=aes(Longitude_pot, Latitude_pot, 
                                        label = Pot_ID, color = AveTemp), size=1.75) +
  scale_colour_viridis_c(begin = 0, end = 0.85, direction = 1, option = "A")
iso.absc.3



### 5) Combine data sets, export data to .csv, and revisualize T/S

# bind rows with common headers
names(adfg.3) == names(absc.3)
dat <- rbind(adfg.3, absc.3)
dat$AveSal[is.nan(dat$AveSal)] <- NA

dim(adfg.3); dim(absc.3); dim(dat)

# write csv of combined data
#write.csv(dat, "./2023_BBRKC_ALL_TEMPS_SAL.csv", row.names=FALSE)

# all temps by station
iso.all.temps <- iso + 
  # annotations
  geom_polygon(rkcsa, mapping=aes(x, y), color = "green", linewidth=1, fill = NA) + 
  geom_polygon(rkcssa, mapping=aes(x, y), color = "green", linewidth=1, fill = NA, linetype = 2) +
  geom_text(data = m2.mooring, mapping=aes(x, y, label = z), colour="blue", size=1.75) +
  annotate(geom="text", x = -165, y = 58, label = "2023 BBRKC Winter/Spring Pot Survey", size=3, hjust = 0) +
  # pot locations
  geom_text(data = dat, mapping=aes(Longitude_pot, Latitude_pot, 
                                       label = Pot_ID, color = AveTemp), size=1.75) +
  scale_colour_viridis_c(begin = 0, end = 0.85, direction = 1, option = "A")
iso.all.temps

iso.all.temps.tile <- iso + 
  # annotations
  geom_polygon(rkcsa, mapping=aes(x, y), color = "green", linewidth=1, fill = NA) + 
  geom_polygon(rkcssa, mapping=aes(x, y), color = "green", linewidth=1, fill = NA, linetype = 2) +
  geom_text(data = m2.mooring, mapping=aes(x, y, label = z), colour="blue", size=1.75) +
  annotate(geom="text", x = -165, y = 58, label = "2023 BBRKC Winter/Spring Pot Survey", size=3, hjust = 0) +
  # pot locations
  geom_tile(data = dat, mapping=aes(Longitude_pot, Latitude_pot, fill = AveTemp), width=0.15, height=0.15) + 
  stat_contour() +
  scale_colour_viridis_c(begin = 0, end = 0.85, direction = 1, option = "A") +
  scale_fill_viridis_c(begin = 0, end = 0.85, direction = 1, option = "A")
iso.all.temps.tile


# all salinities by station
iso.all.sal <- iso + 
  # annotations
  geom_polygon(rkcsa, mapping=aes(x, y), color = "green", linewidth=1, fill = NA) + 
  geom_polygon(rkcssa, mapping=aes(x, y), color = "green", linewidth=1, fill = NA, linetype = 2) +
  geom_text(data = m2.mooring, mapping=aes(x, y, label = z), colour="blue", size=1.75) +
  annotate(geom="text", x = -165, y = 58, label = "2023 BBRKC Winter/Spring Pot Survey", size=3, hjust = 0) +
  # pot locations
  geom_text(data = na.omit(dat), mapping=aes(Longitude_pot, Latitude_pot, 
                                    label = Pot_ID, color = AveSal), size=1.75) +
  scale_colour_viridis_c(begin = 0, end = 0.85, direction = 1, option = "A")
iso.all.sal

iso.all.sal.tile <- iso + 
  # annotations
  geom_polygon(rkcsa, mapping=aes(x, y), color = "green", linewidth=1, fill = NA) + 
  geom_polygon(rkcssa, mapping=aes(x, y), color = "green", linewidth=1, fill = NA, linetype = 2) +
  geom_text(data = m2.mooring, mapping=aes(x, y, label = z), colour="blue", size=1.75) +
  annotate(geom="text", x = -165, y = 58, label = "2023 BBRKC Winter/Spring Pot Survey", size=3, hjust = 0) +
  # pot locations
  geom_tile(data = na.omit(dat), mapping=aes(Longitude_pot, Latitude_pot, fill = AveSal), width=0.15, height=0.15) + 
  stat_contour() +
  scale_colour_viridis_c(begin = 0, end = 0.85, direction = 1, option = "A") +
  scale_fill_viridis_c(begin = 0, end = 0.85, direction = 1, option = "A")
iso.all.sal.tile

# too save any plot...
# set parent working directory to Data Folder
# getwd()
# setwd("V:/BSAI King Crab/2023 BB RKC Tag/Survey_Design_Leah/Z- Jared")
# ggsave("XXX.jpg", dpi = 600, width=8, units="in")
# ggsave("./XXX.tiff", device = "tiff", plot=XXX, dpi = 600, width=8, units="in")


### End of Code ###