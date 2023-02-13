library(data.table)
library(tidyverse)
library(dplyr)

sample<- list.files("./testData/EBS_Data/SAMPLE/")
sample # it reads these in!
sample_values<- list.files("./testData/EBS_Data/SAMPLE_VALUES/TEST")
sample_values # also read in!



# create blank df to hold sample and sample_value
sample_series<- sample_value_series <- data.frame()

############# First process the SAMPLE files....these contain ancillary information 
# necessary to ID data contained in SAMPLE_VALUES files 
# need to set wd before running for loop or it won't read csvs......
# Not ideal...want to figure out alternative once code is running properly

setwd("C:/Users/jon.richar/Work/GitRepos/WinterSurvey2023/testData/EBS_Data/SAMPLE")

for(i in 1:length(sample)){ # start i loop (each file)
  i<-1 #first file for test purposes
  #i <- 1:length(sample) #full data set
  path <- paste("./testData/EBS_Data/SAMPLE/", sample[i], sep="")
  x<-read_csv(sample[i])
  x %>% select(HAUL_ID, CATCH_SAMPLE_ID,SPECIES_CODE,SEX,SAMPLE_MODIFIER)->sample_series
}

############ Now process the SAMPLE_VALUES files ################################
############ These contain the data needed to calculate sampling factors for RKC
############ and total caught for Chionoecetes spp.
#need to change wd or won't read csvs..
# Again. not ideal...want to figure out alternative once code is running properly

setwd("C:/Users/jon.richar/Work/GitRepos/WinterSurvey2023/testData/EBS_Data/SAMPLE_VALUES/TEST")

# For-loop still needs to use KEEP to calc sampling_factors by CATCH_SAMPLE_ID
# CATCH_SAMPLE_ID will then be used to link to dataframe created by sample for-loop

for(i in 1:length(sample_values)){ # start i loop (each file)
  i <- 1 #first file for test purposes
  #i <- 1:length(sample_values)#full data set
  path <- paste("./testData/EBS_Data/SAMPLE_VALUES/", sample[i], sep="")
  x<-read_csv(sample_values[i])
  x %>% select(HAUL_ID, CATCH_SAMPLE_ID,COUNT,KEEP) %>% 
  mutate(KEPT_A = if_else(KEEP=='Y', COUNT,0),TOSSED_A = if_else(KEEP=='N', COUNT,0)) %>%
  select (c(-COUNT,-KEEP)) %>%
  group_by(HAUL_ID, CATCH_SAMPLE_ID) %>%
  summarize(KEPT =sum(KEPT_A),TOSSED=sum(TOSSED_A)) %>%
  mutate(SAMPLING_FACTOR=((KEPT+TOSSED)/(KEPT))) ->sample_value_series
}


sample_value_series#






#### following is backup code for SAMPLE_VALUE for-loop that works
#for(i in 1:length(sample_values)){ # start i loop (each file)
 # i <- 1 #first file for test purposes
  #i <- 1:length(sample_values)#full data set
 # path <- paste("./testData/EBS_Data/SAMPLE_VALUES/", sample[i], sep="")
 # x<-read_csv(sample_values[i])
 # x %>% select(HAUL_ID, CATCH_SAMPLE_ID,COUNT,KEEP)->sample_value_series
#}


#merge outputs
dat<-merge(sample_series, sample_value_series, by = c("CATCH_SAMPLE_ID")) #need to combine output for correct appending to CATCH/SPECIMEN tables
dat #works! YAY!Now just need to finish for-loops

# ADD CATCH_SUMMARY and update

# ADD SPECIMEN TABLE
# USE HAUL_ID and HAUL from CATCH SUMMARY as key to link to correct HAUL in SPECIMEN TABLE
#THEN USE SPECIES_CODE, SEX, and MATURITY STATUS (ALL FROM SAMPLE data to link to correct specimens)
