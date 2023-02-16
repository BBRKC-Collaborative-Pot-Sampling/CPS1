library(data.table)
library(tidyverse)
library(dplyr)

sample<- list.files("./testData/EBS_Data/SAMPLE/")
sample # it reads these in!
sample_values<- list.files("./testData/EBS_Data/SAMPLE_VALUES/TEST")  # Need to edit this to pull from correct folder once done
sample_values # also read in!

catch_summary<-list.files("./testData/EBS_Data/CATCH_SUMMARY")        # Need to edit this to pull from correct folder once done
catch_summary
specimen<-list.files("./testData/EBS_Data/SPECIMEN_TABLE")           # Need to edit this to pull from correct folder once done


# create blank df to hold sample and sample_value
sample_series<- sample_value_series <- data.frame()

############# First process the SAMPLE files....these contain ancillary information 
# necessary to ID data contained in SAMPLE_VALUES files 
# need to set wd before running for loop or it won't read csvs......
# Not ideal...want to figure out alternative once code is running properly
# Add maturity status column for use with specimen table

setwd("C:/Users/jon.richar/Work/GitRepos/WinterSurvey2023/testData/EBS_Data/SAMPLE")# Need to edit this to pull from correct folder once done

for(i in 1:length(sample)){ # start i loop (each file)
  i<-1 #first file for test purposes
  #i <- 1:length(sample) #full data set
  path <- paste("./testData/EBS_Data/SAMPLE/", sample[i], sep="")                    # Need to edit this to pull from correct folder once done
  x<-read_csv(sample[i])
  x %>% select(HAUL_ID, CATCH_SAMPLE_ID,SPECIES_CODE,SEX,SAMPLE_MODIFIER) %>%
  mutate(MATURITY = if_else(SAMPLE_MODIFIER =='Immature', 0, if_else(SAMPLE_MODIFIER =='Mature', 1,if_else(SAMPLE_MODIFIER =='All Sizes', 3,0),0),0))->sample_series
  }


#sample_series
############ Now process the SAMPLE_VALUES files ################################
############ These contain the data needed to calculate sampling factors for RKC
############ and total caught for Chionoecetes spp.
# need to change wd or won't read csvs..
# Again. not ideal...want to figure out alternative once code is running properly

setwd("C:/Users/jon.richar/Work/GitRepos/WinterSurvey2023/testData/EBS_Data/SAMPLE_VALUES/TEST") # Need to edit this to pull from correct folder once done

# For-loop still needs to use KEEP to calculate sampling_factors by CATCH_SAMPLE_ID
# CATCH_SAMPLE_ID will then be used to link to data frame created by sample for-loop

for(i in 1:length(sample_values)){ # start i loop (each file)
  i <- 1 #first file for test purposes
  #i <- 1:length(sample_values)#full data set
  path <- paste("./testData/EBS_Data/SAMPLE_VALUES/", sample[i], sep="")  # Need to edit this to pull from correct folder once done
  x<-read_csv(sample_values[i])
  x %>% select(HAUL_ID, CATCH_SAMPLE_ID,COUNT,KEEP) %>% 
  mutate(KEPT_A = if_else(KEEP=='Y', COUNT,0),TOSSED_A = if_else(KEEP=='N', COUNT,0)) %>%
  select (c(-COUNT,-KEEP)) %>%
  group_by(HAUL_ID, CATCH_SAMPLE_ID) %>%
  summarize(KEPT =sum(KEPT_A),TOSSED=sum(TOSSED_A)) %>%
  mutate(TOTAL =(KEPT+TOSSED), SAMPLING_FACTOR=((KEPT+TOSSED)/(KEPT))) ->sample_values
}


sample_value_series<-sample_values[-c(1)]
sample_value_series
#merge outputs
dat<-as.data.frame(merge(sample_series, sample_value_series, by = c("CATCH_SAMPLE_ID")))    #need to combine output for correct appending to CATCH/SPECIMEN tables

dat%>%
group_by(HAUL_ID,SPECIES_CODE) %>%
summarize(NUMBER_CRAB=sum(TOTAL)) -> crab_caught                                                                         #works! YAY!Now just need to finish for-loops


# ADD CATCH_SUMMARY and update
# need to change wd or won't read csvs..
# Again. not ideal...want to figure out alternative once code is running properly

setwd("C:/Users/jon.richar/Work/GitRepos/WinterSurvey2023/testData/EBS_Data/CATCH_SUMMARY")# Need to edit this to pull from correct folder once done

for(i in 1:length(catch_summary)){ # start i loop (each file)
  i<-1 #first file for test purposes
  #i <- 1:length(catch_summary #full data set
  path <- paste("./testData/EBS_Data/CATCH_SUMMARY/", sample[i], sep="")                    # Need to edit this to pull from correct folder once done
  x<-read_csv(catch_summary[i])
  x %>% select(VESSEL,CRUISE,HAUL, SPECIES_CODE, HAUL_ID)->summary
}



catch_summary<-merge(crab_caught,summary,by=c("HAUL_ID", "SPECIES_CODE"))


# ADD SPECIMEN TABLE
# USE HAUL_ID and HAUL from CATCH SUMMARY as key to link to correct HAUL in SPECIMEN TABLE
# THEN USE SPECIES_CODE, SEX, and MATURITY (ALL FROM SAMPLE data to link to correct specimens)

catch_summary%>%select(HAUL_ID,HAUL) -> haul_key
haul_key

# need to change wd or won't read csvs..
# Again. not ideal...want to figure out alternative once code is running properly

setwd("C:/Users/jon.richar/Work/GitRepos/WinterSurvey2023/testData/EBS_Data/SPECIMEN_TABLE")# Need to edit this to pull from correct folder once done

for(i in 1:length(specimen)){ # start i loop (each file)
  i<-1 #first file for test purposes
  #i <- 1:length(specimen) #full data set
  path <- paste("./testData/EBS_Data/SPECIMEN_TABLE/", sample[i], sep="")                    # Need to edit this to pull from correct folder once done
  x<-read_csv(specimen[i])
  x %>% select(1:21,) ->x              
  merge(x,haul_key, by =c("HAUL"))->x 
  x %>% mutate(MATURITY = if_else(SEX == 2&CLUTCH_SIZE ==0, 0, if_else(SEX==2&CLUTCH_SIZE >=1, 1,if_else(SEX =='1', 3,0),0),0))->x
  specimen_out<-merge(x, dat, by=c("HAUL_ID","SPECIES_CODE","SEX","MATURITY"))
}

specimen_out # With minimal editing data are now ready for Oracle

setwd("C:/Users/jon.richar/Work/GitRepos/WinterSurvey2023/DataForOracle")# Need to edit this to pull from correct folder once done

write.csv(specimen_out, "Processed_Specimen_Data.csv")
write.csv(catch_summary, "Processed_Catch_Summary.csv")