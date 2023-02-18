# PURPOSE
# 1) To automate accessing the raw data folder and haul files (by pot) from the winter survey
# - Each pot is identified by the four digit number after "HAUL" and there can be potentially up to 20 pots a day
# - Each file name is appended by a time/date stamp. If multiple versions are sent of a file, the most recent one
#   needs to be read

# Load packages
library(tidyverse)
library(dplyr)

# Read in data
catch <- list.files("./Data/FTP/Catch - FTP/") %>%
  purrr::map_df(~read.csv(paste0("./Data/FTP/Catch - FTP/", .x))) 

specimen <- list.files("./Data/FTP/Specimen - FTP/") %>%
  purrr::map_df(~read.csv(paste0("./Data/FTP/Specimen - FTP/", .x)))

raw_sample <- list.files("./Data/FTP/RawData - FTP/", pattern = "_SAMPLE_0") %>% # RECORDS of SAMPLE INFO
  purrr::map_df(~read.csv(paste0("./Data/FTP/RawData - FTP/", .x))) #E.G. SEX, SPECIES

raw_sample_values <- list.files("./Data/FTP/RawData - FTP/", pattern = "_SAMPLE_VALUES") %>% #RECORDS OF # TOSSED
  purrr::map_df(~read.csv(paste0("./Data/FTP/RawData - FTP/", .x))) %>%
  mutate(TOSSED = ifelse(is.na(COUNT) == FALSE, COUNT,0)) %>%
  group_by(HAUL_ID, CATCH_SAMPLE_ID) %>%
  dplyr::summarise(TOSSED = sum(TOSSED)) 

raw_specimen <- list.files("./Data/FTP/RawData - FTP/", pattern = "_SPECIMEN_0") %>% 
  purrr::map_df(~read.csv(paste0("./Data/FTP/RawData - FTP/", .x))) %>%
  select(HAUL_ID, SPECIMEN_ID, CATCH_SAMPLE_ID, SPECIES_CODE)

raw_specimen_bio<- list.files("./Data/FTP/RawData - FTP/", pattern = "_SPECIMEN_BIOMETRICS") %>% 
  purrr::map_df(~read.csv(paste0("./Data/FTP/RawData - FTP/", .x))) 

potlifts <- list.files("./Data/FTP/Supplemental - FTP/", pattern = "POTLIFTS") %>% #MAY NEED TO CHANGE
  purrr::map_df(~read.csv(paste0("./Data/FTP/Supplemental - FTP/", .x)))

# Join raw_sample_values and raw_sample to get # tossed by sample modifier (e.g., "All sizes", "Immature")
samples <- right_join(raw_sample, raw_sample_values) %>%
  select(HAUL_ID, CATCH_SAMPLE_ID, SPECIES_CODE, SPECIES_NAME, SEX, TOSSED)

# Expand specimen biometric table, join to raw_specimen table
raw_specimen_bio %>%
  select(HAUL_ID, SPECIMEN_ID, BIOMETRIC_NAME, VALUE) %>%
  pivot_wider(., names_from = "BIOMETRIC_NAME", values_from = "VALUE") %>%
  dplyr::rename(SHELL_CONDITION = CRAB_SHELL_CONDITION, EGG_COLOR = CRAB_EGG_COLOR,
                EGG_CONDITION = CRAB_EGG_CONDITION, CLUTCH_SIZE = CRAB_EGG_CLUTCH_SIZE,
                LENGTH = CARAPACE_LENGTH) %>%
  right_join(., raw_specimen) %>%
  right_join(samples, ., by = c("HAUL_ID", "CATCH_SAMPLE_ID", "SPECIES_CODE", "SEX")) -> specimen_sum

specimen_sum %>%
  group_by(HAUL_ID, CATCH_SAMPLE_ID, SEX) %>%
  summarise(KEPT = n(),
            TOSSED = TOSSED,
            SAMPLING_FACTOR = (KEPT + TOSSED)/KEPT) %>%
  distinct() %>%
  right_join(specimen_sum, by = c("HAUL_ID", "CATCH_SAMPLE_ID", "SEX", "TOSSED")) %>%
  right_join(catch,., by = c("SPECIES_CODE", "HAUL_ID")) %>%
  dplyr::select(!c(DEVICE_RECORD_ID.x, INCLUDE_IN_KEEP_COUNT, RECORDING_DEVICE.x, RECORDING_DEVICE.y, 
                   ID, RECORDER, DATABASE_STATUS,
                   STATUS_TIMESTAMP, DEVICE_RECORD_ID.y, KEPT, TOSSED, NUMBER_CRAB, WEIGHT,
                   HAUL_ID, CATCH_SAMPLE_ID, SPECIMEN_ID, RECORD_TIMESTAMP)) %>%
  dplyr::rename(POT = HAUL) %>% 
  # Join with potlifts to get float #, lat/lon, set and pull time for each pot (MAY NEED TO CHANGE)
  right_join(potlifts, ., by = c("VESSEL", "POT"))  -> catch_summary


# Make maturity/sex and legal/sublegal categories, bind together
catch_summary %>%
  mutate(MAT_SEX = case_when((COMMON_NAME == "red king crab" & SEX == 1 & LENGTH >= 120) ~ "Mature male",
                             (COMMON_NAME == "red king crab" & SEX == 1 & LENGTH < 120) ~ "Immature male",
                             (COMMON_NAME == "red king crab" & SEX == 2 & CLUTCH_SIZE >= 1) ~ "Mature female",
                             (COMMON_NAME == "red king crab" & SEX == 2 & CLUTCH_SIZE == 0) ~ "Immature female")) -> maturity

catch_summary %>%
  mutate(MAT_SEX = case_when((COMMON_NAME == "red king crab" & SEX == 1 & LENGTH >= 135) ~ "Legal male",
                             (COMMON_NAME == "red king crab" & SEX == 1 & LENGTH < 135) ~ "Sublegal male")) %>%
  filter(is.na(MAT_SEX) == "FALSE") -> legal

# Bind to get non-overlapping classes of mat/sex and legal/sublegal  
rbind(maturity, legal) -> mat_spec


# Calculate CPUE per pot
mat_spec %>%
  # Calculate soak time, divide by sampling factor by soak time
  dplyr::mutate(SET_TIME = as.POSIXct(SET_TIME, format = "%m/%d/%Y %H:%M"),
                PULL_TIME = as.POSIXct(PULL_TIME, format = "%m/%d/%Y %H:%M"),
                SOAK_TIME = as.numeric(difftime(PULL_TIME, SET_TIME, units = "mins")),
                CPUE = SAMPLING_FACTOR/SOAK_TIME) %>%
  dplyr::group_by(VESSEL, POT, FLOAT, LATITUDE, LONGITUDE, MAT_SEX, COMMON_NAME) %>%
  dplyr::summarise(COUNT = sum(SAMPLING_FACTOR),
                   CPUE = sum(CPUE)) -> pot_cpue












tt[4,] <- tt[3,]
tt[4,]$RECORD_TIMESTAMP <- yy[10,]$RECORD_TIMESTAMP

tt%>%
  group_by(HAUL_ID, CATCH_SAMPLE_ID) %>%
  mutate(dup = ifelse(n() > 1, "Y", "N")) %>%
  ungroup() %>%
  nest(data = !c(HAUL_ID, CATCH_SAMPLE_ID, RECORD_TIMESTAMP, dup)) %>%
  filter(dup == "Y", max(RECORD_TIMESTAMP)) -> tt
mutate(data = ifelse(dup == "Y", filter(max(RECORD_TIMESTAMP)), data)) -> ll


group_by(HAUL_ID, CATCH_SAMPLE_ID) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  unnest()

%>%
  nest() %>%
  mutate(n = n()) %>%
  unnest() -> pp

pp %>%
  mutate(data = ifelse(duplicated(CATCH_SAMPLE_ID) == TRUE, data %>% filter(max(RECORD_TIMESTAMP)), data))


%>%
  unnest()



ifelse(duplicated(CATCH_SAMPLE_ID) == TRUE, filter(CATCH_SAMPLE_ID,)
       
       
       tt%>%
         group_by(HAUL_ID, CATCH_SAMPLE_ID) %>%
         nest() %>%
         summarise(n = n()) %>%
         unnest()
       
       
       list.files("./Data/RawData", pattern = "SAMPLE_VALUES") %>%
         purrr::map_df(~read.csv(paste0("./Data/RawData/", .x))) %>%
         mutate(KEPT_A = ifelse(KEEP=='Y', COUNT,0),
                TOSSED_A = ifelse(KEEP=='N', COUNT,0)) %>%
         group_by(HAUL_ID, CATCH_SAMPLE_ID) %>%
         summarize(KEPT =sum(KEPT_A),
                   TOSSED=sum(TOSSED_A)) %>%
         mutate(SAMPLING_FACTOR=((KEPT+TOSSED)/(KEPT))) -> sample_value_series
       
       