# PURPOSE -------------------------------------------------------------------------------------------------------------------
# 1) To produce figures and analyses for CPS1 technical report

# Author: Emily Ryznar, NOAA-AFSC

# LOAD PACKAGES -------------------------------------------------------------------------------------------------------------
library(tidyverse)
library(gsubfn)
library(terra)
library(rgdal)
library(colorRamps)
library(sf)
library(viridis)
library(grid)
library(shadowtext)
library(gstat)
library(ggnewscale)
library(purrr)
library(ggpubr)
library(ggh4x)
#remotes::install_github("coolbutuseless/ggpattern")
library("ggpattern")

# LOAD DATA -----------------------------------------------------------------------------------------------------------------
  # Run CPS1 data processing script to get output
    source("./Scripts/CPS1_2023_processing.R")  

  # NMFS EBS survey data for RKC ---
    read.csv("./Data/crabhaul_rkc.csv") %>%
      mutate(START_DATE = lubridate::mdy(START_DATE),
             AKFIN_SURVEY_YEAR = lubridate::year(START_DATE),
             WIDTH_1MM = NA,
             SPECIES_NAME = "Red King Crab") -> catch_rkc #rkc catch
  
    read.csv("./Data/strata_rkc_newtimeseries.csv") -> strata_rkc #rkc strata
    
    read.csv("./Data/retow_sta.csv") -> retow_sta #retow stations by year for mature rkc females
    
    read.csv("./Data/haul_newtimeseries.csv") -> all_hauls
  
# REQUEST 1: GENERATE BBRKC NMFS SURVEY ABUNDANCE TS -------------------------------------------------------------------------------------
  # Load functions
  bio_abund <- function(data_crab, data_strata, stock, years){
    
    #District area lookup table by stock, pull districts by stock specified
    data.frame(STOCK = c("BBRKC", 
                         rep("PribRKC", 2),
                         "NSRKC",
                         rep("PribBKC", 2), 
                         rep("StMattBKC", 2),
                         "BKCNBS",
                         rep("TannerW", 3),
                         "TannerWNBS",
                         "TannerE", 
                         "TannerENBS",
                         rep("Snow", 3), 
                         "SnowNBS",
                         rep("Hybrid", 3), 
                         "HybridNBS",
                         rep("Hair", 4), 
                         "HairNBS",
                         rep("Allstations", 10),
                         rep("NorthRKC",1)),
               DISTRICT = c("Bristol Bay", 
                            "Pribilof MTCA", "Pribilof Single",
                            "Norton Sound",
                            "Pribilof MTCA", "Pribilof Single",
                            "St. Matthew MTCA", "St. Matthew Single",
                            "NBS All",
                            "Pribilof MTCA", "St. Matthew MTCA", "West 166",
                            "NBS All",
                            "East 166",
                            "NBS All",
                            "Pribilof MTCA", "Single", "St. Matthew MTCA",
                            "NBS All",
                            "Pribilof MTCA", "Single", "St. Matthew MTCA",
                            "NBS All",
                            "Bristol Bay", "Northern Unstratified", "Pribilof MTCA", "Pribilof Single",
                            "NBS All",
                            "Bristol Bay", "Northern Unstratified", "Pribilof MTCA","Pribilof Single",      
                            "BKC Unstratified", "St. Matthew MTCA", "St. Matthew Single", "East 166",             
                            "West 166", "Single",("Northern Unstratified"))) %>%
      dplyr::filter(STOCK == stock) %>%
      pull(DISTRICT) -> dist
    
    #Pull stock stations from strata tables using stock districts
    data_strata %>%
      dplyr::filter(DISTRICT %in% dist) %>%
      pull(STATION_ID) %>%
      unique() -> stock_stations
    
    #Specify retow stations for BBRKC, pull by year
    retow_sta <- read.csv("./Data/retow_sta.csv") ### NEED TO ADD RETOW STATIONS EACH YEAR
    
    retow_sta %>%
      dplyr::filter(SURVEY_YEAR %in% years) %>%
      dplyr::select(GIS_STATION) %>%
      pull() -> retow_stations
    
    # NBS stocks that have same mat/sex criteria as EBS stocks
    if(stock == "TannerENBS"){
      stock <- "TannerE"
    } else if(stock == "TannerWNBS"){
      stock <- "TannerW"
    } else if(stock == "HybridNBS"){
      stock <- "Hybrid"
    } else if(stock == "HairNBS"){
      stock <- "Hair"
    } else{
      stock <- stock
    }
    
    #Create lookup table for male cutlines 
    mat_lookup <- tibble(stock = c("BBRKC", "PribRKC", "PribBKC", "StMattBKC", "TannerE", "TannerW","Snow", "Hybrid", "Hair", "NorthRKC","NSRKC", "BKCNBS", "SnowNBS"),
                         cutline = c(120, 120, 120, 105, 113, 103, 95, 95, NA, 120, 94, 105, 68),
                         legal = c(135, 135, 135, 120, 120, 110, 78, 78, 83, 135,104, 120, 78),
                         recruit = c(134, 134, 134, 119, 124, 124, 101, 101, NA, 134, NA,NA,NA))
    
    #Create maturity levels
    if (stock != "Hair"){
      data_crab %>%
        dplyr::filter(SEX %in% 1:2,
                      AKFIN_SURVEY_YEAR %in% years) %>%
        # dplyr::filter(SHELL_CONDITION == 2)%>%   ######Use if only want SC 2 (e.g. for Tanner Snow tables comparing with chela maturity)
        mutate(MAT_SEX = case_when((SEX == 1 & stock %in% c("BBRKC", "PribRKC", "PribBKC", "StMattBKC","NorthRKC","NSRKC","BKCNBS") & LENGTH_1MM >= mat_lookup$cutline[mat_lookup$stock == stock]) ~ "Mature Male",
                                   (SEX == 1 & stock %in% c("BBRKC", "PribRKC", "PribBKC", "StMattBKC","NorthRKC","NSRKC","BKCNBS") & LENGTH_1MM < mat_lookup$cutline[mat_lookup$stock == stock]) ~ "Immature Male",
                                   (SEX == 1 & stock %in% c("TannerE", "TannerW", "Snow", "SnowNBS", "Hybrid") & WIDTH_1MM >= mat_lookup$cutline[mat_lookup$stock == stock]) ~ "Mature Male",
                                   (SEX == 1 & stock %in% c("TannerE", "TannerW", "Snow", "SnowNBS", "Hybrid") & WIDTH_1MM < mat_lookup$cutline[mat_lookup$stock == stock]) ~ "Immature Male",
                                   (SEX == 2 & CLUTCH_SIZE >= 1) ~ "Mature Female",
                                   (SEX == 2 & CLUTCH_SIZE == 0) ~ "Immature Female")) -> mature #for mature/immature males and females
      data_crab %>%  
        dplyr::filter(SEX %in% 1, 
                      AKFIN_SURVEY_YEAR %in% years) %>%
        mutate(MAT_SEX = case_when((stock %in% c("BBRKC", "PribRKC", "PribBKC", "StMattBKC","NorthRKC","NSRKC","BKCNBS") & LENGTH_1MM >= mat_lookup$legal[mat_lookup$stock == stock]) ~ "Legal Male",
                                   (stock %in% c("BBRKC", "PribRKC", "PribBKC", "StMattBKC","NorthRKC","NSRKC","BKCNBS") & LENGTH_1MM <= mat_lookup$recruit[mat_lookup$stock == stock]) ~ "Pre-recruit Male", 
                                   (stock %in% c("TannerE", "TannerW", "Snow", "SnowNBS", "Hybrid") & WIDTH_1MM >= mat_lookup$legal[mat_lookup$stock == stock]) ~ "Legal Male",
                                   (stock %in% c("TannerE", "TannerW", "Snow", "SnowNBS", "Hybrid") & WIDTH_1MM <= mat_lookup$recruit[mat_lookup$stock == stock]) ~ "Pre-recruit Male")) -> leg
      
      data_crab %>%  
        dplyr::filter(SEX %in% 1, 
                      AKFIN_SURVEY_YEAR %in% years) %>%
        mutate(MAT_SEX = case_when((stock == "TannerE" & WIDTH_1MM >= 125) ~ "Industry Preferred",
                                   (stock == "TannerW" & WIDTH_1MM >= 125) ~ "Industry Preferred",
                                   #Snow Crab industry preferred
                                   (stock == "Snow" & WIDTH_1MM >= 102) ~ "Industry Preferred",
                                   (stock == "SnowNBS" & WIDTH_1MM >= 102) ~ "Industry Preferred")) -> ind_pref
      
      rbind(mature, leg, ind_pref) -> mat_spec 
      
    } else{
      data_crab %>%
        dplyr::filter(SEX %in% 1:2,
                      AKFIN_SURVEY_YEAR %in% years) %>%
        mutate(MAT_SEX = case_when((SEX == 1 & LENGTH_1MM >= mat_lookup$legal[mat_lookup$stock == stock]) ~ "Legal Male",
                                   (SEX == 1 & LENGTH_1MM < mat_lookup$legal[mat_lookup$stock == stock]) ~ "Sublegal Male",
                                   (SEX == 2) ~ "Female")) -> mat_spec
    }
    
    #Calculate CPUE by GIS STATION, maturity, and year
    mat_spec %>%
      dplyr::filter(GIS_STATION %in% stock_stations) %>%
      mutate(COUNT = SAMPLING_FACTOR,
             CPUE = SAMPLING_FACTOR/AREA_SWEPT,
             CPUE_KG = (SAMPLING_FACTOR * CALCULATED_WEIGHT_1MM) / AREA_SWEPT / 1000)%>%
      group_by(AKFIN_SURVEY_YEAR, HAUL_TYPE, GIS_STATION, MAT_SEX)%>%
      dplyr::summarise(COUNT = sum(COUNT),
                       CPUE = sum(CPUE),
                       CPUE_KG =sum(CPUE_KG)) %>%
      filter(is.na(MAT_SEX) == FALSE) -> cpue
    
    #Pull out unique haul types by year in order to include in zero-catch station dataframe below
    data_crab%>%
      dplyr::filter(AKFIN_SURVEY_YEAR %in% years) %>%
      distinct(HAUL_TYPE)->HT
    
    #Conditionally specifying maturity/sex combos for each stock to expand grid below
    if (stock %in% c("BBRKC", "PribRKC", "PribBKC", "StMattBKC", "Hybrid", "NorthRKC", "NSRKC","BKCNBS")){
      
      mat_sex_combos <- c("Mature Male", "Immature Male", "Mature Female", 
                          "Immature Female", "Legal Male", "Pre-recruit Male")
      
    } else if (stock %in% c("TannerE", "TannerW", "Snow", "SnowNBS")) {
      
      mat_sex_combos <- c("Mature Male", "Immature Male", "Mature Female", 
                          "Immature Female", "Legal Male", "Industry Preferred")
    } else{
      mat_sex_combos <- c("Sublegal Male", "Legal Male", "Female")
    }
    
    
    #Join to zero catch stations, summarize 
    cpue %>%
      right_join(expand_grid(MAT_SEX = mat_sex_combos,
                             HAUL_TYPE = as.numeric(levels(factor(HT$HAUL_TYPE))),
                             data_strata %>%
                               dplyr::filter(SURVEY_YEAR %in% years,
                                             STATION_ID %in% stock_stations) %>% 
                               distinct(SURVEY_YEAR, STATION_ID, STRATUM, TOTAL_AREA) %>%
                               dplyr::rename(c("GIS_STATION"="STATION_ID", "STRATUM" = "STRATUM", 
                                               "AKFIN_SURVEY_YEAR"="SURVEY_YEAR", "TOTAL_AREA" = "TOTAL_AREA")))) %>%
      replace_na(list(COUNT = 0, CPUE = 0, CPUE_KG=0))%>%
      group_by(AKFIN_SURVEY_YEAR, GIS_STATION, MAT_SEX) %>%
      nest() %>%
      #Females: replacing original stations with resampled stations in retow yrs for BBRKC females
      mutate(data = purrr::map2(data, MAT_SEX, function(data, mat_sex) {
        if(17 %in% data$HAUL_TYPE & stock == "BBRKC" & (mat_sex == "Mature Female"|mat_sex == "Immature Female")
           & GIS_STATION %in% retow_stations) 
        {data %>% dplyr::filter(HAUL_TYPE == 17) -> x} else{x <- data %>% dplyr::filter(HAUL_TYPE != 17)}
        return(x)  
      })) %>%
      unnest() -> station_haul_cpue
    
    #Filtering out GIS_STATION E-11 in year 2000 for BBRKC males because it wasn't sampled in leg 1
    if(stock=="BBRKC"){
      station_haul_cpue %>%
        dplyr::filter(!(AKFIN_SURVEY_YEAR==2000 & GIS_STATION=="E-11" & grepl("Male", MAT_SEX))) -> station_haul_cpue
    } else{
      station_haul_cpue = station_haul_cpue
    }
    
    #Sum across haul, scale abundance, biomass, and variance to strata, then sum across strata and calc CIs
    station_haul_cpue %>%
      group_by(AKFIN_SURVEY_YEAR, GIS_STATION, STRATUM, MAT_SEX, TOTAL_AREA) %>% 
      dplyr::summarise(COUNT = sum(COUNT), CPUE = sum(CPUE), CPUE_KG = sum(CPUE_KG))%>%
      #Scale to abundance by strata
      group_by(AKFIN_SURVEY_YEAR, STRATUM, MAT_SEX) %>%
      dplyr::reframe(AREA = TOTAL_AREA,
                     MEAN_CPUE = mean(CPUE),
                     N_CPUE = n(),
                     VAR_CPUE = (var(CPUE)*(AREA^2))/N_CPUE,
                     MEAN_CPUE_KG = mean(CPUE_KG),
                     N_CPUE_KG = n(),
                     VAR_CPUE_KG = (var(CPUE_KG)*(AREA^2))/N_CPUE_KG,
                     ABUNDANCE = (MEAN_CPUE * AREA),
                     BIOMASS = (MEAN_CPUE_KG * AREA),
                     N_STATIONS = length(unique(GIS_STATION)))%>%
      distinct() %>%
      #Sum across strata
      group_by(AKFIN_SURVEY_YEAR, MAT_SEX) %>%
      dplyr::reframe(AREA=sum(AREA),
                     MEAN_CPUE = sum(MEAN_CPUE),
                     SD_CPUE = sqrt(sum(VAR_CPUE)),
                     N_CPUE = sum(N_CPUE),
                     MEAN_CPUE_KG = sum(MEAN_CPUE_KG),
                     SD_CPUE_KG = sqrt(sum(VAR_CPUE_KG)),
                     N_CPUE_KG = sum(N_CPUE_KG),
                     ABUNDANCE = sum(ABUNDANCE),
                     ABUNDANCE_CI = 1.96*(SD_CPUE),
                     BIOMASS = sum(BIOMASS),
                     BIOMASS_CI = 1.96*(SD_CPUE_KG),
                     N_STATIONS = sum(N_STATIONS)) %>%
      dplyr::mutate(N_STATIONS = ifelse((AKFIN_SURVEY_YEAR == 2000
                                         & stock == "BBRKC"), 135, N_STATIONS)) %>%
      ungroup() -> bio_abund_df
    
    
    return(list(bio_abund_df))
  } # generates biomass and abundance timeseries
  
  table_plot_fun<-function(bio_abund_out_EBS, bio_abund_out_NBS, stock){
    
    
    #Create lookup table for male cutlines 
    mat_lookup <- tibble(stock = c("BBRKC", "PribRKC", "PribBKC", "StMattBKC", "TannerE", "TannerW", "Snow", "Hybrid", "Hair","NorthRKC","NSRKC", "BKCNBS", "SnowNBS"),
                         cutline = c(120, 120, 120, 105, 113, 103, 95, 95, NA, 120,94, 105, 68),
                         legal = c(135, 135, 135, 120, 120, 110, 78, 78, 83, 135,104, 120, 78),
                         recruit = c(134, 134, 134, 119, 124, 124, 101, 101, NA, 134,NA,NA,NA),
                         preferred = c(NA, NA, NA, NA, 125, 125, 102, NA, NA, NA, NA, NA, 102),
                         title =c("Bristol Bay Red King Crab", "Pribilof Islands Red King Crab", "Pribilof Islands Blue King Crab",
                                  "St. Matthew Island Blue King Crab", "Tanner Crab East", "Tanner Crab West", "Snow Crab", 
                                  "Hybrid Tanner-Snow Crab","Hair Crab","Northern District Red King Crab", "Norton Sound Red King Crab","Northern Bering Sea Blue King Crab","Northern Bering Sea Snow Crab"))
 
  
    # Generate abundance plots
    if(is.null(bio_abund_out_EBS) == FALSE){
      bio_abund_out_EBS %>%
        dplyr::select(!N_STATIONS) -> bio_abund_out_FULL
    } else{
      bio_abund_out_FULL <- NULL
    } 
    
    
    bio_abund_out_NBS -> NBS
    
    
    
    #create data frame for 2020 missing data so ggplot plots it
    data.frame(AKFIN_SURVEY_YEAR = rep(2020, length(levels(factor(bio_abund_out_FULL$MAT_SEX)))), 
                MAT_SEX = c(levels(factor(bio_abund_out_FULL$MAT_SEX))),
                MEAN_CPUE = NA, SD_CPUE = NA, N_CPUE = NA, MEAN_CPUE_KG = NA, AREA = NA, SD_CPUE_KG = NA, 
                N_CPUE_KG = NA, ABUNDANCE = NA, ABUNDANCE_CI = NA, BIOMASS = NA, BIOMASS_CI = NA) -> missing
      
    #Bind missing and dummy data frames to abundance/biomass data, filter for maturity-sex categories for plotting
    rbind(missing, bio_abund_out_FULL) %>% #Replace bio_abudnd_out_FULL with incomplete, to only use years with full TS, can add dummy to use 2022 dummy data
        mutate(REGION = "EBS") -> plot_dat  
    
    #specify plot labels
    labs <- c(paste0("Mature male (\u2265", mat_lookup$cutline[mat_lookup$stock == stock], " mm)"),
              "Mature female")
    
    
    names(labs) <- c("Mature Male", "Mature Female") 
    
    # Specify plot title
    title <- mat_lookup$title[mat_lookup$stock == stock]
      
    # Specify plot layers
    if (is.null(NBS) == FALSE & !stock %in% c("BKCNBS", "NSRKC")){
      
      nbs_col <- "#009966"
      
      nbs_points <- geom_point(NBS %>%
                                 filter(MAT_SEX %in% names(labs))%>%
                                 mutate(MAT_SEX = factor(MAT_SEX, levels = names(labs))),
                               mapping = aes(AKFIN_SURVEY_YEAR, ABUNDANCE/1e6, color = "#009E73"), linewidth = 2)
      
      nbs_error <- geom_errorbar(NBS %>%
                                   filter(MAT_SEX %in% names(labs)) %>%
                                   mutate(MAT_SEX = factor(MAT_SEX, levels = names(labs))),
                                 mapping = aes(x = AKFIN_SURVEY_YEAR, ymin = ifelse(((ABUNDANCE-ABUNDANCE_CI)/1e6)<0, 0, 
                                                                                    ((ABUNDANCE-ABUNDANCE_CI)/1e6)), ymax = (ABUNDANCE+ABUNDANCE_CI)/1e6), color = nbs_col, linewidth = 0.75)
      
      ebs_nbs_labs <- scale_color_manual(values = c("#006699", "#009E73"), labels = c("Eastern Bering Sea", 
                                                                                      "Northern Bering Sea")) 
      ebs_line <- geom_line(plot_dat %>%
                              filter(MAT_SEX %in% names(labs)) %>%
                              mutate(MAT_SEX = factor(MAT_SEX, levels = names(labs))),
                            mapping = aes(AKFIN_SURVEY_YEAR, ABUNDANCE/1e6,
                                          color = "#006699"), linewidth = 1) 
      
      ebs_ribbon <- geom_ribbon(plot_dat %>%
                                  filter(MAT_SEX %in% names(labs)) %>%
                                  mutate(MAT_SEX = factor(MAT_SEX, levels = names(labs))),
                                mapping = aes(x = AKFIN_SURVEY_YEAR, 
                                              ymin = ifelse(((ABUNDANCE-ABUNDANCE_CI)/1e6)<0, 0,
                                                            ((ABUNDANCE-ABUNDANCE_CI)/1e6)), 
                                              ymax = (ABUNDANCE+ABUNDANCE_CI)/1e6), 
                                alpha = 0.2, fill = "#006699")
      
      
    } else if (is.null(NBS) == FALSE & stock %in% c("BKCNBS", "NSRKC")){
      
      plot_dat <- NULL
      
      nbs_points <- geom_point(NBS %>%
                                 filter(MAT_SEX %in% names(labs))%>%
                                 mutate(MAT_SEX = factor(MAT_SEX, levels = names(labs))),
                               mapping = aes(AKFIN_SURVEY_YEAR, ABUNDANCE/1e6),color = "#009E73", linewidth = 2)
      
      nbs_error <-  geom_errorbar(NBS %>%
                                    filter(MAT_SEX %in% names(labs)) %>%
                                    mutate(MAT_SEX = factor(MAT_SEX, levels = names(labs))) -> tt,
                                  mapping = aes(x = AKFIN_SURVEY_YEAR, ymin = ((ABUNDANCE-ABUNDANCE_CI)/1e6),
                                                ymax = ((ABUNDANCE+ABUNDANCE_CI)/1e6)), 
                                  color = "#009E73", linewidth = 0.75)
      
      
      
      ebs_nbs_labs <- NULL
      ebs_line <- NULL
      ebs_ribbon <- NULL
      plot_dat <- NBS
      
    } else{
      nbs_col <- NA
      
      nbs_points <- NULL
      nbs_error <- NULL
      
      ebs_nbs_labs <- NULL
      
      ebs_line <- geom_line(plot_dat %>%
                              filter(MAT_SEX %in% names(labs)) %>%
                              mutate(MAT_SEX = factor(MAT_SEX, levels = names(labs))),
                            mapping = aes(AKFIN_SURVEY_YEAR, ABUNDANCE/1e6),
                            color = "#006699", linewidth = 1) 
      
      ebs_ribbon <- geom_ribbon(plot_dat %>%
                                  filter(MAT_SEX %in% names(labs)) %>%
                                  mutate(MAT_SEX = factor(MAT_SEX, levels = names(labs))),
                                mapping = aes(x = AKFIN_SURVEY_YEAR, ymin = ifelse(((ABUNDANCE-ABUNDANCE_CI)/1e6)<0, 0,
                                                                                   ((ABUNDANCE-ABUNDANCE_CI)/1e6)), 
                                              ymax = (ABUNDANCE+ABUNDANCE_CI)/1e6), 
                                alpha = 0.2, fill = "#006699")
      
    }
    
    #Extra step to change how geomribbon handles NAs
    GeomRibbon$handle_na <- function(data, params) {  data }
    
    #Specify color blind color palette
    cb <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") #color blind colors
    
    #Plot abundance timeseries save figure
    ggplot()+
      ebs_line+
      ebs_ribbon+
      nbs_points +
      nbs_error +
      ebs_nbs_labs+
      facet_wrap(~MAT_SEX, scales = "free_y", ncol = 1, nrow = 2, labeller = labeller(MAT_SEX = labs)) +
      scale_x_continuous(breaks= seq(min(plot_dat$AKFIN_SURVEY_YEAR), max(plot_dat$AKFIN_SURVEY_YEAR), by=5),minor_breaks = NULL)+
      labs(title=title, y = "Abundance (millions)") +
      theme_bw()+
      theme(legend.position = ifelse(is.null(NBS) == FALSE, "bottom", "none"), legend.title = element_blank(),
            axis.title.x = element_blank()) -> abundance_plot
    
  
    return(list(abundance_plot = abundance_plot))
  } # timeseries figures
  
  # Run functions, save figure
  c(1979:2019, 2021:2023) %>%
  map_df(~bio_abund(catch_rkc, strata_rkc, "BBRKC", .x)) -> bio_abund_BBRKC
  
  table_plot_fun(bio_abund_BBRKC, NULL, "BBRKC") -> plot_out

  ggsave(plot = plot_out$abundance_plot,
         paste0("./Figures/CPS1 Report/Abundance_timeseries", "BBRKC", ".png"),
         height=7.5, width=10, units="in")

# REQUEST 2: MAP OF TRAWL CLOSURE AREAS W/CPS1 SURVEY AREA --------------------------------------------------------------------------------
  # Load boundaries
  st_read(survey_gdb,layer="Area516") ->  area516
  
  st_read(survey_gdb,layer="BycatchZone1") ->  zone1
  
  st_read(survey_gdb,layer="NBBTCA") ->  nbbtca
  st_read(survey_gdb,layer="TogiakTrawlArea") ->  togtrawl
  st_read("./Data/BBDistrict.gdb", layer = "BB_District") -> bb_dist
  
  # Transform plot boundary
  plot.boundary.untrans <- data.frame(y = c(54, 59.5), 
                                      x = c(-168, -158)) 
  
  plot.boundary <-  plot.boundary.untrans %>%
    sf::st_as_sf(coords = c(x = "x", y = "y"), crs = sf::st_crs(4326)) %>%
    sf::st_transform(crs = map.crs) %>%
    sf::st_coordinates() %>%
    as.data.frame() %>%
    dplyr::rename(x = X, y = Y) # plot boundary projected
  
  breaks.x <- map_layers$lon.breaks[(map_layers$lon.breaks >= plot.boundary.untrans$x[1] &  # set lon breaks
                                       map_layers$lon.breaks < plot.boundary.untrans$x[2]) == TRUE]
  
  breaks.y <- map_layers$lat.breaks[(map_layers$lat.breaks > plot.boundary.untrans$y[1] & # set lat breaks
                                       map_layers$lat.breaks < plot.boundary.untrans$y[2]) == TRUE]
 # Plot
 ggplot() +
   geom_sf(data = st_transform(map_layers$bathymetry, map.crs), color=alpha("grey70")) +
   geom_sf(data = st_as_sf(RKCSA_sub), mapping = aes(fill = "firebrick3"), color = "black", alpha= 0.5, linewidth = 1) +
   geom_sf(data = st_as_sf(RKCSA), aes(fill = "firebrick2"),  color = "black", alpha =0.5, linewidth = 0.5) +
   geom_sf(data = st_as_sf(bb_dist), fill = NA, aes(colour = "black"), linewidth = 1) +
   #geom_sf(data = st_as_sf(zone1$Shape), fill = "blue", alpha = 0.15, aes(color = "black"), linewidth = 0.5) +
   geom_sf_pattern(data = st_as_sf(zone1$Shape),
                   aes(pattern_type = "stripe", pattern_angle = 30), fill = NA, color = "black", pattern_alpha = 0.15)+
   geom_sf(data = st_as_sf(area516$Shape), fill = NA, aes(color = "darkblue"), linewidth = 1) +
   geom_sf(data = st_as_sf(togtrawl$Shape), aes(fill = "gold2"), alpha = 0.5) +
   geom_sf_pattern(data = st_as_sf(nbbtca$Shape),
                   aes(pattern_type = "stripe", pattern_angle = 120), fill = NA, color = "black", pattern_alpha = 0.15)+
   geom_sf(data = st_as_sf(CPS1_bound), fill = NA, aes(color = "olivedrab2"), linewidth = 1)+
   geom_sf(data = st_transform(map_layers$akland, map.crs), fill = "grey80") +
   scale_x_continuous(breaks = c(-165, -160), labels = paste0(c(165, 160), "°W"))+
   scale_y_continuous(breaks = c(56, 58), labels = paste0(c(56, 58), "°N"))+
   labs(title = "2023 BBRKC Collaborative Pot Sampling", subtitle = "Survey extent and closure areas")+
   scale_color_manual(values = c("black", "darkblue", "olivedrab2"), 
                      labels = c("Bristol Bay management boundary", "Area 516", "CPS1 survey extent"),
                      name = "") +
   scale_fill_manual(values = c(alpha("firebrick3", 0.5), alpha("indianred1", 0.25), alpha("gold2", 0.5)),
                     labels = c("Red King Crab Savings Area", "Red King Crab Savings Subarea", "Nearshore Bristol Bay Trawl Area"),
                     name = "")+

   scale_pattern_type_manual(values = c("stripe", "stripe"),
                        labels = c("Bycatch Limitation Zone 1", "Nearshore Bristol Bay Trawl Closure Area"),
                        name = "")+
   scale_pattern_angle_continuous(range = c(30, 120),
                                  breaks = c(30, 120),
                              labels = c("Bycatch Limitation Zone 1", "Nearshore Bristol Bay Trawl Closure Area"),
                              name = "")+
   coord_sf(xlim = plot.boundary$x,
            ylim = plot.boundary$y) +
   geom_sf_text(sf::st_as_sf(data.frame(lab= c("50m", "100m"), 
                                        x = c(-165.8, -166.2), y = c(58.3, 56.5)),
                             coords = c(x = "x", y = "y"), crs = sf::st_crs(4326)) %>%
                  sf::st_transform(crs = map.crs),
                mapping = aes(label = lab))+
   guides(color = guide_legend(nrow = 3), fill = guide_legend(nrow = 3), pattern_angle = guide_legend(nrow = 2), 
          pattern_type = "none") +
   theme_bw() +
   theme(axis.title = element_blank(),
         axis.text = element_text(size = 10),
         legend.text = element_text(size = 10),
         legend.title = element_text(size = 10),
         legend.position = "bottom",
         legend.direction = "horizontal",
         plot.title = element_text(face = "bold", size = 15),
         plot.subtitle = element_text(size = 12),
         panel.grid.major = element_blank()) -> CPS1extent_closures
 
 ggsave(plot = CPS1extent_closures, "./Figures/CPS1 Report/CPS1extent_closures.png", height=7, width=10, units="in")
 
 
# REQUEST 3: MAP OF PLANNED VS. SURVEYED CPS1 STATIONS ----------------------------------------------------------------------
 # Load planned pot IDS and coordinates csv, transform lat and lon to mapping crs
 surv_effort <- read.csv("./Data/Pot IDs and Coordinates.csv") %>%
   dplyr::filter(!Pot.ID %in% unique((pot_cpue$POT_ID))) %>%
   rename(POT_ID = Pot.ID, LAT_DD = Latitude, LON_DD = Longitude) %>%
   mutate(VESSEL = "Not surveyed", LON_DD = LON_DD*-1) %>%
   rbind(pot_cpue %>% dplyr::select(c(VESSEL, POT_ID, LAT_DD, LON_DD))) %>%
   sf::st_as_sf(coords = c(x = "LON_DD", y = "LAT_DD"), crs = sf::st_crs(4326)) %>%
   sf::st_transform(crs = map.crs) %>%
   distinct()

 # Transform plot boundary
 plot.boundary.untrans <- data.frame(y = c(54.5, 58.5), 
                                     x = c(-164.8, -159))
 
 plot.boundary <-  plot.boundary.untrans %>%
   sf::st_as_sf(coords = c(x = "x", y = "y"), crs = sf::st_crs(4326)) %>%
   sf::st_transform(crs = map.crs) %>%
   sf::st_coordinates() %>%
   as.data.frame() %>%
   dplyr::rename(x = X, y = Y) # plot boundary projected
 
 breaks.x <- map_layers$lon.breaks[(map_layers$lon.breaks >= plot.boundary.untrans$x[1] &  # set lon breaks
                                      map_layers$lon.breaks < plot.boundary.untrans$x[2]) == TRUE]
 
 breaks.y <- map_layers$lat.breaks[(map_layers$lat.breaks > plot.boundary.untrans$y[1] & # set lat breaks
                                      map_layers$lat.breaks < plot.boundary.untrans$y[2]) == TRUE]
 
 # Make row and column labels
 sta_labs <- data.frame(lab = c("A", "B", "C", "D", "E", "F", "G", "H", "I",
                                "J", "K", "10", "20", "30", "40", "50", "60", "70", "80"),
                        x = c(-164.2729, -163.8046, -163.3358, -162.8671, -162.3984, -161.9298, -161.4614,-160.9909, -160.5225,
                              -160.0544, -159.583, -164.6662, -164.6588, -164.6516, -164.6446, -164.6376, -164.6308, -164.6241,
                              -164.6175),
                        y = c(57.99736, 58.01606, 58.01643, 58.01504, 58.0119, 58.007, 58.00035, 58.02509, 58.01492,
                              58.003, 58.02244, 57.49888, 57.16714, 56.83529, 56.50332, 56.17121, 55.83898, 55.50659,
                              55.17405)) %>%
             sf::st_as_sf(coords = c(x = "x", y = "y"), crs = sf::st_crs(4326)) %>%
             sf::st_transform(crs = map.crs)
           
 
 # Plot
 ggplot() +
    #geom_sf(data = st_as_sf(CPS1_bound), fill = NA, aes(color = "olivedrab2"), linewidth = 1)+
   geom_sf(data = st_transform(map_layers$akland, map.crs), fill = "grey80") +
   geom_sf(data = surv_effort,
           mapping = aes(fill = VESSEL), shape = 21, size = 3)+   
   scale_x_continuous(breaks = c(-165, -160), labels = paste0(c(165, 160), "°W"))+
   scale_y_continuous(breaks = c(56, 58), labels = paste0(c(56, 58), "°N"))+
   labs(title = "2023 BBRKC Collaborative Pot Sampling", subtitle = "Survey vessel effort")+
   coord_sf(xlim = plot.boundary$x,
            ylim = plot.boundary$y) +
   scale_fill_manual(values = c("white", "cyan4", "darkgoldenrod2"),
                     labels = c("Not surveyed", "Silver Spray", "Summer Bay"),
                      name = "")+
   geom_sf_text(data = sta_labs,
                mapping = aes(label = lab))+
   theme_bw() +
   theme(axis.title = element_blank(),
         axis.text = element_text(size = 10),
         legend.text = element_text(size = 10),
         legend.title = element_text(size = 10),
         legend.position = "bottom",
         legend.direction = "horizontal",
         plot.title = element_text(face = "bold", size = 15),
         plot.subtitle = element_text(size = 12),
         panel.grid.major = element_blank()) -> survey_effort
 
 ggsave(plot = survey_effort, "./Figures/CPS1 Report/CPS1surveyeffort.png", height=7, width=10, units="in")

 
# REQUEST 4: SOAK TIME FREQUENCY DISTRIBUTION -------------------------------------------------------------------------------
 ggplot()+
   geom_histogram(potlifts, mapping = aes(x = SOAK_TIME), fill = "cyan4", bins = 35, color = "black")+
   theme_bw() +
   ylab("Frequency")+
   xlab("Soak time (hours)")+
   scale_x_continuous(breaks = seq(0,max(potlifts$SOAK_TIME), 2)) +
   labs(title = "2023 BBRKC Collaborative Pot Sampling", subtitle = "Soak time")+
   theme(axis.title = element_text(size = 12),
         axis.text = element_text(size = 12),
         legend.text = element_text(size = 10),
         legend.title = element_text(size = 10),
         legend.position = "bottom",
         legend.direction = "horizontal",
         plot.title = element_text(face = "bold", size = 15),
         plot.subtitle = element_text(size = 12),
         panel.grid.major = element_blank()) -> soaktime_freq
 
 ggsave(plot = soaktime_freq, "./Figures/CPS1 Report/CPS1soaktime_freq.png", height=7, width=10, units="in")

# REQUEST 5: TEMPERATURE MAP (2 EBS MAPS FROM NMFS SURVEY + CPS1 MAP) --------------------------------------------------------

 # Load function to generate NMFS temperature maps for Bristol Bay
 temp_map_ebs_nbs <- function(haul_ebs, years){
  
   haul_ebs %>%
       dplyr::filter(SURVEY_YEAR %in% years, HAUL_TYPE != 17) %>%
       dplyr::select("REGION", "SURVEY_YEAR", "GIS_STATION", 
                     "MID_LATITUDE", "MID_LONGITUDE", "GEAR_TEMPERATURE") %>%
       dplyr::filter(!is.na(GEAR_TEMPERATURE), 
                     !is.na(MID_LATITUDE), 
                     !is.na(MID_LONGITUDE)) -> temp
   
   # load EBS-NBS survey extent for masking
   interpolation.crs <- map.crs
  
   # Make raster for interpolation
   cell.resolution = 1000
   in.crs = "+proj=longlat +datum=NAD83"
   extrap.box = c(xmn = -165, xmx = -157, ymn = 50, ymx = 59)
                                                                 
   # Transform plot boundary
   plot.boundary.untrans <- data.frame(y = c(50, 59), 
                                       x = c(-165, -157))
   
   plot.boundary <-  plot.boundary.untrans %>%
     sf::st_as_sf(coords = c(x = "x", y = "y"), crs = sf::st_crs(4326)) %>%
     sf::st_transform(crs = map.crs) %>%
     sf::st_coordinates() %>%
     as.data.frame() %>%
     dplyr::rename(x = X, y = Y) # plot boundary projected
   
   breaks.x <- map_layers$lon.breaks[(map_layers$lon.breaks >= plot.boundary.untrans$x[1] &  # set lon breaks
                                        map_layers$lon.breaks < plot.boundary.untrans$x[2]) == TRUE]
   
   breaks.y <- map_layers$lat.breaks[(map_layers$lat.breaks > plot.boundary.untrans$y[1] & # set lat breaks
                                        map_layers$lat.breaks < plot.boundary.untrans$y[2]) == TRUE]
   
   
   n_dim <- floor(abs(plot.boundary$x[1] - plot.boundary$x[2]))/cell.resolution
   
   sp_interp.raster <- raster::raster(xmn = plot.boundary$x[1], 
                                      xmx = plot.boundary$x[2], 
                                      ymn = plot.boundary$y[1], 
                                      ymx = plot.boundary$y[2], 
                                      nrow = n_dim, 
                                      ncol = n_dim)
   
   raster::projection(sp_interp.raster) <- interpolation.crs
   
   # Transform data for interpolation ----
   sp_interp.df <- unique(temp)
   sp::coordinates(sp_interp.df) <- c(x = "MID_LONGITUDE", y = "MID_LATITUDE")
   sp::proj4string(sp_interp.df) <- sp::CRS(in.crs)
   sp_interp.df <- sp::spTransform(sp_interp.df, sp::CRS(interpolation.crs))
   
   # Set up a new IDW for ordinary kriging ----
   idw_vgm_fit <- gstat::gstat(formula = GEAR_TEMPERATURE ~ 1, 
                               locations = sp_interp.df, 
                               nmax = Inf)
   
   # Ordinary Kriging: Stein's Matern VGM----
   ste.vgfit <- gstat::fit.variogram(variogram(idw_vgm_fit), 
                                     vgm(c("Ste")))
   
   ste_fit <- gstat::gstat(formula = GEAR_TEMPERATURE ~ 1, 
                           locations = sp_interp.df, 
                           model = ste.vgfit, 
                           nmax = Inf)
   
   ste.predict <- predict(ste_fit, as(sp_interp.raster, "SpatialGrid"))
   
    # write unmasked surfaces to raster, stacked by year
   ste.predict %>%
     raster::raster(.) %>%
     mask(st_transform(map_layers$survey.area, map.crs)) ->  temp_rast
   
   # extract interpolated data from raster to data frame
   coords<-coordinates(temp_rast)
   
   temp_df_nmfs <-na.omit(data.frame(coords, temperature = temp_rast@data@values, year = years))
   
   
   temp_breaks <- c(-Inf, seq(-1,8,1), Inf)
   viridis_option <- "H" # viridis turbo palette
   n_temp_breaks <- length(temp_breaks)-1
   
   # Transform plot boundary
   plot.boundary.untrans <- data.frame(y = c(54.5, 58.5), 
                                       x = c(-164.8, -159))
   
   plot.boundary <-  plot.boundary.untrans %>%
     sf::st_as_sf(coords = c(x = "x", y = "y"), crs = sf::st_crs(4326)) %>%
     sf::st_transform(crs = map.crs) %>%
     sf::st_coordinates() %>%
     as.data.frame() %>%
     dplyr::rename(x = X, y = Y) # plot boundary projected
   
   breaks.x <- map_layers$lon.breaks[(map_layers$lon.breaks >= plot.boundary.untrans$x[1] &  # set lon breaks
                                        map_layers$lon.breaks < plot.boundary.untrans$x[2]) == TRUE]
   
   breaks.y <- map_layers$lat.breaks[(map_layers$lat.breaks > plot.boundary.untrans$y[1] & # set lat breaks
                                        map_layers$lat.breaks < plot.boundary.untrans$y[2]) == TRUE]
   # Year labels
   sf::st_as_sf(data.frame(lab= paste("", years, "\nNMFS"), 
                           x = c(-160), y = c(55.2)),
                coords = c(x = "x", y = "y"), crs = sf::st_crs(4326)) %>%
     sf::st_transform(crs = map.crs) %>%
     cbind(years, st_coordinates(.)) -> year_lab
   
   # Plot interpolated data
   ggplot2::ggplot() +
     ggplot2::geom_tile(data = temp_df_nmfs, 
                        aes(x = x, 
                            y = y,
                            fill = cut(temperature, 
                                       breaks = temp_breaks))) +
     geom_sf(data = st_transform(map_layers$survey.area, map.crs), fill = NA, linewidth = 0.5) +
     geom_sf(data = st_as_sf(CPS1_bound), fill = NA, color = "black", linewidth = 1)+
     geom_sf(data = st_transform(map_layers$akland, map.crs), fill = "grey80") +
     scale_x_continuous(breaks = c(-165, -160), labels = paste0(c(165, 160), "°W"))+
     scale_y_continuous(breaks = c(56, 58), labels = paste0(c(56, 58), "°N"))+
     coord_sf(xlim = plot.boundary$x,
              ylim = plot.boundary$y) +
     ggplot2::scale_fill_manual(name = "Temperature (°C)", values = viridis_pal(option = viridis_option)(n_temp_breaks),
                                labels = c(expression(""<=-1), "-0.9-0", "0.1-1", "1.1-2", "2.1-3", "3.1-4",
                                           "4.1-5", "5.1-6", "6.1-7", "7.1-8", ">8.1"), drop = FALSE) +
     geom_shadowtext(year_lab,
                  mapping = aes(label = lab, x = X, y = Y), size = 8,  color = "black", bg.color = "white")+
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.position = "bottom",
        legend.direction = "horizontal",
        plot.title = element_text(face = "bold", size = 15),
        plot.subtitle = element_text(size = 12))  -> temp_map
   return(temp_map)
   
 }
 
 2022:2023 %>%
   map(~temp_map_ebs_nbs(all_hauls, .x)) -> nmfs_temp_out
 
 ggsave(plot = nmfs_temp_out[[1]], "./Figures/CPS1 Report/NMFStemp_22.png", height=7, width=10, units="in")
 ggsave(plot = nmfs_temp_out[[2]], "./Figures/CPS1 Report/NMFStemp_23.png", height=7, width=10, units="in")
 
 
 # Generate temperature map for CPS1
 # Year labels
 sf::st_as_sf(data.frame(lab= paste("", 2023, "\nCPS1"), 
                         x = c(-160), y = c(55.2)),
              coords = c(x = "x", y = "y"), crs = sf::st_crs(4326)) %>%
   sf::st_transform(crs = map.crs) %>%
   cbind(2023, st_coordinates(.)) -> CPS1year_lab
 
 # Temp breaks
 temp_breaks <- c(-Inf, seq(-1,8,1), Inf)
 viridis_option <- "H" # viridis turbo palette
 n_temp_breaks <- length(temp_breaks)-1
 
 # Plot
 ggplot() +
   ggplot2::geom_tile(data = temp_df, 
                      aes(x = x, 
                          y = y,
                          fill = cut(temperature, 
                                     breaks = temp_breaks))) +
   geom_sf(data = st_as_sf(CPS1_bound), fill = NA, color = "black", linewidth = 1)+
   geom_sf(data = CPS1_bathy, color=alpha("white")) +
   #geom_sf(data = st_as_sf(BB_strata), fill = NA, mapping = aes(color = "black"), linewidth = 1) +
   geom_sf(data = temploggers, color = "black", size = 1)+
   geom_sf(data = st_transform(map_layers$akland, map.crs), fill = "grey80") +
   scale_x_continuous(breaks = c(-165, -160), labels = paste0(c(165, 160), "°W"))+
   scale_y_continuous(breaks = c(56, 58), labels = paste0(c(56, 58), "°N"))+
   coord_sf(xlim = plot.boundary$x,
            ylim = plot.boundary$y) +
   ggplot2::scale_fill_manual(name = "Temperature (°C)", values = viridis_pal(option = viridis_option)(n_temp_breaks),
    labels = c(expression(""<=-1), "-0.9-0", "0.1-1", "1.1-2", "2.1-3", "3.1-4",
    "4.1-5", "5.1-6", "6.1-7", "7.1-8", ">8.1"), drop = FALSE) +
   guides(color = guide_legend(nrow = 2), fill = guide_legend(nrow = 2)) +
   geom_shadowtext(CPS1year_lab,
                   mapping = aes(label = lab, x = X, y = Y), size = 8,  color = "black", bg.color = "white")+
   geom_shadowtext(data = (sf::st_as_sf(data.frame(lab= c("35m", "45m", "55m", "65m", "75m", "85m"), 
                                                   x = c(-160.6, -161.6, -161.6, -161.2, -163.5, -163.6), 
                                                   y = c(56.1, 56, 57.5, 57.2, 56.76, 56.38)),
                                        coords = c(x = "x", y = "y"), crs = sf::st_crs(4326)) %>%
                             sf::st_transform(crs = map.crs) %>%
                             cbind(st_coordinates(.))),
                   mapping = aes(label = lab, x = X, y = Y), color = "black", bg.color = "white", size = 3.5)+
   theme_bw() +
   theme(axis.title = element_blank(),
         axis.text = element_text(size = 15),
         legend.key.width = unit(12, "mm"),
         legend.position = "bottom",
         legend.direction = "horizontal",
         plot.title = element_text(face = "bold", size = 15),
         plot.subtitle = element_text(size = 12)) -> CPS1_tempdiscrete
 
 ggsave(plot = CPS1_tempdiscrete, "./Figures/CPS1 Report/ CPS1_tempdiscrete.png", height=7, width=10, units="in")

# REQUEST 6: TEMPERATURE MAPS WITH CRAB (2 EBS MAPS FROM NMFS SURVEY + CPS1 MAP) --------------------------------------------------
 
 # Load and run function to calculate BBRKC cpue for NMFS summer survey
 cpue_calc <- function(data_crab_EBS, data_crab_NBS, calc_factor, year){
   
   if (is.null(data_crab_EBS) == "FALSE" & is.null(data_crab_NBS) == "TRUE"){
     data_crab_EBS$REGION = "EBS"
     data_crab <- data_crab_EBS
     
   } else if (is.null(data_crab_EBS) == "TRUE" & is.null(data_crab_NBS) == "FALSE"){
     data_crab_NBS$REGION = "NBS"
     data_crab <- data_crab_NBS
     
   } else {
     data_crab_EBS$REGION = "EBS"
     data_crab_NBS$REGION = "NBS"
     data_crab <- rbind(data_crab_EBS[,-c(28,40)], data_crab_NBS)#for hybrid crab the EBS files don't have extra columns, so 28 and 40 &41 must be removed
   }
   
   species <- unique(data_crab$SPECIES_NAME)
   
   #Management area lookup table (file output by stock_stations.R script generated directly from strata tables)
   sta <- read.csv("./Data/stock_stations_updated.csv", na.strings = "") %>% as_tibble()
   stations <- na.omit(pull(sta, calc_factor))
   
   #Create lookup table for male cutlines 
   mat_lookup <- tibble(stock = c("BBRKC", "PribRKC", "PribBKC", "StMattBKC", "TannerE", "TannerW", "Snow", "Hybrid", "Hair", "NSRKC", "BKCNBS", "SnowNBS"),
                        cutline = c(120, 120, 120, 105, 113, 103, 95, 95, NA, 94, 105, 68),
                        legal = c(135, 135, 135, 120, 120, 110, 78, 78, 83, 104, 120, 78)) 
   
   #Add STOCK column to data_crab based on stock stations
   data_crab %>%
     mutate(STOCK = case_when((SPECIES_NAME == "Red King Crab" & GIS_STATION %in% sta$BBRKC)~ "BBRKC",
                              (SPECIES_NAME == "Red King Crab" & GIS_STATION %in% sta$PribRKC) ~ "PribRKC",
                              (SPECIES_NAME == "Red King Crab" & GIS_STATION %in% sta$NSRKC) ~ "NSRKC",
                              (SPECIES_NAME == "Blue King Crab" & GIS_STATION %in% sta$PribBKC) ~ "PribBKC",
                              (SPECIES_NAME == "Blue King Crab" & GIS_STATION %in% sta$StMattBKC) ~"StMattBKC",
                              (SPECIES_NAME == "Blue King Crab" & GIS_STATION %in% sta$BKCNBS) ~ "BKCNBS",
                              (SPECIES_NAME == "Bairdi Tanner Crab" & GIS_STATION %in% sta$TannerE) ~ "TannerE",
                              (SPECIES_NAME == "Bairdi Tanner Crab" & GIS_STATION %in% sta$TannerW) ~ "TannerW",
                              (SPECIES_NAME == "Opilio Crab" & GIS_STATION %in% sta$EBS_NBS) ~ "Snow",
                              (SPECIES_NAME == "Chionoecetes Hybrid Crab" & GIS_STATION %in% sta$EBS_NBS) ~ "Hybrid",
                              (SPECIES_NAME == "Erimacrus Hair Crab" & GIS_STATION %in% sta$EBS_NBS) ~ "Hair",
                              TRUE ~ "Unstratified")) -> data_crab2
   
   #Create maturity levels
   if (species != "Erimacrus Hair Crab"){
     data_crab2 %>%
       filter(SEX %in% 1:2,
              AKFIN_SURVEY_YEAR %in% year) %>%
       #Red King Crab male maturity
       mutate(MAT_SEX = case_when((SEX == 1 & SPECIES_NAME == "Red King Crab" & REGION == "EBS" & LENGTH_1MM >= mat_lookup$cutline[mat_lookup$stock == "BBRKC"]) ~ "Mature Male",
                                  (SEX == 1 & SPECIES_NAME == "Red King Crab" & REGION == "EBS" & LENGTH_1MM < mat_lookup$cutline[mat_lookup$stock == "BBRKC"]) ~ "Immature Male",
                                  (SEX == 1 & SPECIES_NAME == "Red King Crab" & REGION == "NBS" & LENGTH_1MM >= mat_lookup$cutline[mat_lookup$stock == "NSRKC"]) ~ "Mature Male",
                                  (SEX == 1 & SPECIES_NAME == "Red King Crab" & REGION == "NBS" & LENGTH_1MM < mat_lookup$cutline[mat_lookup$stock == "NSRKC"]) ~ "Immature Male",
                                  #Blue King Crab male maturity
                                  (SEX == 1 & SPECIES_NAME == "Blue King Crab" & REGION == "EBS" & MID_LATITUDE <= 58.7 & LENGTH_1MM >= mat_lookup$cutline[mat_lookup$stock == "PribBKC"]) ~ "Mature Male",
                                  (SEX == 1 & SPECIES_NAME == "Blue King Crab" & REGION == "EBS" & MID_LATITUDE <= 58.7 & LENGTH_1MM < mat_lookup$cutline[mat_lookup$stock == "PribBKC"]) ~ "Immature Male",
                                  (SEX == 1 & SPECIES_NAME == "Blue King Crab" & REGION == "EBS" & MID_LATITUDE > 58.7 & LENGTH_1MM >= mat_lookup$cutline[mat_lookup$stock == "StMattBKC"]) ~ "Mature Male",
                                  (SEX == 1 & SPECIES_NAME == "Blue King Crab" & REGION == "EBS" & MID_LATITUDE > 58.7 & LENGTH_1MM < mat_lookup$cutline[mat_lookup$stock == "StMattBKC"]) ~ "Immature Male",
                                  (SEX == 1 & SPECIES_NAME == "Blue King Crab" & REGION == "NBS" & LENGTH_1MM >= mat_lookup$cutline[mat_lookup$stock == "BKCNBS"]) ~ "Mature Male",
                                  (SEX == 1 & SPECIES_NAME == "Blue King Crab" & REGION == "NBS" & LENGTH_1MM < mat_lookup$cutline[mat_lookup$stock == "BKCNBS"]) ~ "Immature Male",
                                  #Tanner Crab male maturity
                                  (SEX == 1 & STOCK == "TannerE" & REGION == "EBS" & WIDTH_1MM >= mat_lookup$cutline[mat_lookup$stock == "TannerE"]) ~ "Mature Male",
                                  (SEX == 1 & STOCK == "TannerE" & REGION == "EBS" & WIDTH_1MM < mat_lookup$cutline[mat_lookup$stock == "TannerE"]) ~ "Immature Male",
                                  (SEX == 1 & STOCK == "TannerW" & REGION == "EBS" & WIDTH_1MM >= mat_lookup$cutline[mat_lookup$stock == "TannerW"]) ~ "Mature Male",
                                  (SEX == 1 & STOCK == "TannerW" & REGION == "EBS" & WIDTH_1MM < mat_lookup$cutline[mat_lookup$stock == "TannerW"]) ~ "Immature Male",
                                  (SEX == 1 & SPECIES_NAME == "Bairdi Tanner Crab" & REGION == "NBS" & MID_LONGITUDE > -166 & WIDTH_1MM >= mat_lookup$cutline[mat_lookup$stock == "TannerE"]) ~ "Mature Male",
                                  (SEX == 1 & SPECIES_NAME == "Bairdi Tanner Crab" & REGION == "NBS" & MID_LONGITUDE > -166 & WIDTH_1MM < mat_lookup$cutline[mat_lookup$stock == "TannerE"]) ~ "Immature Male",
                                  (SEX == 1 & SPECIES_NAME == "Bairdi Tanner Crab" & REGION == "NBS" & MID_LONGITUDE < -166 & WIDTH_1MM >= mat_lookup$cutline[mat_lookup$stock == "TannerW"]) ~ "Mature Male",
                                  (SEX == 1 & SPECIES_NAME == "Bairdi Tanner Crab" & REGION == "NBS" & MID_LONGITUDE < -166 & WIDTH_1MM < mat_lookup$cutline[mat_lookup$stock == "TannerW"]) ~ "Immature Male",
                                  #Snow Crab male maturity
                                  (SEX == 1 & STOCK == "Snow" & REGION == "EBS" & WIDTH_1MM >= mat_lookup$cutline[mat_lookup$stock == "Snow"]) ~ "Mature Male",
                                  (SEX == 1 & STOCK == "Snow" & REGION == "EBS" & WIDTH_1MM < mat_lookup$cutline[mat_lookup$stock == "Snow"]) ~ "Immature Male",
                                  (SEX == 1 & STOCK == "Snow" & REGION == "NBS" & WIDTH_1MM >= mat_lookup$cutline[mat_lookup$stock == "SnowNBS"]) ~ "Mature Male",
                                  (SEX == 1 & STOCK == "Snow" & REGION == "NBS" & WIDTH_1MM < mat_lookup$cutline[mat_lookup$stock == "SnowNBS"]) ~ "Immature Male",
                                  #Hybrid Crab male maturity
                                  (SEX == 1 & STOCK == "Hybrid" & WIDTH_1MM >= mat_lookup$cutline[mat_lookup$stock == "Hybrid"]) ~ "Mature Male",
                                  (SEX == 1 & STOCK == "Hybrid" & WIDTH_1MM < mat_lookup$cutline[mat_lookup$stock == "Hybrid"]) ~ "Immature Male",
                                  #All species female maturity
                                  (SEX == 2 & CLUTCH_SIZE >= 1) ~ "Mature Female",
                                  (SEX == 2 & CLUTCH_SIZE == 0) ~ "Immature Female")) -> mature #for mature/immature males and females
     data_crab2 %>%  
       dplyr::filter(SEX %in% 1, 
                     AKFIN_SURVEY_YEAR %in% year) %>%
       #Red King Crab legal
       mutate(MAT_SEX = case_when((SPECIES_NAME == "Red King Crab" & REGION == "EBS" & LENGTH_1MM >= mat_lookup$legal[mat_lookup$stock == "BBRKC"]) ~ "Legal Male",
                                  (SPECIES_NAME == "Red King Crab" & REGION == "NBS" & LENGTH_1MM >= mat_lookup$legal[mat_lookup$stock == "NSRKC"]) ~ "Legal Male",
                                  #Blue King Crab legal
                                  (SPECIES_NAME == "Blue King Crab" & REGION == "EBS" & MID_LATITUDE <= 58.7 & LENGTH_1MM >= mat_lookup$legal[mat_lookup$stock == "PribBKC"]) ~"Legal Male",
                                  (SPECIES_NAME == "Blue King Crab" & REGION == "EBS" & MID_LATITUDE > 58.7 & LENGTH_1MM >= mat_lookup$legal[mat_lookup$stock == "StMattBKC"]) ~"Legal Male",
                                  (SPECIES_NAME == "Blue King Crab" & REGION == "NBS" & LENGTH_1MM >= mat_lookup$legal[mat_lookup$stock == "BKCNBS"]) ~ "Legal Male",
                                  #Tanner Crab legal
                                  (STOCK == "TannerE" & REGION == "EBS" & WIDTH_1MM >= mat_lookup$legal[mat_lookup$stock == "TannerE"]) ~ "Legal Male",
                                  (STOCK == "TannerW" & REGION == "EBS" & WIDTH_1MM >= mat_lookup$legal[mat_lookup$stock == "TannerW"]) ~ "Legal Male",
                                  (SPECIES_NAME == "Bairdi Tanner Crab" & REGION == "NBS" & MID_LONGITUDE < -166 & WIDTH_1MM >= mat_lookup$legal[mat_lookup$stock == "TannerW"]) ~ "Legal Male",
                                  (SPECIES_NAME == "Bairdi Tanner Crab" & REGION == "NBS" & MID_LONGITUDE > -166 & WIDTH_1MM >= mat_lookup$legal[mat_lookup$stock == "TannerE"]) ~ "Legal Male",
                                  #Snow Crab legal
                                  (STOCK == "Snow" & REGION == "EBS" & WIDTH_1MM >= mat_lookup$legal[mat_lookup$stock == "Snow"]) ~ "Legal Male",
                                  (STOCK == "Snow" & REGION == "NBS" & WIDTH_1MM >= mat_lookup$legal[mat_lookup$stock == "SnowNBS"]) ~ "Legal Male",
                                  #Hybrid Crab legal
                                  (STOCK == "Hybrid" & WIDTH_1MM >= mat_lookup$legal[mat_lookup$stock == "Hybrid"]) ~ "Legal Male")) -> leg
     
     data_crab2 %>%  
       dplyr::filter(SEX %in% 1, 
                     AKFIN_SURVEY_YEAR %in% year) %>%
       #Red King Crab legal
       mutate(MAT_SEX = case_when((SPECIES_NAME == "Bairdi Tanner Crab" & WIDTH_1MM >= 125) ~ "Industry Preferred",
                                  #Snow Crab legal and industry preferred
                                  (STOCK == "Snow" & REGION == "EBS" & WIDTH_1MM >= 102) ~ "Industry Preferred",
                                  (STOCK == "Snow" & REGION == "NBS" & WIDTH_1MM >= 102) ~ "Industry Preferred")) -> ind_pref
     
     rbind(mature, leg, ind_pref) -> mat_spec 
     #join "mature" dataframe to legal/prerecruit dataframe to create overlapping size classes
     
   } else{
     data_crab2 %>%
       dplyr::filter(SEX %in% 1:2,
                     AKFIN_SURVEY_YEAR %in% year) %>%
       mutate(MAT_SEX = case_when((SEX == 1 & LENGTH_1MM >= mat_lookup$legal[mat_lookup$stock == "Hair"]) ~ "Legal Male",
                                  (SEX == 1 & LENGTH_1MM < mat_lookup$legal[mat_lookup$stock == "Hair"]) ~ "Sublegal Male",
                                  (SEX == 2) ~ "Female")) -> mat_spec
   }
   
   #Compute recent year per-station cpue by maturity within mgmt area 
   mat_spec %>%
     dplyr::filter(GIS_STATION %in% stations) %>%
     dplyr::mutate(CPUE = SAMPLING_FACTOR/AREA_SWEPT) %>%
     dplyr::group_by(AKFIN_SURVEY_YEAR, GIS_STATION, AREA_SWEPT, HAUL_TYPE, MAT_SEX, MID_LATITUDE, MID_LONGITUDE)%>%
     dplyr::summarise(CPUE = sum(CPUE)) %>%
     na.omit() -> cpue
   
   #Conditionally specifying maturity/sex combos for each stock to expand grid below
   if (species %in% c("Red King Crab", "Blue King Crab", "Chionoecetes Hybrid Crab")) {
     
     mat_sex_combos <- c("Mature Male", "Immature Male", "Mature Female", "Immature Female", "Legal Male") 
     
   } else if (species %in% c("Opilio Crab", "Bairdi Tanner Crab")) {
     
     mat_sex_combos <- c("Mature Male", "Immature Male", "Mature Female", 
                         "Immature Female", "Legal Male", "Industry Preferred")
   } else{
     
     mat_sex_combos <- c("Sublegal Male", "Legal Male", "Female")
     
   }
   
   #Create dataframe with all stations to retain zero catch stations
   expand_grid(MAT_SEX = mat_sex_combos,
               data_crab2 %>%
                 dplyr::filter(GIS_STATION %in% stations,
                               AKFIN_SURVEY_YEAR %in% year) %>% 
                 dplyr::distinct(AKFIN_SURVEY_YEAR, GIS_STATION, 
                                 HAUL_TYPE, AREA_SWEPT, MID_LATITUDE, MID_LONGITUDE)) -> haul_sex_combos
   
   
   #Join to positive catch stations and replace female data with retow data (if applicable) 
   cpue %>%
     right_join(haul_sex_combos,
                by = c("AKFIN_SURVEY_YEAR", "GIS_STATION", "AREA_SWEPT", "HAUL_TYPE", "MAT_SEX", "MID_LATITUDE", "MID_LONGITUDE")) %>%
     replace_na(list(CPUE = 0)) %>%
     group_by(AKFIN_SURVEY_YEAR, GIS_STATION, MAT_SEX) %>%
     nest() %>%
     #Females: replacing original stations with resampled stations in BB retow yrs
     dplyr::mutate(data = purrr::map2(data, MAT_SEX, function(data, mat_sex) {
       if(17 %in% data$HAUL_TYPE & species == "Red King Crab" & (mat_sex == "Mature Female"|mat_sex == "Immature Female")) 
       {data %>% filter(HAUL_TYPE == 17) -> x} else{x <- data %>% filter(HAUL_TYPE != 17)}
       return(x)  
     })) %>%
     unnest() -> station_cpue
   
   
   #Create objects for each maturity/sex category for mapping purposes
   
   #Mature male CPUE
   if (species != "Erimacrus Hair Crab"){
     station_cpue %>%
       filter(MAT_SEX == "Mature Male") ->  mat_male_cpue
   } else{
     mat_male_cpue <- NA
   }
   
   #Legal male CPUE
   station_cpue %>%
     filter(MAT_SEX == "Legal Male") -> legal_male_cpue
   
   #Sublegal male (for Hair Crab)
   if (species == "Erimacrus Hair Crab"){
     station_cpue %>%
       filter(MAT_SEX == "Sublegal Male") -> sublegal_male_cpue
   } else{
     sublegal_male_cpue <- NA
   }
   
   #Industry preferred CPUE (for snow and tanner)
   if (species %in% c("Bairdi Tanner Crab", "Opilio Crab")){
     station_cpue %>%
       filter(MAT_SEX == "Industry Preferred") -> indpref_male_cpue
   } else{
     indpref_male_cpue <- NA
   }
   
   #Immature male CPUE
   if (species != "Erimacrus Hair Crab"){
     station_cpue %>%
       filter(MAT_SEX == "Immature Male") -> imm_male_cpue
   } else{
     imm_male_cpue <- NA
   }
   
   #Mature female CPUE
   if (species != "Erimacrus Hair Crab"){
     station_cpue %>%
       filter(MAT_SEX == "Mature Female") -> mat_fem_cpue
   } else{
     station_cpue %>%
       filter(MAT_SEX == "Female") -> mat_fem_cpue
   }
   
   
   # Immature female CPUE
   if (species != "Erimacrus Hair Crab"){
     station_cpue %>%
       filter(MAT_SEX == "Immature Female") -> imm_fem_cpue 
   } else{
     imm_fem_cpue <- NA
   }
   
   
   
   return(list(mat_male_cpue=mat_male_cpue, legal_male_cpue=legal_male_cpue, sublegal_male_cpue=sublegal_male_cpue, indpref_male_cpue=indpref_male_cpue,
               imm_male_cpue=imm_male_cpue, mat_fem_cpue=mat_fem_cpue, imm_fem_cpue=imm_fem_cpue))
 } 
 
 cpue_calc(catch_rkc, NULL, "BBRKC", 2022:2023) -> NMFScpue_out
 
 # Sum cpue across all mat-sex categories for all crab cpue
 rbind(NMFScpue_out$mat_male_cpue, NMFScpue_out$imm_male_cpue, NMFScpue_out$mat_fem_cpue,NMFScpue_out$imm_fem_cpue) %>%
   group_by(AKFIN_SURVEY_YEAR, GIS_STATION, MID_LATITUDE, MID_LONGITUDE) %>%
   dplyr::reframe(CPUE = sum(CPUE)) -> NMFScpue_sum
 
 nmfs_temp_df <- function(haul_ebs, years){
   
   haul_ebs %>%
     dplyr::filter(SURVEY_YEAR %in% years, HAUL_TYPE != 17) %>%
     dplyr::select("REGION", "SURVEY_YEAR", "GIS_STATION", 
                   "MID_LATITUDE", "MID_LONGITUDE", "GEAR_TEMPERATURE") %>%
     dplyr::filter(!is.na(GEAR_TEMPERATURE), 
                   !is.na(MID_LATITUDE), 
                   !is.na(MID_LONGITUDE)) -> temp
   
   # load EBS-NBS survey extent for masking
   interpolation.crs <- map.crs
   
   # Make raster for interpolation
   cell.resolution = 1000
   in.crs = "+proj=longlat +datum=NAD83"
   extrap.box = c(xmn = -165, xmx = -157, ymn = 50, ymx = 59)
   
   # Transform plot boundary
   plot.boundary.untrans <- data.frame(y = c(50, 59), 
                                       x = c(-165, -157))
   
   plot.boundary <-  plot.boundary.untrans %>%
     sf::st_as_sf(coords = c(x = "x", y = "y"), crs = sf::st_crs(4326)) %>%
     sf::st_transform(crs = map.crs) %>%
     sf::st_coordinates() %>%
     as.data.frame() %>%
     dplyr::rename(x = X, y = Y) # plot boundary projected
   
   breaks.x <- map_layers$lon.breaks[(map_layers$lon.breaks >= plot.boundary.untrans$x[1] &  # set lon breaks
                                        map_layers$lon.breaks < plot.boundary.untrans$x[2]) == TRUE]
   
   breaks.y <- map_layers$lat.breaks[(map_layers$lat.breaks > plot.boundary.untrans$y[1] & # set lat breaks
                                        map_layers$lat.breaks < plot.boundary.untrans$y[2]) == TRUE]
   
   
   n_dim <- floor(abs(plot.boundary$x[1] - plot.boundary$x[2]))/cell.resolution
   
   sp_interp.raster <- raster::raster(xmn = plot.boundary$x[1], 
                                      xmx = plot.boundary$x[2], 
                                      ymn = plot.boundary$y[1], 
                                      ymx = plot.boundary$y[2], 
                                      nrow = n_dim, 
                                      ncol = n_dim)
   
   raster::projection(sp_interp.raster) <- interpolation.crs
   
   # Transform data for interpolation ----
   sp_interp.df <- unique(temp)
   sp::coordinates(sp_interp.df) <- c(x = "MID_LONGITUDE", y = "MID_LATITUDE")
   sp::proj4string(sp_interp.df) <- sp::CRS(in.crs)
   sp_interp.df <- sp::spTransform(sp_interp.df, sp::CRS(interpolation.crs))
   
   # Set up a new IDW for ordinary kriging ----
   idw_vgm_fit <- gstat::gstat(formula = GEAR_TEMPERATURE ~ 1, 
                               locations = sp_interp.df, 
                               nmax = Inf)
   
   # Ordinary Kriging: Stein's Matern VGM----
   ste.vgfit <- gstat::fit.variogram(variogram(idw_vgm_fit), 
                                     vgm(c("Ste")))
   
   ste_fit <- gstat::gstat(formula = GEAR_TEMPERATURE ~ 1, 
                           locations = sp_interp.df, 
                           model = ste.vgfit, 
                           nmax = Inf)
   
   ste.predict <- predict(ste_fit, as(sp_interp.raster, "SpatialGrid"))
   
   # write unmasked surfaces to raster, stacked by year
   ste.predict %>%
     raster::raster(.) %>%
     mask(st_transform(map_layers$survey.area, map.crs)) ->  temp_rast
   
   # extract interpolated data from raster to data frame
   coords<-coordinates(temp_rast)
   
   temp_df_nmfs <-na.omit(data.frame(coords, temperature = temp_rast@data@values, year = years))
   
   return(temp_df_nmfs)
 }
 
 # NMFS TEMP DATA
 nmfs_temp_df(all_hauls, 2022) -> NMFS_temp22
 nmfs_temp_df(all_hauls, 2023) -> NMFS_temp23
 
 # SET UP PLOTTING FEATURES
 
 temp_cpue_map <- function(temp_data, cpue_data, years, mat_sex){
   temp_breaks <- c(-Inf, seq(-1,8,1), Inf)
   viridis_option <- "H" # viridis turbo palette
   n_temp_breaks <- length(temp_breaks)-1
   
   # Transform plot boundary
   plot.boundary.untrans <- data.frame(y = c(54.5, 58.5), 
                                       x = c(-164.8, -159))
   
   plot.boundary <-  plot.boundary.untrans %>%
     sf::st_as_sf(coords = c(x = "x", y = "y"), crs = sf::st_crs(4326)) %>%
     sf::st_transform(crs = map.crs) %>%
     sf::st_coordinates() %>%
     as.data.frame() %>%
     dplyr::rename(x = X, y = Y) # plot boundary projected
   
   breaks.x <- map_layers$lon.breaks[(map_layers$lon.breaks >= plot.boundary.untrans$x[1] &  # set lon breaks
                                        map_layers$lon.breaks < plot.boundary.untrans$x[2]) == TRUE]
   
   breaks.y <- map_layers$lat.breaks[(map_layers$lat.breaks > plot.boundary.untrans$y[1] & # set lat breaks
                                        
                                        map_layers$lat.breaks < plot.boundary.untrans$y[2]) == TRUE]
   
   # Year labels
   sf::st_as_sf(data.frame(lab= paste("", years, "\nNMFS"), 
                           x = c(-160), y = c(55.2)),
                coords = c(x = "x", y = "y"), crs = sf::st_crs(4326)) %>%
     sf::st_transform(crs = map.crs) %>%
     cbind(years, st_coordinates(.)) -> year_lab
   
   if(mat_sex == "All crab"){
     as.data.frame(cpue_data) %>%
       filter(AKFIN_SURVEY_YEAR==years)%>%
       sf::st_as_sf(coords = c(x = "MID_LONGITUDE", y = "MID_LATITUDE"), 
                    crs = sf::st_crs(4326)) %>% 
       sf::st_transform(crs = map.crs) -> cpue_data2
   } else{
     as.data.frame(cpue_data) %>%
       filter(AKFIN_SURVEY_YEAR==years, MAT_SEX == mat_sex)%>%
       sf::st_as_sf(coords = c(x = "MID_LONGITUDE", y = "MID_LATITUDE"), 
                    crs = sf::st_crs(4326)) %>% 
       sf::st_transform(crs = map.crs) -> cpue_data2
   }
   
   
   # Number of size breaks
   if(!mat_sex %in% c("Mature Female", "All crab")){
     num = 500
   } else{
     num = 1500
   }
   
   # Plot interpolated data
   ggplot() +
     ggplot2::geom_tile(data = temp_data, 
                        aes(x = x, 
                            y = y,
                            fill = cut(temperature, 
                                       breaks = temp_breaks))) +
     geom_sf(data = st_as_sf(CPS1_bound), fill = NA, color = "black", linewidth = 0.5)+
     geom_sf(data = st_transform(map_layers$akland, map.crs), fill = "grey80") +
     geom_sf(data = cpue_data2,
             mapping = aes(size=CPUE, shape = CPUE == 0), alpha = 0.5, fill = "black")+
     scale_shape_manual(values = c('TRUE' = 4, 'FALSE' = 21), guide = "none")+
     geom_sf(data = st_transform(map_layers$survey.area, map.crs), fill = NA, linewidth = 0.5) +
     scale_x_continuous(breaks = c(-165, -160), labels = paste0(c(165, 160), "°W"))+
     scale_y_continuous(breaks = c(56, 58), labels = paste0(c(56, 58), "°N"))+
     labs(title = paste(mat_sex))+
     coord_sf(xlim = plot.boundary$x,
              ylim = plot.boundary$y) +
     ggplot2::scale_fill_manual(name = "Temperature (°C)", values = viridis_pal(option = viridis_option)(n_temp_breaks),
                                labels = c(expression(""<=-1), "-0.9-0", "0.1-1", "1.1-2", "2.1-3", "3.1-4",
                                           "4.1-5", "5.1-6", "6.1-7", "7.1-8", ">8.1"), drop = FALSE) +
     geom_shadowtext(year_lab,
                     mapping = aes(label = lab, x = X, y = Y), size = 8,  color = "black", bg.color = "white")+
     scale_size_continuous(range = c(2, 10), limits = c(0, max(cpue_data$CPUE)), 
                           breaks =seq(0, max(cpue_data$CPUE), by = num))+
     guides(size = guide_legend(title.position = "top",ncol = 1, order = 1, 
           override.aes = list(shape = c(4, rep(21, length(seq(0, max(cpue_data$CPUE), by = num))-1)))),
            fill = "none")+
     theme_bw() +
     theme(axis.title = element_blank(),
           axis.text = element_text(size = 15),
           legend.text = element_text(size = 15),
           legend.title = element_text(size = 15),
           legend.position = "right",
           legend.direction = "horizontal",
           plot.title = element_text(face = "bold", size = 20),
           plot.subtitle = element_text(size = 12))  -> temp_cpue_map
   
   ggsave(plot = temp_cpue_map, paste0("./Figures/CPS1 Report/NMFStempcpue_", years, "_", mat_sex, ".png"), 
          height=7, width=10, units="in")
   
   return(temp_cpue_map)
 }
 
 # ALL CRAB
 temp_cpue_map(NMFS_temp22, NMFScpue_sum, 2022, "All crab") -> out
 temp_cpue_map(NMFS_temp23, NMFScpue_sum, 2023, "All crab") -> out
 
 # LEGAL MALES
 temp_cpue_map(NMFS_temp22, NMFScpue_out$legal_male_cpue, 2022, "Legal Male") -> out
 temp_cpue_map(NMFS_temp23, NMFScpue_out$legal_male_cpue, 2023, "Legal Male") -> out
 
 # MATURE MALES
 temp_cpue_map(NMFS_temp22, NMFScpue_out$mat_male_cpue, 2022, "Mature Male") -> out
 temp_cpue_map(NMFS_temp23, NMFScpue_out$mat_male_cpue, 2023, "Mature Male") -> out
 
 # IMMATURE MALES
 temp_cpue_map(NMFS_temp22, NMFScpue_out$imm_male_cpue, 2022, "Immature Male") -> out
 temp_cpue_map(NMFS_temp23, NMFScpue_out$imm_male_cpue, 2023, "Immature Male") -> out
 
 # MATURE FEMALES
 temp_cpue_map(NMFS_temp22, NMFScpue_out$mat_fem_cpue, 2022, "Mature Female") -> out
 temp_cpue_map(NMFS_temp23, NMFScpue_out$mat_fem_cpue, 2023, "Mature Female") -> out
 
 # IMMATURE FEMALES
 temp_cpue_map(NMFS_temp22, NMFScpue_out$imm_fem_cpue, 2022, "Immature Female") -> out
 temp_cpue_map(NMFS_temp23, NMFScpue_out$imm_fem_cpue, 2023, "Immature Female") -> out
 

 
 # CPS1 

  # Year labels
 sf::st_as_sf(data.frame(lab= paste("", 2023, "\nCPS1"), 
                         x = c(-160), y = c(55.2)),
              coords = c(x = "x", y = "y"), crs = sf::st_crs(4326)) %>%
   sf::st_transform(crs = map.crs) %>%
   cbind(2023, st_coordinates(.)) -> CPS1year_lab
 
 # Temp breaks
 temp_breaks <- c(-Inf, seq(-1,8,1), Inf)
 viridis_option <- "H" # viridis turbo palette
 n_temp_breaks <- length(temp_breaks)-1
 
 # Transform plot boundary
 plot.boundary.untrans <- data.frame(y = c(54.5, 58.5), 
                                     x = c(-164.8, -159))
 
 plot.boundary <-  plot.boundary.untrans %>%
   sf::st_as_sf(coords = c(x = "x", y = "y"), crs = sf::st_crs(4326)) %>%
   sf::st_transform(crs = map.crs) %>%
   sf::st_coordinates() %>%
   as.data.frame() %>%
   dplyr::rename(x = X, y = Y) # plot boundary projected
 
 breaks.x <- map_layers$lon.breaks[(map_layers$lon.breaks >= plot.boundary.untrans$x[1] &  # set lon breaks
                                      map_layers$lon.breaks < plot.boundary.untrans$x[2]) == TRUE]
 
 breaks.y <- map_layers$lat.breaks[(map_layers$lat.breaks > plot.boundary.untrans$y[1] & # set lat breaks
                                      map_layers$lat.breaks < plot.boundary.untrans$y[2]) == TRUE]
 # ALL CRAB
 ggplot() +
   ggplot2::geom_tile(data = temp_df, 
                      aes(x = x, 
                          y = y,
                          fill = cut(temperature, 
                                     breaks = temp_breaks))) +
   geom_sf(data = st_as_sf(CPS1_bound), fill = NA, color = "black", linewidth = 1)+
   #geom_sf(data = CPS1_bathy, color=alpha("white")) +
   #geom_sf(data = st_as_sf(BB_strata), fill = NA, mapping = aes(color = "black"), linewidth = 1) +
   geom_sf(data = st_transform(map_layers$akland, map.crs), fill = "grey80") +
   geom_sf(data = pot_cpue_allcrabmapdat,
           mapping = aes(size=COUNT, shape = COUNT == 0), alpha = 0.5, fill = "black")+
   scale_shape_manual(values = c('TRUE' = 4, 'FALSE' = 21), guide = "none")+
   coord_sf(xlim = plot.boundary$x,
            ylim = plot.boundary$y) +
   scale_x_continuous(breaks = c(-165, -160), labels = paste0(c(165, 160), "°W"))+
   scale_y_continuous(breaks = c(56, 58), labels = paste0(c(56, 58), "°N"))+
   scale_size_continuous(range = c(2, 10), limits = c(0, max(pot_cpue_allcrabmapdat$COUNT)), 
                         breaks =seq(0, max(pot_cpue_allcrabmapdat$COUNT), by = 50))+ 
   
   ggplot2::scale_fill_manual(name = "Temperature (°C)", values = viridis_pal(option = viridis_option)(n_temp_breaks),
                              labels = c(expression(""<=-1), "-0.9-0", "0.1-1", "1.1-2", "2.1-3", "3.1-4",
                                         "4.1-5", "5.1-6", "6.1-7", "7.1-8", ">8.1"), drop = FALSE) +
   guides(color = guide_legend(nrow = 2), size = guide_legend(title.position = "top",ncol = 1, order = 1, 
          override.aes = list(shape = c(4, rep(21, 6)))), fill = "none") +
   geom_shadowtext(CPS1year_lab,
                   mapping = aes(label = lab, x = X, y = Y), size = 8,  color = "black", bg.color = "white")+
   theme_bw() +
   theme(axis.title = element_blank(),
         axis.text = element_text(size = 15),
         legend.text = element_text(size = 15),
         legend.title = element_text(size = 15),
         legend.position = "right",
         legend.direction = "horizontal",
         plot.title = element_text(face = "bold", size = 20),
         plot.subtitle = element_text(size = 12)) -> CPS1_tempcpue_Allcrab
 
 ggsave(plot = CPS1_tempcpue_Allcrab, "./Figures/CPS1 Report/CPS1_tempcpuediscrete_Allcrab.png", height=7, width=10, units="in")
 
 # BY SEX-MATURITY CATEGORY
 mat_sex_combos %>%
   purrr::map(~ggplot() +
         ggplot2::geom_tile(data = temp_df, 
                            aes(x = x, 
                                y = y,
                                fill = cut(temperature, 
                                           breaks = temp_breaks))) +
         geom_sf(data = st_as_sf(CPS1_bound), fill = NA, color = "black", linewidth = 1)+
         #geom_sf(data = CPS1_bathy, color=alpha("white")) +
         #geom_sf(data = st_as_sf(BB_strata), fill = NA, mapping = aes(color = "black"), linewidth = 1) +
         geom_sf(data = st_transform(map_layers$akland, map.crs), fill = "grey80") +
         geom_sf(data = filter(pot_cpue_mapdat, MAT_SEX == .x),
                 mapping = aes(size=COUNT, shape = COUNT == 0), alpha = 0.5, fill = "black")+
         scale_shape_manual(values = c('TRUE' = 4, 'FALSE' = 21), guide = "none")+
         coord_sf(xlim = plot.boundary$x,
                  ylim = plot.boundary$y) +
         scale_x_continuous(breaks = c(-165, -160), labels = paste0(c(165, 160), "°W"))+
         scale_y_continuous(breaks = c(56, 58), labels = paste0(c(56, 58), "°N"))+
         scale_size_continuous(range = c(2, 10), limits = c(0, max(pot_cpue_mapdat$COUNT)), 
                               breaks =seq(0, max(pot_cpue_mapdat$COUNT), by = 50))+ 
         ggplot2::scale_fill_manual(name = "Temperature (°C)", values = viridis_pal(option = viridis_option)(n_temp_breaks),
                                    labels = c(expression(""<=-1), "-0.9-0", "0.1-1", "1.1-2", "2.1-3", "3.1-4",
                                               "4.1-5", "5.1-6", "6.1-7", "7.1-8", ">8.1"), drop = FALSE) +
         guides(color = guide_legend(nrow = 2), size = guide_legend(title.position = "top",ncol = 1, order = 1, 
                override.aes = list(shape = c(4, rep(21, length(seq(0, max(pot_cpue_mapdat$COUNT), by = 50))-1)))), 
                fill = "none") +
         geom_shadowtext(CPS1year_lab,
                         mapping = aes(label = lab, x = X, y = Y), size = 8,  color = "black", bg.color = "white")+
         theme_bw() +
           theme(axis.title = element_blank(),
                 axis.text = element_text(size = 15),
                 legend.text = element_text(size = 15),
                 legend.title = element_text(size = 15),
                 legend.position = "right",
                 legend.direction = "horizontal",
                 plot.title = element_text(face = "bold", size = 20),
                 plot.subtitle = element_text(size = 12))) -> CPS1_tempcount
 
 ggsave(plot = CPS1_tempcount[[1]], "./Figures/CPS1 Report/ CPS1_tempcpuediscrete_matmale.png", height=7, width=10, units="in")
 ggsave(plot = CPS1_tempcount[[2]], "./Figures/CPS1 Report/ CPS1_tempcpuediscrete_immmale.png", height=7, width=10, units="in")
 ggsave(plot = CPS1_tempcount[[3]], "./Figures/CPS1 Report/ CPS1_tempcpuediscrete_matfem.png", height=7, width=10, units="in")
 ggsave(plot = CPS1_tempcount[[4]], "./Figures/CPS1 Report/ CPS1_tempcpuediscrete_immfem.png", height=7, width=10, units="in")
 ggsave(plot = CPS1_tempcount[[5]], "./Figures/CPS1 Report/ CPS1_tempcpuediscrete_legmale.png", height=7, width=10, units="in")
 
# REQUEST 8A: NMFS and CPS1 FEMALE SHELL, CLUTCH, EGG BY 1MM BINS -----------------------------------------------------------------
 
  # Load function to generate NMFS female data
  fem_clutch_shell_data <- function(data_crab, data_strata, stock, years){
   
   #District area lookup table by stock, pull districts by stock specified
   data.frame(STOCK = c("BBRKC", 
                        rep("PribRKC", 2),
                        "NSRKC",
                        rep("PribBKC", 2), 
                        rep("StMattBKC", 2),
                        "BKCNBS",
                        rep("TannerW", 3), 
                        "TannerE", 
                        rep("Snow", 3), 
                        "SnowNBS",
                        rep("Hybrid", 3), 
                        rep("Hair", 4), 
                        rep("Allstations", 10),
                        rep("NorthRKC",1)),
              DISTRICT = c("Bristol Bay", 
                           "Pribilof MTCA", "Pribilof Single",
                           "Norton Sound",
                           "Pribilof MTCA", "Pribilof Single",
                           "St. Matthew MTCA", "St. Matthew Single",
                           "NBS All",
                           "Pribilof MTCA", "St. Matthew MTCA", "West 166",
                           "East 166",
                           "Pribilof MTCA", "Single", "St. Matthew MTCA",
                           "NBS All",
                           "Pribilof MTCA", "Single", "St. Matthew MTCA",
                           "Bristol Bay", "Northern Unstratified", "Pribilof MTCA", "Pribilof Single",
                           "Bristol Bay", "Northern Unstratified", "Pribilof MTCA","Pribilof Single",      
                           "BKC Unstratified", "St. Matthew MTCA", "St. Matthew Single", "East 166",             
                           "West 166", "Single",("Northern Unstratified"))) %>%
     dplyr::filter(STOCK == stock) %>%
     pull(DISTRICT) -> dist
   
   #Pull stock stations from strata tables using stock districts
   data_strata %>%
     dplyr::filter(DISTRICT %in% dist) %>%
     pull(STATION_ID) %>%
     unique() -> stock_stations
   
   #Specify retow stations for BBRKC, pull by year
   retow_sta <- read.csv("./Data/retow_sta.csv")
   
   retow_sta %>%
     dplyr::filter(SURVEY_YEAR %in% years) %>%
     dplyr::select(GIS_STATION) %>%
     pull() -> retow_stations
   
   #Create female shell condition, clutch & egg levels 
   data_crab %>%
     dplyr::filter(SEX == 2,
                   AKFIN_SURVEY_YEAR %in% years)  %>%
     mutate(SHELL_TEXT = case_when(SHELL_CONDITION %in% 0:1 ~ "Soft Molting",
                                   SHELL_CONDITION == 2 ~ "New Hard",
                                   SHELL_CONDITION == 3 ~ "Old",
                                   SHELL_CONDITION %in% 4:5 ~ "Very Old"),
            EGG_CONDITION_TEXT = case_when(EGG_CONDITION == 0 ~ "No Eggs",
                                           EGG_CONDITION == 1 ~ "Uneyed Eggs",
                                           EGG_CONDITION == 2 ~ "Eyed Eggs",
                                           EGG_CONDITION == 3 ~ "Dead Eggs",
                                           EGG_CONDITION == 4 ~ "Empty Egg Cases",
                                           EGG_CONDITION == 5 ~ "Hatching",
                                           ((CLUTCH_SIZE == 999 & EGG_CONDITION == 0)| is.na(EGG_CONDITION) == TRUE) ~ "Unknown"),
            CLUTCH_TEXT = case_when(CLUTCH_SIZE == 0 ~ "Immature",
                                    CLUTCH_SIZE == 1 ~ "Mature Barren",
                                    CLUTCH_SIZE == 2 ~ "Trace",
                                    CLUTCH_SIZE == 3 ~ "Quarter Full",
                                    CLUTCH_SIZE == 4 ~ "Half Full",
                                    CLUTCH_SIZE == 5 ~ "Three Quarter Full",
                                    CLUTCH_SIZE == 6 ~ "Full",
                                    ((CLUTCH_SIZE == 999) | is.na(CLUTCH_SIZE) == TRUE) ~ "Unknown")) -> female_ebs_spec
   
   
   #Compute 6-year mean cpue by 1mm size bin for shell/clutch categories  
   female_ebs_spec %>%
     dplyr::filter(GIS_STATION %in% stock_stations) %>%
     #create 1mm size bin
     mutate(SIZE_BIN = case_when(SPECIES_CODE %in% c(69322,69323) ~ floor(LENGTH_1MM),  
                                 SPECIES_CODE %in% c(68560,68580) ~ floor(WIDTH_1MM)),
            CPUE = SAMPLING_FACTOR/AREA_SWEPT) %>% 
     group_by(AKFIN_SURVEY_YEAR, HAUL_TYPE, GIS_STATION, SIZE_BIN, SHELL_TEXT, EGG_CONDITION_TEXT, CLUTCH_TEXT) %>%
     dplyr::summarise(CPUE = sum(CPUE)) -> female_cpue 
   
   
   #Pull out unique haul types by year in order to include in zero-catch station dataframe below
   data_crab%>%
     dplyr::filter(AKFIN_SURVEY_YEAR %in% years) %>%
     distinct(HAUL_TYPE)->HT
   
   #Join to zero catch stations
   female_cpue %>%
     right_join(expand_grid(SIZE_BIN = 1:180, 
                            HAUL_TYPE = as.numeric(levels(factor(HT$HAUL_TYPE))),
                            data_strata %>%
                              dplyr::filter(SURVEY_YEAR %in% years,
                                            STATION_ID %in% stock_stations) %>% 
                              distinct(SURVEY_YEAR, STATION_ID, STRATUM, TOTAL_AREA) %>%
                              rename_all(~c("AKFIN_SURVEY_YEAR", "GIS_STATION",
                                            "STRATUM", "TOTAL_AREA"))) %>%
                  expand_grid(female_ebs_spec %>%
                                distinct(SHELL_TEXT, EGG_CONDITION_TEXT, CLUTCH_TEXT))) %>%
     replace_na(list(CPUE = 0)) %>%
     group_by(AKFIN_SURVEY_YEAR, GIS_STATION) %>%
     nest() %>%
     #Replacing original stations with resampled stations in BB retow yrs
     mutate(data = purrr::map(data, function(x) {
       if(17 %in% x$HAUL_TYPE & stock == "BBRKC" & GIS_STATION %in% retow_stations)  
       {x %>% dplyr::filter(HAUL_TYPE == 17) -> x} else {x %>% dplyr::filter(HAUL_TYPE != 17) -> x}
       return(x)  
     })) %>%
     
     unnest()%>%
     #sum across haul types 
     group_by(AKFIN_SURVEY_YEAR, GIS_STATION, STRATUM, SIZE_BIN, TOTAL_AREA, SHELL_TEXT, EGG_CONDITION_TEXT, CLUTCH_TEXT) %>%
     dplyr::reframe(CPUE = sum(CPUE)) %>%
     
     #Scale to abundance by strata
     group_by(AKFIN_SURVEY_YEAR, GIS_STATION, STRATUM, SIZE_BIN, TOTAL_AREA, SHELL_TEXT, EGG_CONDITION_TEXT, CLUTCH_TEXT) %>%
     dplyr::reframe(CPUE = sum(CPUE)) %>%
     group_by(AKFIN_SURVEY_YEAR, STRATUM, SIZE_BIN, TOTAL_AREA, SHELL_TEXT, EGG_CONDITION_TEXT, CLUTCH_TEXT) %>%
     dplyr::reframe( MEAN_CPUE = mean(CPUE, na.rm = T),
                     ABUNDANCE = (MEAN_CPUE * TOTAL_AREA)) %>%
     distinct() %>%
     group_by(AKFIN_SURVEY_YEAR, SIZE_BIN, SHELL_TEXT, EGG_CONDITION_TEXT, CLUTCH_TEXT) %>%
     #Sum across strata
     dplyr::reframe(ABUNDANCE_MIL = sum(ABUNDANCE)/1e6) %>%
     dplyr::filter(is.na(CLUTCH_TEXT) == FALSE) -> female_abundance
   
   return(female_abundance = female_abundance)
 }
 
  # Load function to generate NMFS plots
  fem_shell_clutch_tabplot <- function(female_abundance, years, stock) {
     
     #Create lookup table for plots
     plot_lookup <- tibble(stock = c("BBRKC", "PribRKC", "PribBKC", "StMattBKC", "TannerE", "TannerW", "Snow","NSRKC", "BKCNBS", "SnowNBS"),
                           title = c("Female Bristol Bay Red King Crab", "Female Pribilof Islands Red King Crab", "Female Pribilof Islands Blue King Crab",
                                     "Female St. Matthew Island Blue King Crab", "Female Tanner Crab East", "Female Tanner Crab West", "Female Snow Crab",
                                     "Female Norton Sound Red King Crab","Female Northern Bering Sea Blue King Crab","Female Northern Bering Sea Snow Crab"),
                           xlab = c("Carapace length (mm)","Carapace length (mm)","Carapace length (mm)","Carapace length (mm)",
                                    "Carapace width (mm)","Carapace width (mm)","Carapace width (mm)","Carapace length (mm)","Carapace length (mm)","Carapace width (mm)"),
                           XMax = c(190,190,175,135,125,115,90,115,135,75)) 
     
     design <- "
    AB
   "
       
     #Specify plot labels
     labs <- paste("NMFS", years)
     
     names(labs) <- years
     
     
     #Female Mature vs Immature
     female_abundance%>%
       dplyr::filter(AKFIN_SURVEY_YEAR %in% years, ABUNDANCE_MIL>0) %>%
       mutate(Maturity = ifelse(CLUTCH_TEXT == "Unknown", "Unknown", 
                                ifelse(CLUTCH_TEXT=="Immature","Immature","Mature"))) %>%
       mutate(Maturity= factor(Maturity, levels = c("Immature", "Mature", "Unknown"))) %>%
       
       ggplot(aes(x=SIZE_BIN, y=ABUNDANCE_MIL, fill=Maturity)) +
       geom_bar(position="stack", stat="identity") +
       theme_bw() +
       scale_x_continuous(expand = c(0.01,0.01), limits = c(30, 175))+ # same as CPS1/NMFS
       labs(fill = "Maturity", title=plot_lookup$title[plot_lookup$stock == stock]) +
       theme(plot.title=element_text(hjust=0.5)) +
       ylab("Abundance (millions)") +
       xlab(plot_lookup$xlab[plot_lookup$stock == stock]) +
       scale_fill_manual(name = "Maturity", values = c("#5ab4ac","#d8b365","#252525"),
                         labels = c("Immature", "Mature", "Unknown"), drop = FALSE)+
       facet_manual(~AKFIN_SURVEY_YEAR, design, scales = "free_y", labeller = labeller(AKFIN_SURVEY_YEAR = labs)) +
       theme(legend.position = "bottom", legend.background=element_rect(fill="white", color="black"))->maturity_plot
     
       
     
     
     ggsave(paste0("./Figures/CPS1 Report/NMFSfemale_Maturity_", stock, ".png"), plot=maturity_plot, height=8.5, width=7.5, units="in")
     
     
     #Shell condition figure
     female_abundance %>%
       dplyr::filter(AKFIN_SURVEY_YEAR %in% years, ABUNDANCE_MIL>0,CLUTCH_TEXT!="Immature") %>% 
       mutate(SHELL_TEXT  = factor(SHELL_TEXT , levels = c('Very Old', 'Old', 'New Hard', 'Soft Molting'))) %>%
       ggplot(aes(x=SIZE_BIN, y=ABUNDANCE_MIL, fill=SHELL_TEXT)) +
       geom_bar(position="stack", stat="identity") +
       theme_bw() +
       scale_x_continuous(expand = c(0.01,0.01), limits = c(80, 175))+ # same as CPS1/NMFS
       labs(fill = "Shell Condition", title=paste("Mature", plot_lookup$title[plot_lookup$stock == stock])) +
       theme(plot.title=element_text(hjust=0.5)) +
       ylab("Abundance (millions)") +
       xlab(plot_lookup$xlab[plot_lookup$stock == stock]) +
       # scale_fill_manual(name = "Shell Condition", values = c("#999999", "#E69F00", "#56B4E9", "#009E73"),
       #                  labels = c("Very Old", "Old", "New Hard", "Soft Molting"), drop = FALSE)+
       scale_fill_manual(name = "Shell Condition", values = c("#018571", "#80cdc1", "#dfc27d", "#a6611a"),
                         labels = c("Very Old", "Old", "New Hard", "Soft Molting"), drop = FALSE)+
       facet_manual(~AKFIN_SURVEY_YEAR, design, scales = "free_y", labeller = labeller(AKFIN_SURVEY_YEAR = labs)) +
       theme(legend.position = "bottom", legend.background=element_rect(fill="white", color="black"))->shell_plot
     
     
     
     ggsave(paste0("./Figures/CPS1 Report/NMFSfemale_SC_", stock, ".png"), plot=shell_plot, height=8.5, width=7.5, units="in")
     
     #Egg condition figure
     female_abundance %>%
       dplyr::filter(AKFIN_SURVEY_YEAR %in% years, ABUNDANCE_MIL>0, CLUTCH_TEXT!="Immature") %>%
       mutate(EGG_CONDITION_TEXT = ifelse(is.na(EGG_CONDITION_TEXT) == "TRUE", "Unknown", EGG_CONDITION_TEXT)) %>%
       mutate(EGG_CONDITION_TEXT = factor(EGG_CONDITION_TEXT, levels = c('Uneyed Eggs','No Eggs','Dead Eggs','Hatching','Empty Egg Cases', 
                                                                         'Eyed Eggs', "Unknown"))) %>%
       ggplot(aes(x=SIZE_BIN, y=ABUNDANCE_MIL, fill=EGG_CONDITION_TEXT)) +
       geom_bar(position="stack", stat="identity") +
       theme_bw() +
       scale_x_continuous(expand = c(0.01,0.01), limits = c(80, 175))+ # same as CPS1/NMFS
       labs(fill = "Egg Condition",title=paste("Mature", plot_lookup$title[plot_lookup$stock == stock])) +
       theme(plot.title=element_text(hjust=0.5)) +
       ylab("Abundance (millions)") +
       xlab(plot_lookup$xlab[plot_lookup$stock == stock]) +
       scale_fill_manual(name = "Egg Condition", values=c("#01665e","#5ab4ac","#c7eae5","#f6e8c3","#d8b365","#8c510a", 
                                                          "#252525"), labels = c('Uneyed Eggs','No Eggs','Dead Eggs','Hatching','Empty Egg Cases', 
                                                                                 'Eyed Eggs', "Unknown"), drop = FALSE) +
       facet_manual(~AKFIN_SURVEY_YEAR, design, scales = "free_y", labeller = labeller(AKFIN_SURVEY_YEAR = labs)) +
       theme(legend.position = "bottom", legend.background=element_rect(fill="white", color="black"))-> egg_plot
     
     ggsave(paste0("./Figures/CPS1 Report/NMFSfemale_Egg_", stock, ".png"), plot=egg_plot, height=8.5, width=7.5, units="in")
     
     #Clutch size figure
     female_abundance %>%
       dplyr::filter(AKFIN_SURVEY_YEAR %in% years, ABUNDANCE_MIL>0,CLUTCH_TEXT!="Immature") %>%
       mutate(CLUTCH_TEXT = ifelse(is.na(CLUTCH_TEXT) == "TRUE", "Unknown", CLUTCH_TEXT)) %>%
       mutate(CLUTCH_TEXT= factor(CLUTCH_TEXT, levels = c('Full', 'Three Quarter Full', 'Half Full', 
                                                          'Quarter Full', 'Trace', 'Mature Barren', "Unknown"))) %>%
       ggplot(aes(x=SIZE_BIN, y=ABUNDANCE_MIL, fill=CLUTCH_TEXT)) +
       scale_x_continuous(expand = c(0.01,0.01), limits = c(80, 175))+ # same as CPS1/NMFS
       geom_bar(position="stack", stat="identity") +
       theme_bw() +
       labs(fill = "Clutch Size",title=paste("Mature", plot_lookup$title[plot_lookup$stock == stock])) +
       theme(plot.title=element_text(hjust=0.5)) +
       ylab("Abundance (millions)") +
       xlab(plot_lookup$xlab[plot_lookup$stock == stock]) +
       # scale_fill_manual(name = "Clutch Size", values=c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", 
       #                  "#003366", "#CC79A7"), labels =c('Immature', "Full", 'Three Quarter Full', 'Half Full', 
       #                  'Quarter Full', 'Trace', 'Mature Barren', "Unknown"), drop = FALSE) +
       scale_fill_manual(name = "Clutch Size", values=c("#01665e", "#5ab4ac", "#c7eae5", "#f6e8c3", "#d8b365", 
                                                        "#8c510a", "#252525"), labels =c("Full", 'Three Quarter Full', 'Half Full', 
                                                                                         'Quarter Full', 'Trace', 'Mature Barren', "Unknown"), drop = FALSE) +
       facet_manual(~AKFIN_SURVEY_YEAR, design, scales = "free_y", labeller = labeller(AKFIN_SURVEY_YEAR = labs)) +
       theme(legend.position = "bottom", legend.background=element_rect(fill="white", color="black"))-> clutch_plot
     
     ggsave(paste0("./Figures/CPS1 Report/Female_Clutch_", stock, ".png"), plot=clutch_plot, height=8.5, width=7.5, units="in")
     
     return(list(female_abundance = female_abundance, shell_plot = shell_plot, egg_plot = egg_plot, clutch_plot = clutch_plot, maturity_plot=maturity_plot))
  }
  
  # Run functions for NMFS female data and plots
  2022:2023 %>%
    purrr::map_df(~fem_clutch_shell_data(catch_rkc, strata_rkc, "BBRKC", .x)) -> female_abundance
  
  fem_shell_clutch_tabplot(female_abundance, 2022:2023, "BBRKC") -> NMFS_fem_out
  
  # Calculate CPS1 female clutch, shell, egg condition by 1mm bin 
    #Create female shell condition, clutch & egg levels 
    specimen_table %>%
      mutate(YEAR = 2023) %>%
      dplyr::filter(SEX == 2) %>%
      mutate(SHELL_TEXT = case_when(SHELL_CONDITION %in% 0:1 ~ "Soft Molting",
                                    SHELL_CONDITION == 2 ~ "New Hard",
                                    SHELL_CONDITION == 3 ~ "Old",
                                    SHELL_CONDITION %in% 4:5 ~ "Very Old"),
             EGG_CONDITION_TEXT = case_when(EGG_CONDITION == 0 ~ "No Eggs",
                                            EGG_CONDITION == 1 ~ "Uneyed Eggs",
                                            EGG_CONDITION == 2 ~ "Eyed Eggs",
                                            EGG_CONDITION == 3 ~ "Dead Eggs",
                                            EGG_CONDITION == 4 ~ "Empty Egg Cases",
                                            EGG_CONDITION == 5 ~ "Hatching",
                                            ((CLUTCH_SIZE == 999 & EGG_CONDITION == 0)| is.na(EGG_CONDITION) == TRUE) ~ "Unknown"),
             CLUTCH_TEXT = case_when(CLUTCH_SIZE == 0 ~ "Immature",
                                     CLUTCH_SIZE == 1 ~ "Mature Barren",
                                     CLUTCH_SIZE == 2 ~ "Trace",
                                     CLUTCH_SIZE == 3 ~ "Quarter Full",
                                     CLUTCH_SIZE == 4 ~ "Half Full",
                                     CLUTCH_SIZE == 5 ~ "Three Quarter Full",
                                     CLUTCH_SIZE == 6 ~ "Full",
                                     ((CLUTCH_SIZE == 999) | is.na(CLUTCH_SIZE) == TRUE) ~ "Unknown")) -> female_CPS1_spec
    #Compute cpue by 1mm size bin for shell/clutch categories  
    female_CPS1_spec %>%
      #create 1mm size bin
      mutate(SIZE_BIN = floor(LENGTH),
             CATCH_PER_HOUR = SAMPLING_FACTOR/SOAK_TIME) %>% 
      dplyr::group_by(YEAR, VESSEL, SPN, POT_ID, BUOY, LAT_DD, LON_DD, SIZE_BIN, 
                      SHELL_TEXT, EGG_CONDITION_TEXT, CLUTCH_TEXT) %>%
      dplyr::reframe(COUNT = sum(SAMPLING_FACTOR),
                     CATCH_PER_HOUR = sum(CATCH_PER_HOUR)) -> positive_fem_catch
    
    
    #Join to zero catch stations
    positive_fem_catch %>%
      dplyr::right_join(expand_grid(SIZE_BIN = 1:max(positive_fem_catch$SIZE_BIN),
                                    potlifts) %>%
                          expand_grid(female_CPS1_spec %>%
                                        distinct(SHELL_TEXT, EGG_CONDITION_TEXT, CLUTCH_TEXT))) %>%
      replace_na(list(COUNT = 0, CATCH_PER_HOUR = 0)) %>%
      group_by(YEAR, SIZE_BIN, SHELL_TEXT, EGG_CONDITION_TEXT, CLUTCH_TEXT) %>%
      reframe(COUNT = sum(COUNT)) %>%
      mutate(YEAR = 2023) -> CPS1_fem_1mm
    
    #Specify plot labels
    labs <- paste("CPS1", 2023)
    
    names(labs) <- 2023
    
    
    #Female Mature vs Immature
    CPS1_fem_1mm %>%
      dplyr::filter(COUNT>0) %>%
      mutate(Maturity = ifelse(CLUTCH_TEXT == "Unknown", "Unknown", 
                               ifelse(CLUTCH_TEXT=="Immature","Immature","Mature"))) %>%
      mutate(Maturity= factor(Maturity, levels = c("Immature", "Mature", "Unknown"))) %>%
      
      ggplot(aes(x=SIZE_BIN, y=COUNT, fill=Maturity)) +
      geom_bar(position="stack", stat="identity") +
      theme_bw() +
      #xlim(0,plot_lookup$XMax[plot_lookup$stock == stock])+
      scale_x_continuous(expand = c(0.01,0.01), limits = c(30, 175))+ # same as CPS1/NMFS
      ylab("Abundance") +
      xlab("Carapace length (mm) ") +
      scale_fill_manual(name = "Maturity", values = c("#5ab4ac","#d8b365","#252525"),
                        labels = c("Immature", "Mature", "Unknown"), drop = FALSE)+
      facet_manual(~YEAR, "A", scales = "free_y", labeller = labeller(YEAR = labs)) +
      theme(legend.position = "bottom", 
            legend.background=element_rect(fill="white", color="black"))-> CPS1_maturity_plot
    
    # Shell condition
    CPS1_fem_1mm %>%
      dplyr::filter(COUNT>0,CLUTCH_TEXT!="Immature") %>% 
      mutate(SHELL_TEXT  = factor(SHELL_TEXT , levels = c('Very Old', 'Old', 'New Hard', 'Soft Molting'))) %>%
      ggplot(aes(x=SIZE_BIN, y=COUNT, fill=SHELL_TEXT))+
      geom_bar(position="stack", stat="identity") +
      theme_bw() +
      #xlim(0,plot_lookup$XMax[plot_lookup$stock == stock])+
      scale_x_continuous(expand = c(0.01,0.01), limits = c(80, 175))+ # same as CPS1/NMFS
      ylab("Abundance") +
      xlab("Carapace length (mm) ") +
      scale_fill_manual(name = "Shell Condition", values = c("#018571", "#80cdc1", "#dfc27d", "#a6611a"),
                        labels = c("Very Old", "Old", "New Hard", "Soft Molting"), drop = FALSE)+
      facet_manual(~YEAR, "A", scales = "free_y", labeller = labeller(YEAR = labs)) +
      theme(legend.position = "bottom", 
            legend.background=element_rect(fill="white", color="black"))-> CPS1_shell_plot
    
    # Egg condition
    CPS1_fem_1mm %>%
      dplyr::filter(COUNT>0, CLUTCH_TEXT!="Immature") %>%
      mutate(EGG_CONDITION_TEXT = ifelse(is.na(EGG_CONDITION_TEXT) == "TRUE", "Unknown", EGG_CONDITION_TEXT)) %>%
      mutate(EGG_CONDITION_TEXT = factor(EGG_CONDITION_TEXT, levels = c('Uneyed Eggs','No Eggs','Dead Eggs','Hatching','Empty Egg Cases', 
                                                                        'Eyed Eggs', "Unknown"))) %>%
      ggplot(aes(x=SIZE_BIN, y=COUNT, fill=EGG_CONDITION_TEXT)) +
      geom_bar(position="stack", stat="identity") +
      theme_bw() +
      #xlim(0,plot_lookup$XMax[plot_lookup$stock == stock])+
      scale_x_continuous(expand = c(0.01,0.01), limits = c(80, 175))+ # same as CPS1/NMFS
      ylab("Abundance") +
      xlab("Carapace length (mm) ") +
      scale_fill_manual(name = "Egg Condition", values=c("#01665e","#5ab4ac","#c7eae5","#f6e8c3","#d8b365","#8c510a", 
                                                         "#252525"), labels = c('Uneyed Eggs','No Eggs','Dead Eggs','Hatching','Empty Egg Cases', 
                                                                                'Eyed Eggs', "Unknown"), drop = FALSE) +
      facet_manual(~YEAR, "A", scales = "free_y", labeller = labeller(YEAR = labs)) +
      theme(legend.position = "bottom", 
            legend.background=element_rect(fill="white", color="black"))-> CPS1_egg_plot
    
    # Clutch code
    CPS1_fem_1mm %>%
      dplyr::filter(COUNT>0,CLUTCH_TEXT!="Immature") %>%
      mutate(CLUTCH_TEXT = ifelse(is.na(CLUTCH_TEXT) == "TRUE", "Unknown", CLUTCH_TEXT)) %>%
      mutate(CLUTCH_TEXT= factor(CLUTCH_TEXT, levels = c('Full', 'Three Quarter Full', 'Half Full', 
                                                         'Quarter Full', 'Trace', 'Mature Barren', "Unknown"))) %>%
      ggplot(aes(x=SIZE_BIN, y=COUNT, fill=CLUTCH_TEXT)) +
      geom_bar(position="stack", stat="identity") +
      theme_bw() +
      #xlim(0,plot_lookup$XMax[plot_lookup$stock == stock])+
      scale_x_continuous(expand = c(0.01,0.01), limits = c(80, 175))+ # same as CPS1/NMFS
      ylab("Abundance") +
      xlab("Carapace length (mm) ") +
      scale_fill_manual(name = "Clutch Size", values=c("#01665e", "#5ab4ac", "#c7eae5", "#f6e8c3", "#d8b365", 
                                                       "#8c510a", "#252525"), labels =c("Full", 'Three Quarter Full', 'Half Full', 
                                                                                        'Quarter Full', 'Trace', 'Mature Barren', "Unknown"), drop = FALSE) +
      facet_manual(~YEAR, "A", scales = "free_y", labeller = labeller(YEAR = labs)) +
      theme(legend.position = "bottom", 
            legend.background=element_rect(fill="white", color="black"))-> CPS1_clutch_plot
    
  
  
 
  # Combine NMFS/CPS1 plots, save 
    # Maturity
      ggarrange(ggarrange(NMFS_fem_out$maturity_plot+ rremove("xlab"),
                          legend = "none"),
                CPS1_maturity_plot, heights = c(1,1), nrow= 2) 
      
      ggsave("./Figures/CPS1 Report/CPS1_maturity_panels.png",
             height=8.5, width=8.5, units="in")
      
    # Shell condition
      ggarrange(ggarrange(NMFS_fem_out$shell_plot+ rremove("xlab"),
                          legend = "none"),
                CPS1_shell_plot, heights = c(1,1), nrow= 2) 
      
      ggsave("./Figures/CPS1 Report/CPS1_shell_panels.png",
             height=8.5, width=8.5, units="in")
      
    # Egg condition
      ggarrange(ggarrange(NMFS_fem_out$egg_plot+ rremove("xlab"),
                          legend = "none"),
                CPS1_egg_plot, heights = c(1,1), nrow= 2) 
      
      ggsave("./Figures/CPS1 Report/CPS1_egg_panels.png",
             height=8.5, width=8.5, units="in")
      
    # Clutch code
      ggarrange(ggarrange(NMFS_fem_out$clutch_plot+ rremove("xlab"),
                          legend = "none"),
                CPS1_clutch_plot, heights = c(1,1), nrow= 2) 
      
      ggsave("./Figures/CPS1 Report/CPS1_clutch_panels.png",
             height=8.5, width=8.5, units="in")
    
    
    
    
    
# REQUEST 8B: NMFS and CPS1 MALE SHELL CONDITION BY 1MM BINS ------------------------------------------------------------------
      
  # Load function to generate NMFS male data
    male_shell_data <- function(data_crab, data_strata, stock, years){
        
        #District area lookup table by stock, pull districts by stock specified
        data.frame(STOCK = c("BBRKC", 
                             rep("PribRKC", 2),
                             "NSRKC",
                             rep("PribBKC", 2), 
                             rep("StMattBKC", 2),
                             "BKCNBS",
                             rep("TannerW", 3), 
                             "TannerE", 
                             rep("Snow", 3), 
                             "SnowNBS",
                             rep("Hybrid", 3), 
                             rep("Hair", 4), 
                             rep("Allstations", 10),
                             rep("NorthRKC",1)),
                   DISTRICT = c("Bristol Bay", 
                                "Pribilof MTCA", "Pribilof Single",
                                "Norton Sound",
                                "Pribilof MTCA", "Pribilof Single",
                                "St. Matthew MTCA", "St. Matthew Single",
                                "NBS All",
                                "Pribilof MTCA", "St. Matthew MTCA", "West 166",
                                "East 166",
                                "Pribilof MTCA", "Single", "St. Matthew MTCA",
                                "NBS All",
                                "Pribilof MTCA", "Single", "St. Matthew MTCA",
                                "Bristol Bay", "Northern Unstratified", "Pribilof MTCA", "Pribilof Single",
                                "Bristol Bay", "Northern Unstratified", "Pribilof MTCA","Pribilof Single",      
                                "BKC Unstratified", "St. Matthew MTCA", "St. Matthew Single", "East 166",             
                                "West 166", "Single",("Northern Unstratified"))) %>%
          dplyr::filter(STOCK == stock) %>%
          pull(DISTRICT) -> dist
        
        #Pull stock stations from strata tables using stock districts
        data_strata %>%
          dplyr::filter(DISTRICT %in% dist) %>%
          pull(STATION_ID) %>%
          unique() -> stock_stations
        
        #Specify retow stations for BBRKC, pull by year
        retow_sta <- read.csv("./Data/retow_sta.csv")
        
        retow_sta %>%
          dplyr::filter(SURVEY_YEAR %in% years) %>%
          dplyr::select(GIS_STATION) %>%
          pull() -> retow_stations
        
        #Create male shell condition levels 
        data_crab %>%
          filter(SEX == 1,
                 AKFIN_SURVEY_YEAR %in% years) %>%
          mutate(SHELL_TEXT = case_when(SHELL_CONDITION %in% 0:1 ~ "Soft Molting",
                                        SHELL_CONDITION == 2 ~ "New Hard",
                                        SHELL_CONDITION == 3 ~ "Old",
                                        SHELL_CONDITION %in% 4:5 ~ "Very Old"))-> male_ebs_spec
        
        #Compute 6-year mean cpue by 1mm size bin for shell categories  
        male_ebs_spec %>%
          filter(GIS_STATION %in% stock_stations) %>%
          #create 1mm size bin
          mutate(SIZE_BIN = case_when(SPECIES_CODE %in% c(69322,69323) ~ floor(LENGTH_1MM),  
                                      SPECIES_CODE %in% c(68560,68580) ~ floor(WIDTH_1MM)),
                 CPUE = SAMPLING_FACTOR/AREA_SWEPT) %>% 
          group_by(AKFIN_SURVEY_YEAR, GIS_STATION, HAUL_TYPE, SIZE_BIN, SHELL_TEXT) %>%
          dplyr::summarise(CPUE = sum(CPUE)) -> male_cpue 
        
        #Pull out unique haul types by year in order to include in zero-catch station dataframe below
        data_crab%>%
          dplyr::filter(AKFIN_SURVEY_YEAR %in% years) %>%
          
          distinct(HAUL_TYPE)  ->HT
        
        #Join to zero catch stations
        male_cpue %>%
          right_join(expand_grid(SIZE_BIN = 1:max(male_cpue$SIZE_BIN),
                                 
                                 HAUL_TYPE = as.numeric(levels(factor(HT$HAUL_TYPE))),
                                 data_strata %>%
                                   filter(SURVEY_YEAR %in% years,
                                          STATION_ID %in% stock_stations) %>% 
                                   distinct(SURVEY_YEAR, STATION_ID, STRATUM, TOTAL_AREA) %>%
                                   rename_all(~c("AKFIN_SURVEY_YEAR", "GIS_STATION",
                                                 "STRATUM", "TOTAL_AREA"))) %>%
                       expand_grid(male_ebs_spec %>%
                                     distinct(SHELL_TEXT))) %>%
          replace_na(list(CPUE = 0)) %>%
          
          filter(HAUL_TYPE != 17) %>%
          
          group_by(AKFIN_SURVEY_YEAR, GIS_STATION, STRATUM, SIZE_BIN, TOTAL_AREA, SHELL_TEXT) %>%
          dplyr::summarise(CPUE = sum(CPUE)) %>%
          #Scale to abundance by strata
          group_by(AKFIN_SURVEY_YEAR, STRATUM, SIZE_BIN, TOTAL_AREA, SHELL_TEXT) %>%
          summarise(MEAN_CPUE = mean(CPUE, na.rm = T),
                    ABUNDANCE = (MEAN_CPUE * mean(TOTAL_AREA))) %>%
          group_by(AKFIN_SURVEY_YEAR, SIZE_BIN, SHELL_TEXT) %>%
          #Sum across strata
          summarise(ABUNDANCE_MIL = sum(ABUNDANCE)/1e6) -> male_abundance
        
        return(male_abundance = male_abundance)
      }
      
  # Load function to generate NMFS plots
   male_shell_tabplot <- function(male_abundance, years, stock) {
        
        #Create lookup table for plots
        plot_lookup <- tibble(stock = c("BBRKC", "PribRKC", "PribBKC", "StMattBKC", "TannerE", "TannerW", "Snow","NSRKC", "BKCNBS", "SnowNBS"),
                              title = c("Male Bristol Bay Red King Crab", "Male Pribilof Islands Red King Crab", "Male Pribilof Islands Blue King Crab",
                                        "Male St. Matthew Island Blue King Crab", "Male Tanner Crab East", "Male Tanner Crab West", "Male Snow Crab",
                                        "Male Norton Sound Red King Crab","Male Northern Bering Sea Blue King Crab","Male Northern Bering Sea Snow Crab"),
                              xlab = c("Carapace length (mm)","Carapace length (mm)","Carapace length (mm)","Carapace length (mm)",
                                       "Carapace width (mm)","Carapace width (mm)","Carapace width (mm)","Carapace length (mm)","Carapace length (mm)","Carapace width (mm)"),
                              XMax = c(220,220,175,175,190,190,140,165,150,115)) 
        
        design <- "
    AB
   "
        
        #Specify plot labels
        labs <- paste("NMFS", years)
        
        names(labs) <- years
        
        
        #Shell condition figure
        male_abundance %>%
          filter(AKFIN_SURVEY_YEAR %in% years, ABUNDANCE_MIL>0) %>%
          mutate(SHELL_TEXT  = factor(SHELL_TEXT , levels = c('Very Old', 'Old', 'New Hard', 'Soft Molting'))) %>%
          ggplot(aes(x=SIZE_BIN, y=ABUNDANCE_MIL, fill=SHELL_TEXT)) +
          geom_bar(position="stack", stat="identity") +
          theme_bw() +
          labs(fill = "Shell Condition", title=plot_lookup$title[plot_lookup$stock == stock]) +
          #xlim(0,plot_lookup$XMax[plot_lookup$stock == stock])+
          scale_x_continuous(expand = c(0.01,0.01), limits = c(16, 192))+ #same as CPS1/NMFS
          theme(plot.title=element_text(hjust=0.5)) +
          ylab("Abundance (millions)") +
          xlab(plot_lookup$xlab[plot_lookup$stock == stock]) +
          scale_fill_manual(name = "Shell Condition", values = c("#018571", "#80cdc1", "#dfc27d", "#a6611a"),
                            labels = c("Very Old", "Old", "New Hard", "Soft Molting"), drop = FALSE) +
          facet_manual(~AKFIN_SURVEY_YEAR, design, scales = "free_y", labeller = labeller(AKFIN_SURVEY_YEAR = labs)) +
          theme(legend.position = "bottom", legend.background=element_rect(fill="white", color="black"))-> shell_plot_male
        
        ggsave(paste0("./Figures/Male_SC_", stock, ".png"), plot=shell_plot_male, height=8.5, width=7.5, units="in")
        
        return(list(male_abundance = male_abundance, shell_plot = shell_plot_male))
      }
      
  # Run functions for NMFS female data and plots
    2022:2023 %>%
      purrr::map_df(~male_shell_data(catch_rkc, strata_rkc, "BBRKC", .x)) -> male_abundance
      
    male_shell_tabplot(male_abundance, 2022:2023, "BBRKC") -> NMFS_male_out
      
  # Calculate CPS1 male shell condition by 1mm bin 
    #Create male shell condition, clutch & egg levels 
      specimen_table %>%
        mutate(YEAR = 2023)%>%
        filter(SEX == 1) %>%
        mutate(SHELL_TEXT = case_when(SHELL_CONDITION %in% 0:1 ~ "Soft Molting",
                                      SHELL_CONDITION == 2 ~ "New Hard",
                                      SHELL_CONDITION == 3 ~ "Old",
                                      SHELL_CONDITION %in% 4:5 ~ "Very Old"))-> male_CPS1_spec
      
    #Compute cpue by 1mm size bin for shell/clutch categories  
      male_CPS1_spec %>%
        #create 1mm size bin
        mutate(SIZE_BIN = floor(LENGTH),
               CATCH_PER_HOUR = SAMPLING_FACTOR/SOAK_TIME) %>% 
        dplyr::group_by(YEAR, VESSEL, SPN, POT_ID, BUOY, LAT_DD, LON_DD, SIZE_BIN, 
                        SHELL_TEXT) %>%
        dplyr::reframe(COUNT = sum(SAMPLING_FACTOR),
                       CATCH_PER_HOUR = sum(CATCH_PER_HOUR)) -> positive_male_catch
      
      
    #Join to zero catch stations
      positive_male_catch %>%
        dplyr::right_join(expand_grid(SIZE_BIN = 1:max(positive_male_catch$SIZE_BIN),
                                      potlifts) %>%
                            expand_grid(male_CPS1_spec %>%
                                          distinct(SHELL_TEXT))) %>%
        replace_na(list(COUNT = 0, CATCH_PER_HOUR = 0)) %>%
        group_by(YEAR, SIZE_BIN, SHELL_TEXT) %>%
        reframe(COUNT = sum(COUNT)) %>%
        mutate(YEAR = 2023) -> CPS1_male_1mm
      
    #Specify plot labels
      labs <- paste("CPS1", 2023)
      
      names(labs) <- 2023
      
      
    # Plot shell condition
      CPS1_male_1mm %>%
        dplyr::filter(COUNT>0) %>% 
        mutate(SHELL_TEXT  = factor(SHELL_TEXT , levels = c('Very Old', 'Old', 'New Hard', 'Soft Molting'))) %>%
        ggplot(aes(x=SIZE_BIN, y=COUNT, fill=SHELL_TEXT))+
        geom_bar(position="stack", stat="identity") +
        theme_bw() +
        #xlim(0,plot_lookup$XMax[plot_lookup$stock == stock])+
        scale_x_continuous(expand = c(0.01,0.01), limits = c(16, 192))+ #same as CPS1/NMFS
        ylab("Abundance") +
        xlab("Carapace length (mm) ") +
        scale_fill_manual(name = "Shell Condition", values = c("#018571", "#80cdc1", "#dfc27d", "#a6611a"),
                          labels = c("Very Old", "Old", "New Hard", "Soft Molting"), drop = FALSE)+
        facet_manual(~YEAR, "A", scales = "free_y", labeller = labeller(YEAR = labs)) +
        theme(legend.position = "bottom", 
              legend.background=element_rect(fill="white", color="black"))-> CPS1_maleshell_plot
      
    # Combine NMFS/CPS1 plots, save 
      # Shell condition
      ggarrange(ggarrange(NMFS_male_out$shell_plot+ rremove("xlab"),
                          legend = "none"),
                CPS1_maleshell_plot, heights = c(1,1), nrow= 2) 
      
      ggsave("./Figures/CPS1 Report/CPS1_maleshell_panels.png",
             height=8.5, width=8.5, units="in")
      
# REQUEST 9: BYCATCH AND TEMPERATURE MAPS -------------------------------------------------------------------------------------
  # WITHOUT TEMP
    # Make map labels
      bc_labs <- data.frame(labs = c("Pacific Cod", "Yellowfin Sole"),
                            name = c("PacificCod", "YellowfinSole"))
      
    #Plot YFS
      ggplot() +
        geom_sf(data = st_transform(map_layers$bathymetry, map.crs), color=alpha("grey70")) +
        #geom_sf(data = st_as_sf(BB_strata), fill = NA, mapping = aes(color = "black"), linewidth = 1) +
        geom_sf(data = st_as_sf(RKCSA_sub), mapping = aes(color = "red"), fill = NA, alpha= 0.9, linewidth = 1) +
        geom_sf(data = st_as_sf(RKCSA), fill = NA,  color = "red", alpha =0.5, linewidth = 1) +
        geom_sf(data = st_as_sf(CPS1_bound), fill = NA, aes(color = "black"), linewidth = 1)+
        geom_sf(data = st_transform(map_layers$akland, map.crs), fill = "grey80") +
        geom_sf(data = bycatch,
                mapping = aes(size=YellowfinSole, fill = YellowfinSole, shape = YellowfinSole == 0), 
                alpha = 0.5, colour = "black")+
        scale_shape_manual(values = c('TRUE' = 4, 'FALSE' = 21), guide = "none") +
        scale_color_manual(values = c("black", "red"), 
                           labels = c("CPS1 survey boundary", "Red King Crab Savings Area"),
                           name = "") +
        scale_size_continuous(range = c(2, 10))+ 
        scale_fill_gradientn(colors = c("gray", rev(pal[5:length(pal)]))) +
        scale_x_continuous(breaks = c(-165, -160), labels = paste0(c(165, 160), "°W"))+
        scale_y_continuous(breaks = c(56, 58), labels = paste0(c(56, 58), "°N"))+
        labs(title = "2023 BBRKC Collaborative Pot Sampling", subtitle = paste(filter(bc_labs, name == "YellowfinSole")$lab))+
        guides(size = guide_legend(title.position = "top", title = "COUNT", nrow = 2, override.aes = list(shape = c(4, rep(21, 4)))),
               fill = guide_legend(title = "COUNT"),
               color = guide_legend(nrow = 2))+
        coord_sf(xlim = plot.boundary$x,
                 ylim = plot.boundary$y) +
        geom_sf_text(sf::st_as_sf(data.frame(lab= c("50m", "100m"), 
                                             x = c(-161.5, -165), y = c(58.3, 56.1)),
                                  coords = c(x = "x", y = "y"), crs = sf::st_crs(4326)) %>%
                       sf::st_transform(crs = map.crs),
                     mapping = aes(label = lab))+
        theme_bw() +
        theme(axis.title = element_blank(),
              axis.text = element_text(size = 10),
              legend.text = element_text(size = 10),
              legend.title = element_text(size = 10),
              legend.position = "bottom",
              legend.direction = "horizontal",
              plot.title = element_text(face = "bold", size = 15),
              plot.subtitle = element_text(size = 12)) -> yfs_map
      
    # Plot Pcod
      ggplot() +
        geom_sf(data = st_transform(map_layers$bathymetry, map.crs), color=alpha("grey70")) +
        #geom_sf(data = st_as_sf(BB_strata), fill = NA, mapping = aes(color = "black"), linewidth = 1) +
        geom_sf(data = st_as_sf(RKCSA_sub), mapping = aes(color = "red"), fill = NA, alpha= 0.9, linewidth = 1) +
        geom_sf(data = st_as_sf(RKCSA), fill = NA,  color = "red", alpha =0.5, linewidth = 1) +
        geom_sf(data = st_as_sf(CPS1_bound), fill = NA, aes(color = "black"), linewidth = 1)+
        geom_sf(data = st_transform(map_layers$akland, map.crs), fill = "grey80") +
        geom_sf(data = bycatch,
                mapping = aes(size=PacificCod, fill = PacificCod, shape = PacificCod == 0), 
                alpha = 0.5, colour = "black")+
        scale_shape_manual(values = c('TRUE' = 4, 'FALSE' = 21), guide = "none") +
        scale_color_manual(values = c("black", "red"), 
                           labels = c("CPS1 survey boundary", "Red King Crab Savings Area"),
                           name = "") +
        scale_size_continuous(range = c(2, 10))+ 
        scale_fill_gradientn(colors = c("gray", rev(pal[5:length(pal)]))) +
        scale_x_continuous(breaks = c(-165, -160), labels = paste0(c(165, 160), "°W"))+
        scale_y_continuous(breaks = c(56, 58), labels = paste0(c(56, 58), "°N"))+
        labs(title = "2023 BBRKC Collaborative Pot Sampling", subtitle = paste(filter(bc_labs, name == "PacificCod")$lab))+
        guides(size = guide_legend(title.position = "top",title = "COUNT", nrow = 2, override.aes = list(shape = c(4, rep(21, 3)))),
               fill = guide_legend(title = "COUNT"),
               color = guide_legend(nrow = 2))+
        coord_sf(xlim = plot.boundary$x,
                 ylim = plot.boundary$y) +
        geom_sf_text(sf::st_as_sf(data.frame(lab= c("50m", "100m"), 
                                             x = c(-161.5, -165), y = c(58.3, 56.1)),
                                  coords = c(x = "x", y = "y"), crs = sf::st_crs(4326)) %>%
                       sf::st_transform(crs = map.crs),
                     mapping = aes(label = lab))+
        theme_bw() +
        theme(axis.title = element_blank(),
              axis.text = element_text(size = 10),
              legend.text = element_text(size = 10),
              legend.title = element_text(size = 10),
              legend.position = "bottom",
              legend.direction = "horizontal",
              plot.title = element_text(face = "bold", size = 15),
              plot.subtitle = element_text(size = 12)) -> pcod_map
      
      # Plot Tanner 
      ggplot() +
        geom_sf(data = st_transform(map_layers$bathymetry, map.crs), color=alpha("grey70")) +
        #geom_sf(data = st_as_sf(BB_strata), fill = NA, mapping = aes(color = "black"), linewidth = 1) +
        geom_sf(data = st_as_sf(RKCSA_sub), mapping = aes(color = "red"), fill = NA, alpha= 0.9, linewidth = 1) +
        geom_sf(data = st_as_sf(RKCSA), fill = NA,  color = "red", alpha =0.5, linewidth = 1) +
        geom_sf(data = st_as_sf(CPS1_bound), fill = NA, aes(color = "black"), linewidth = 1)+
        geom_sf(data = st_transform(map_layers$akland, map.crs), fill = "grey80") +
        geom_sf(data = bycatch,
                mapping = aes(size= (MaleTanner + FemaleTanner), fill = (MaleTanner+FemaleTanner), 
                              shape = (MaleTanner+FemaleTanner) == 0), 
                alpha = 0.5, colour = "black")+
        scale_shape_manual(values = c('TRUE' = 4, 'FALSE' = 21), guide = "none") +
        scale_color_manual(values = c("black", "red"), 
                           labels = c("CPS1 survey boundary", "Red King Crab Savings Area"),
                           name = "") +
        scale_size_continuous(range = c(2, 10))+ 
        scale_fill_gradientn(colors = c("gray", rev(pal[5:length(pal)]))) +
        scale_x_continuous(breaks = c(-165, -160), labels = paste0(c(165, 160), "°W"))+
        scale_y_continuous(breaks = c(56, 58), labels = paste0(c(56, 58), "°N"))+
        labs(title = "2023 BBRKC Collaborative Pot Sampling", subtitle = "Tanner crab")+
        guides(size = guide_legend(title.position = "top",title = "COUNT", nrow = 2, override.aes = list(shape = c(4, rep(21, 3)))),
               fill = guide_legend(title = "COUNT"),
               color = guide_legend(nrow = 2))+
        coord_sf(xlim = plot.boundary$x,
                 ylim = plot.boundary$y) +
        geom_sf_text(sf::st_as_sf(data.frame(lab= c("50m", "100m"), 
                                             x = c(-161.5, -165), y = c(58.3, 56.1)),
                                  coords = c(x = "x", y = "y"), crs = sf::st_crs(4326)) %>%
                       sf::st_transform(crs = map.crs),
                     mapping = aes(label = lab))+
        theme_bw() +
        theme(axis.title = element_blank(),
              axis.text = element_text(size = 10),
              legend.text = element_text(size = 10),
              legend.title = element_text(size = 10),
              legend.position = "bottom",
              legend.direction = "horizontal",
              plot.title = element_text(face = "bold", size = 15),
              plot.subtitle = element_text(size = 12)) -> tanner_map
      
      # Plot Snow
      ggplot() +
        geom_sf(data = st_transform(map_layers$bathymetry, map.crs), color=alpha("grey70")) +
        #geom_sf(data = st_as_sf(BB_strata), fill = NA, mapping = aes(color = "black"), linewidth = 1) +
        geom_sf(data = st_as_sf(RKCSA_sub), mapping = aes(color = "red"), fill = NA, alpha= 0.9, linewidth = 1) +
        geom_sf(data = st_as_sf(RKCSA), fill = NA,  color = "red", alpha =0.5, linewidth = 1) +
        geom_sf(data = st_as_sf(CPS1_bound), fill = NA, aes(color = "black"), linewidth = 1)+
        geom_sf(data = st_transform(map_layers$akland, map.crs), fill = "grey80") +
        geom_sf(data = bycatch,
                mapping = aes(size= (MaleSnow + FemaleSnow), fill = (MaleSnow+FemaleSnow), 
                              shape = (MaleSnow+FemaleSnow) == 0), 
                alpha = 0.5, colour = "black")+
        scale_shape_manual(values = c('TRUE' = 4, 'FALSE' = 21), guide = "none") +
        scale_color_manual(values = c("black", "red"), 
                           labels = c("CPS1 survey boundary", "Red King Crab Savings Area"),
                           name = "") +
        scale_size_continuous(range = c(2, 10))+ 
        scale_fill_gradientn(colors = c("gray", rev(pal[5:length(pal)]))) +
        scale_x_continuous(breaks = c(-165, -160), labels = paste0(c(165, 160), "°W"))+
        scale_y_continuous(breaks = c(56, 58), labels = paste0(c(56, 58), "°N"))+
        labs(title = "2023 BBRKC Collaborative Pot Sampling", subtitle = "Snow crab")+
        guides(size = guide_legend(title.position = "top",title = "COUNT", nrow = 2, override.aes = list(shape = c(4, rep(21, 6)))),
               fill = guide_legend(title = "COUNT"),
               color = guide_legend(nrow = 2))+
        coord_sf(xlim = plot.boundary$x,
                 ylim = plot.boundary$y) +
        geom_sf_text(sf::st_as_sf(data.frame(lab= c("50m", "100m"), 
                                             x = c(-161.5, -165), y = c(58.3, 56.1)),
                                  coords = c(x = "x", y = "y"), crs = sf::st_crs(4326)) %>%
                       sf::st_transform(crs = map.crs),
                     mapping = aes(label = lab))+
        theme_bw() +
        theme(axis.title = element_blank(),
              axis.text = element_text(size = 10),
              legend.text = element_text(size = 10),
              legend.title = element_text(size = 10),
              legend.position = "bottom",
              legend.direction = "horizontal",
              plot.title = element_text(face = "bold", size = 15),
              plot.subtitle = element_text(size = 12)) -> snow_map
      
      # Save figures
      ggsave(plot = yfs_map, "./Figures/CPS1 Report/yfs_map.png", height=7, width=10, units="in")
      ggsave(plot = pcod_map, "./Figures/CPS1 Report/pcod_map.png", height=7, width=10, units="in")
      ggsave(plot = tanner_map, "./Figures/CPS1 Report/tanner_map.png", height=7, width=10, units="in")
      ggsave(plot = snow_map, "./Figures/CPS1 Report/snow_map.png", height=7, width=10, units="in")
  
 # WITH TEMPERATURE
      
      #Plot YFS
      ggplot() +
        ggplot2::geom_tile(data = temp_df, 
                           aes(x = x, 
                               y = y,
                               fill = cut(temperature, 
                                          breaks = temp_breaks))) +
        geom_sf(data = st_transform(map_layers$bathymetry, map.crs), color=alpha("grey70")) +
        #geom_sf(data = st_as_sf(BB_strata), fill = NA, mapping = aes(color = "black"), linewidth = 1) +
        #geom_sf(data = st_as_sf(RKCSA_sub), mapping = aes(color = "red"), fill = NA, alpha= 0.9, linewidth = 1) +
        #geom_sf(data = st_as_sf(RKCSA), fill = NA,  color = "red", alpha =0.5, linewidth = 1) +
        geom_sf(data = st_as_sf(CPS1_bound), fill = NA, color = "black", linewidth = 1)+
        geom_sf(data = st_transform(map_layers$akland, map.crs), fill = "grey80") +
        geom_sf(data = bycatch,
                mapping = aes(size=YellowfinSole, shape = YellowfinSole == 0), 
                alpha = 0.5, colour = "black", fill = "black")+
        scale_shape_manual(values = c('TRUE' = 4, 'FALSE' = 21), guide = "none") +
        #scale_color_manual(values = c("black", "red"), 
                           #labels = c("CPS1 survey boundary", "Red King Crab Savings Area"),
                           #name = "") +
        scale_size_continuous(range = c(2, 10))+ 
        scale_x_continuous(breaks = c(-165, -160), labels = paste0(c(165, 160), "°W"))+
        scale_y_continuous(breaks = c(56, 58), labels = paste0(c(56, 58), "°N"))+
        labs(title = "2023 BBRKC Collaborative Pot Sampling", subtitle = paste(filter(bc_labs, name == "YellowfinSole")$lab))+
        guides(size = guide_legend(title.position = "top", title = "COUNT", nrow = 2, override.aes = list(shape = c(4, rep(21, 4)))),
               fill = guide_legend(nrow = 3),
               color = guide_legend(nrow = 2))+
        coord_sf(xlim = plot.boundary$x,
                 ylim = plot.boundary$y) +
        ggplot2::scale_fill_manual(name = "Temperature (°C)", values = viridis_pal(option = viridis_option)(n_temp_breaks),
                                   labels = c(expression(""<=-1), "-0.9-0", "0.1-1", "1.1-2", "2.1-3", "3.1-4",
                                              "4.1-5", "5.1-6", "6.1-7", "7.1-8", ">8.1"), drop = FALSE) +
        geom_sf_text(sf::st_as_sf(data.frame(lab= c("50m", "100m"), 
                                             x = c(-161.5, -165), y = c(58.3, 56.1)),
                                  coords = c(x = "x", y = "y"), crs = sf::st_crs(4326)) %>%
                       sf::st_transform(crs = map.crs),
                     mapping = aes(label = lab))+
        theme_bw() +
        theme(axis.title = element_blank(),
              axis.text = element_text(size = 10),
              legend.text = element_text(size = 10),
              legend.title = element_text(size = 10),
              legend.position = "bottom",
              legend.direction = "horizontal",
              plot.title = element_text(face = "bold", size = 15),
              plot.subtitle = element_text(size = 12)) -> yfs.temp_map
      
      #Plot Pcod
      ggplot() +
        ggplot2::geom_tile(data = temp_df, 
                           aes(x = x, 
                               y = y,
                               fill = cut(temperature, 
                                          breaks = temp_breaks))) +
        geom_sf(data = st_transform(map_layers$bathymetry, map.crs), color=alpha("grey70")) +
        #geom_sf(data = st_as_sf(BB_strata), fill = NA, mapping = aes(color = "black"), linewidth = 1) +
        #geom_sf(data = st_as_sf(RKCSA_sub), mapping = aes(color = "red"), fill = NA, alpha= 0.9, linewidth = 1) +
        #geom_sf(data = st_as_sf(RKCSA), fill = NA,  color = "red", alpha =0.5, linewidth = 1) +
        geom_sf(data = st_as_sf(CPS1_bound), fill = NA, color = "black", linewidth = 1)+
        geom_sf(data = st_transform(map_layers$akland, map.crs), fill = "grey80") +
        geom_sf(data = bycatch,
                mapping = aes(size=PacificCod, shape = PacificCod == 0), 
                alpha = 0.5, colour = "black", fill = "black")+
        scale_shape_manual(values = c('TRUE' = 4, 'FALSE' = 21), guide = "none") +
        #scale_color_manual(values = c("black", "red"), 
        #labels = c("CPS1 survey boundary", "Red King Crab Savings Area"),
        #name = "") +
        scale_size_continuous(range = c(2, 10))+ 
        scale_x_continuous(breaks = c(-165, -160), labels = paste0(c(165, 160), "°W"))+
        scale_y_continuous(breaks = c(56, 58), labels = paste0(c(56, 58), "°N"))+
        labs(title = "2023 BBRKC Collaborative Pot Sampling", subtitle = paste(filter(bc_labs, name == "PacificCod")$lab))+
        guides(size = guide_legend(title.position = "top", title = "COUNT", nrow = 2, override.aes = list(shape = c(4, rep(21, 3)))),
               fill = guide_legend(nrow = 3),
               color = guide_legend(nrow = 2))+
        coord_sf(xlim = plot.boundary$x,
                 ylim = plot.boundary$y) +
        ggplot2::scale_fill_manual(name = "Temperature (°C)", values = viridis_pal(option = viridis_option)(n_temp_breaks),
                                   labels = c(expression(""<=-1), "-0.9-0", "0.1-1", "1.1-2", "2.1-3", "3.1-4",
                                              "4.1-5", "5.1-6", "6.1-7", "7.1-8", ">8.1"), drop = FALSE) +
        geom_sf_text(sf::st_as_sf(data.frame(lab= c("50m", "100m"), 
                                             x = c(-161.5, -165), y = c(58.3, 56.1)),
                                  coords = c(x = "x", y = "y"), crs = sf::st_crs(4326)) %>%
                       sf::st_transform(crs = map.crs),
                     mapping = aes(label = lab))+
        theme_bw() +
        theme(axis.title = element_blank(),
              axis.text = element_text(size = 10),
              legend.text = element_text(size = 10),
              legend.title = element_text(size = 10),
              legend.position = "bottom",
              legend.direction = "horizontal",
              plot.title = element_text(face = "bold", size = 15),
              plot.subtitle = element_text(size = 12)) -> pcod.temp_map
      
    #Plot snow
      ggplot() +
        ggplot2::geom_tile(data = temp_df, 
                           aes(x = x, 
                               y = y,
                               fill = cut(temperature, 
                                          breaks = temp_breaks))) +
        geom_sf(data = st_transform(map_layers$bathymetry, map.crs), color=alpha("grey70")) +
        #geom_sf(data = st_as_sf(BB_strata), fill = NA, mapping = aes(color = "black"), linewidth = 1) +
        #geom_sf(data = st_as_sf(RKCSA_sub), mapping = aes(color = "red"), fill = NA, alpha= 0.9, linewidth = 1) +
        #geom_sf(data = st_as_sf(RKCSA), fill = NA,  color = "red", alpha =0.5, linewidth = 1) +
        geom_sf(data = st_as_sf(CPS1_bound), fill = NA, color = "black", linewidth = 1)+
        geom_sf(data = st_transform(map_layers$akland, map.crs), fill = "grey80") +
        geom_sf(data = bycatch,
                mapping = aes(size= (MaleSnow + FemaleSnow), 
                              shape = (MaleSnow+FemaleSnow) == 0), 
                alpha = 0.5, colour = "black", fill = "black")+
        scale_shape_manual(values = c('TRUE' = 4, 'FALSE' = 21), guide = "none") +
        #scale_color_manual(values = c("black", "red"), 
        #labels = c("CPS1 survey boundary", "Red King Crab Savings Area"),
        #name = "") +
        scale_size_continuous(range = c(2, 10))+ 
        scale_x_continuous(breaks = c(-165, -160), labels = paste0(c(165, 160), "°W"))+
        scale_y_continuous(breaks = c(56, 58), labels = paste0(c(56, 58), "°N"))+
        labs(title = "2023 BBRKC Collaborative Pot Sampling", subtitle = "Snow crab")+
        guides(size = guide_legend(title.position = "top", title = "COUNT", nrow = 2, override.aes = list(shape = c(4, rep(21, 6)))),
               fill = guide_legend(nrow = 3),
               color = guide_legend(nrow = 2))+
        coord_sf(xlim = plot.boundary$x,
                 ylim = plot.boundary$y) +
        ggplot2::scale_fill_manual(name = "Temperature (°C)", values = viridis_pal(option = viridis_option)(n_temp_breaks),
                                   labels = c(expression(""<=-1), "-0.9-0", "0.1-1", "1.1-2", "2.1-3", "3.1-4",
                                              "4.1-5", "5.1-6", "6.1-7", "7.1-8", ">8.1"), drop = FALSE) +
        geom_sf_text(sf::st_as_sf(data.frame(lab= c("50m", "100m"), 
                                             x = c(-161.5, -165), y = c(58.3, 56.1)),
                                  coords = c(x = "x", y = "y"), crs = sf::st_crs(4326)) %>%
                       sf::st_transform(crs = map.crs),
                     mapping = aes(label = lab))+
        theme_bw() +
        theme(axis.title = element_blank(),
              axis.text = element_text(size = 10),
              legend.text = element_text(size = 10),
              legend.title = element_text(size = 10),
              legend.position = "bottom",
              legend.direction = "horizontal",
              plot.title = element_text(face = "bold", size = 15),
              plot.subtitle = element_text(size = 12)) -> snow.temp_map
      
    #Plot Tanner
      ggplot() +
        ggplot2::geom_tile(data = temp_df, 
                           aes(x = x, 
                               y = y,
                               fill = cut(temperature, 
                                          breaks = temp_breaks))) +
        geom_sf(data = st_transform(map_layers$bathymetry, map.crs), color=alpha("grey70")) +
        #geom_sf(data = st_as_sf(BB_strata), fill = NA, mapping = aes(color = "black"), linewidth = 1) +
        #geom_sf(data = st_as_sf(RKCSA_sub), mapping = aes(color = "red"), fill = NA, alpha= 0.9, linewidth = 1) +
        #geom_sf(data = st_as_sf(RKCSA), fill = NA,  color = "red", alpha =0.5, linewidth = 1) +
        geom_sf(data = st_as_sf(CPS1_bound), fill = NA, color = "black", linewidth = 1)+
        geom_sf(data = st_transform(map_layers$akland, map.crs), fill = "grey80") +
        geom_sf(data = bycatch,
                mapping = aes(size= (MaleTanner + FemaleTanner), 
                              shape = (MaleTanner+FemaleTanner) == 0), 
                alpha = 0.5, colour = "black", fill = "black")+
        scale_shape_manual(values = c('TRUE' = 4, 'FALSE' = 21), guide = "none") +
        #scale_color_manual(values = c("black", "red"), 
        #labels = c("CPS1 survey boundary", "Red King Crab Savings Area"),
        #name = "") +
        scale_size_continuous(range = c(2, 10))+ 
        scale_x_continuous(breaks = c(-165, -160), labels = paste0(c(165, 160), "°W"))+
        scale_y_continuous(breaks = c(56, 58), labels = paste0(c(56, 58), "°N"))+
        labs(title = "2023 BBRKC Collaborative Pot Sampling", subtitle = "Tanner crab")+
        guides(size = guide_legend(title.position = "top", title = "COUNT", nrow = 2, override.aes = list(shape = c(4, rep(21, 3)))),
               fill = guide_legend(nrow = 3),
               color = guide_legend(nrow = 2))+
        coord_sf(xlim = plot.boundary$x,
                 ylim = plot.boundary$y) +
        ggplot2::scale_fill_manual(name = "Temperature (°C)", values = viridis_pal(option = viridis_option)(n_temp_breaks),
                                   labels = c(expression(""<=-1), "-0.9-0", "0.1-1", "1.1-2", "2.1-3", "3.1-4",
                                              "4.1-5", "5.1-6", "6.1-7", "7.1-8", ">8.1"), drop = FALSE) +
        geom_sf_text(sf::st_as_sf(data.frame(lab= c("50m", "100m"), 
                                             x = c(-161.5, -165), y = c(58.3, 56.1)),
                                  coords = c(x = "x", y = "y"), crs = sf::st_crs(4326)) %>%
                       sf::st_transform(crs = map.crs),
                     mapping = aes(label = lab))+
        theme_bw() +
        theme(axis.title = element_blank(),
              axis.text = element_text(size = 10),
              legend.text = element_text(size = 10),
              legend.title = element_text(size = 10),
              legend.position = "bottom",
              legend.direction = "horizontal",
              plot.title = element_text(face = "bold", size = 15),
              plot.subtitle = element_text(size = 12)) -> tanner.temp_map
      
      
      # Save figures
      ggsave(plot = yfs.temp_map, "./Figures/CPS1 Report/yfs.temp_map.png", height=7, width=10, units="in")
      ggsave(plot = pcod.temp_map, "./Figures/CPS1 Report/pcod.temp_map.png", height=7, width=10, units="in")
      ggsave(plot = tanner.temp_map, "./Figures/CPS1 Report/tanner.temp_map.png", height=7, width=10, units="in")
      ggsave(plot = snow.temp_map, "./Figures/CPS1 Report/snow.temp_map.png", height=7, width=10, units="in")
      
# REQUEST 10: # CRAB CAUGHT IN CLOSURE AREAS ----------------------------------------------------------------------------------
 # Pull out column letters and column numbers into separate columns  
 pot_cpue %>%
   mutate(COL = substr(POT_ID, 1, 1), COL_NUM = substr(POT_ID, 2, 3)) -> pot_cpue_cols
 
 # East of RKCSA
 pot_cpue_cols %>%
   filter((COL == "F" & COL_NUM %in% 1:48)|
            (COL == "G" & COL_NUM %in% 1:56)|
            (COL == "H" & COL_NUM %in% 1:54)|
            (COL == "I" & COL_NUM %in% 1:54)|
            (COL == "J" & COL_NUM %in% 1:42)|
            (COL == "K" & COL_NUM %in% 1:37)) %>%
   group_by(MAT_SEX) %>%
   reframe(SUM_COUNT = sum(COUNT)) -> RKCSA_E
 
 # All RKCSA
 pot_cpue_cols %>%
   filter((COL == "B" & COL_NUM %in% 26:55)|
            (COL == "C" & COL_NUM %in% 26:55)|
            (COL == "D" & COL_NUM %in% 26:55)|
            (COL == "E" & COL_NUM %in% 26:55)|
            (COL == "F" & COL_NUM %in% 49:55)) %>%
   group_by(MAT_SEX) %>%
   reframe(SUM_COUNT = sum(COUNT)) -> RKCSA_all
 
 # RKCSA subarea
 pot_cpue_cols %>%
   filter((COL == "B" & COL_NUM %in% 51:55)|
            (COL == "C" & COL_NUM %in% 51:55)|
            (COL == "D" & COL_NUM %in% 51:55)|
            (COL == "E" & COL_NUM %in% 51:55)|
            (COL == "F" & COL_NUM %in% 51:55)) %>%
   group_by(MAT_SEX) %>%
   reframe(SUM_COUNT = sum(COUNT)) -> RKCSA_subarea
 
 # South of RKCSA
 pot_cpue_cols %>%
   filter((COL == "A" & COL_NUM %in% 56:85)|
            (COL == "B" & COL_NUM %in% 56:83)|
            (COL == "C" & COL_NUM %in% 56:81)|
            (COL == "D" & COL_NUM %in% 56:74)|
            (COL == "E" & COL_NUM %in% 56:67)|
            (COL == "F" & COL_NUM %in% 56:60)) %>%
   group_by(MAT_SEX) %>%
   reframe(SUM_COUNT = sum(COUNT)) -> RKCSA_S
 
 # West of RKCSA
 pot_cpue_cols %>%
   filter((COL == "A" & COL_NUM %in% 26:55)) %>%
   group_by(MAT_SEX) %>%
   reframe(SUM_COUNT = sum(COUNT)) -> RKCSA_W
 
 # North of RKCSA
 pot_cpue_cols %>%
   filter((COL == "A" & COL_NUM %in% 1:25)|
            (COL == "B" & COL_NUM %in% 1:25)|
            (COL == "C" & COL_NUM %in% 1:25)|
            (COL == "D" & COL_NUM %in% 1:25)|
            (COL == "E" & COL_NUM %in% 1:25)) %>%
   group_by(MAT_SEX) %>%
   reframe(SUM_COUNT = sum(COUNT)) -> RKCSA_N
 