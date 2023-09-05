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
         paste0("./Figures/CPS1 Report/Abundance_timeseries", stock, ".png"),
         height=7.5, width=10, units="in")

# REQUEST 2: MAP OF TRAWL CLOSURE AREAS W/CPS1 SURVEY AREA --------------------------------------------------------------------------------
  # Load boundaries
  st_read(survey_gdb,layer="Area516") ->  area516
  
  st_read(survey_gdb,layer="BycatchZone1") ->  zone1
  
  st_read(survey_gdb,layer="NBBTCA") ->  nbbtca
  st_read(survey_gdb,layer="TogiakTrawlArea") ->  togtrawl
  
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
   geom_sf(data = st_as_sf(BB_strata), fill = NA, aes(color = "black"), linewidth = 1) +
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
   scale_fill_manual(values = c(alpha("firebrick3", 0.5), alpha("firebrick2", 0.25), alpha("gold2", 0.5)),
                     labels = c("Red King Crab Savings Area", "Red King Crab Savings Subarea", "Togiak Trawl Area"),
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
 sta_labs <- data.frame(lab = c("A", "B"),
                        x = c(-164.2729, -163.8046),
                                  y = c(57.99736, 58.01606)) %>%
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
 
# REQUEST 9: # CRAB CAUGHT IN CLOSURE AREAS ----------------------------------------------------------------------------------
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
 