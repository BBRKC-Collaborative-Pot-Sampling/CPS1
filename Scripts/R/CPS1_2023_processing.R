# PURPOSE -------------------------------------------------------------------------------------------------------------------
  # 1) To automate processing specimen tables and catch summaries from Collaborative Pot Sampling 1 (CPS1) 2023 for 
  #    BBRKC
  # 2) To run error checks on processed specimen and catch summaries
  # 3) To calculate and map cpue by pot and mat/sex category for BBRKC

# INSTALL PACKAGES ----------------------------------------------------------------------------------------------------------
  #install.packages(c("tidyverse", "gsubfn", "terra", "rgdal", "colorRamps", "sf", "viridis", "grid", "shadowtext",
#                     "ggnewscale"))

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

# LOAD DATA -----------------------------------------------------------------------------------------------------------------
  
  # Load summary catch and specimen tables
      catch <- list.files("./Data/Catch - FTP/") %>%
        purrr::map_df(~read.csv(paste0("./Data/Catch - FTP/", .x))) 
      
      specimen <- list.files("./Data/Specimen - FTP/") %>%
        purrr::map_df(~read.csv(paste0("./Data/Specimen - FTP/", .x))%>% mutate(STATION = paste0("X", STATION))) 
      
  # Load raw data for processing below   
      raw_sample <- list.files("./Data/Raw Data - FTP/", pattern = "_SAMPLE_0") %>% # RECORDS of SAMPLE INFO
        purrr::map_df(~read.csv(paste0("./Data/Raw Data - FTP/", .x))) #E.G. SEX, SPECIES
      
      raw_sample_values <- list.files("./Data/Raw Data - FTP/", pattern = "_SAMPLE_VALUES") %>% #RECORDS OF # TOSSED
        purrr::map_df(~read.csv(paste0("./Data/Raw Data - FTP/", .x))) %>%
        mutate(TOSSED = ifelse(is.na(COUNT) == FALSE, COUNT,0)) %>%
        group_by(HAUL_ID, CATCH_SAMPLE_ID) %>%
        dplyr::reframe(TOSSED = sum(TOSSED)) 
      
      raw_specimen <- list.files("./Data/Raw Data - FTP/", pattern = "_SPECIMEN_0") %>% 
        purrr::map_df(~read.csv(paste0("./Data/Raw Data - FTP/", .x))) %>%
        dplyr::select(HAUL_ID, SPECIMEN_ID, CATCH_SAMPLE_ID, SPECIES_CODE)
      
      raw_specimen_bio<- list.files("./Data/Raw Data - FTP/", pattern = "_SPECIMEN_BIOMETRICS") %>% 
        purrr::map_df(~read.csv(paste0("./Data/Raw Data - FTP/", .x))) %>%
        right_join(., catch %>% select(HAUL, HAUL_ID, RECORDING_DEVICE))
      
  # Read in potlifts and tagging data
      SB_potlifts <- read.csv("./Data/SB_POTLIFTS .csv") #Summer Bay
      
      SS_potlifts <- read.csv("./Data/SS_POTLIFTS.csv") %>%
                      dplyr::filter(SPN %in% 1:339) %>% #Hardcoding SS data to removing hotspot tagging sites and Gear code 42
                      dplyr::mutate(SPN = as.integer(SPN), GEAR_CODE = ifelse(GEAR_CODE == 42, 42, "")) %>%
                      dplyr::select(!c(ADFG.Logger, Cory.temp, X.1, X.2, X.3, X)) # Silver Spray
      
      potlifts <- rbind(SB_potlifts, SS_potlifts) %>%
        filter(DATE_HAUL != "", is.na(VESSEL) == "FALSE" & is.na(GEAR_CODE) == TRUE | GEAR_CODE == "") %>% 
        mutate(BUOY = paste0("X", BUOY)) 

    
      tagging <- list.files("./Data/", pattern = "TAGGING", ignore.case = TRUE) %>% #MAY NEED TO CHANGE
        purrr::map_df(~read.csv(paste0("./Data/", .x)) %>% select(!CAPTURE_SPN)) %>%
        filter(is.na(LON_MIN) == FALSE) 
      
      #hardcoding common release point for these hotspots for mapping purposes;
      #actual release points are slightly different and can be found in tagging file
      
      tagging$LAT_MIN[which(tagging$DISK %in% 957:961)] <- 3.18
      tagging$LAT_MIN[which(tagging$DISK %in% c(963:966, 956))] <- 53.29
      
      tagging$LON_MIN[which(tagging$DISK %in% 957:961)] <- 26.9
      tagging$LON_MIN[which(tagging$DISK %in% c(963:966, 956))] <- 20.42
      
  # Read in spatial layers for mapping purposes 
      # Set crs
      map.crs <- "EPSG:3338"
      
      # Read in Bristol Bay management shapefile
      survey_gdb<- "./Data/SAP_layers.gdb"
      
      st_read(survey_gdb,layer="BristolBaySurveyStrata") %>%
        st_transform(map.crs) %>%
        vect() -> BB_strata
      
      # Read in RKCSA + subarea shapefiles
      st_read("./Data/RKCSA layers/RKCSA_sub.shp") %>%
        st_transform(map.crs) %>%
        vect() -> RKCSA_sub
      
      st_read("./Data/RKCSA layers/RKCSA.shp") %>%
        st_transform(map.crs) %>%
        vect() -> RKCSA
      
      st_read("./Data/CPS1_survey_boundary.shp") %>%
        st_transform(map.crs) %>%
        vect() -> CPS1_bound
      
      st_read("./Data/CPS1_project_depth_contour.shp") %>%
          st_transform(map.crs)  -> CPS1_bathy
      
# PROCESS DATA ----------------------------------------------------------------------------------------------------------------
    
  # Calculate soak time and lat/lon in degrees decimal for all potlifts, omit bad or gear testing potlifts based on gear code
    potlifts %>%
      dplyr::mutate(DATETIME_SET = as.POSIXct(paste(DATE_SET, TIME_SET), format = "%m/%d/%Y %H:%M"),
                    DATETIME_HAUL = as.POSIXct(paste(DATE_HAUL, TIME_HAUL), format = "%m/%d/%Y %H:%M"),
                    SOAK_TIME = as.numeric(difftime(DATETIME_HAUL, DATETIME_SET, units = "hours")),
                    LAT_DD = LAT_DEG + LAT_MIN/60,
                    LON_DD = (LON_DEG + LON_MIN/60)*-1) %>%
      dplyr::select(!c(DATETIME_SET, DATETIME_HAUL)) -> potlifts
    
  # Calculate lat/lon in degrees decimal for tagging release points and # of tags released per station
    tagging %>%
      dplyr::mutate(LAT_DD = LAT_DEG + LAT_MIN/60,
                    LON_DD = (LON_DEG + LON_MIN/60)*-1) %>%
      dplyr::filter(is.na(VESSEL) == "FALSE") %>%
      dplyr::group_by(VESSEL, LON_DD, LAT_DD) %>%
      dplyr::reframe(LAT_DD = LAT_DD, LON_DD = LON_DD, N = n()) %>%
      dplyr::distinct()-> tagging
    
    
  # Join raw_sample_values and raw_sample to get # tossed per haul, sex, and catch sample id
      samples <- right_join(raw_sample, raw_sample_values) %>%
        right_join(., catch %>% select(HAUL, HAUL_ID, RECORDING_DEVICE)) %>%
        dplyr::select(HAUL, HAUL_ID, CATCH_SAMPLE_ID, SPECIES_CODE, SPECIES_NAME, SEX, TOSSED, RECORDING_DEVICE)
    
  # Expand specimen biometric table, join to raw_specimen table to get catch sample ID, join with samples file to get 
  # number tossed
      raw_specimen_bio %>%
        dplyr::select(HAUL, HAUL_ID, SPECIMEN_ID, BIOMETRIC_NAME, VALUE, RECORDING_DEVICE) %>%
        pivot_wider(., id_cols = c(HAUL, HAUL_ID, SPECIMEN_ID, RECORDING_DEVICE), 
                    names_from = "BIOMETRIC_NAME", values_from = "VALUE") %>%
        dplyr::rename(SHELL_CONDITION = CRAB_SHELL_CONDITION, EGG_COLOR = CRAB_EGG_COLOR,
                      EGG_CONDITION = CRAB_EGG_CONDITION, CLUTCH_SIZE = CRAB_EGG_CLUTCH_SIZE,
                      LENGTH = CARAPACE_LENGTH) %>%
        #right_join(., raw_specimen) %>%
        right_join(samples, ., by = c("HAUL", "HAUL_ID", "SEX", "RECORDING_DEVICE"), 
                   multiple = "all") -> specimen_sum
      
  # Calculate sampling factor from specimen summary table, join back with specimen_sum file to 
  # get specimen information, join with catch file to get vessel and pot #s, join with potlifts
  # file to get lat/lon, set/haul date and time for each pot (with positive catch)
      specimen_sum %>%
        dplyr::group_by(HAUL, HAUL_ID, CATCH_SAMPLE_ID, SEX, RECORDING_DEVICE) %>%
        dplyr::reframe(KEPT = n(),
                  TOSSED = TOSSED,
                  SAMPLING_FACTOR = (KEPT + TOSSED)/KEPT) %>%
        distinct() %>%
        dplyr::right_join(specimen_sum, by = c("HAUL", "HAUL_ID", "CATCH_SAMPLE_ID", "SEX", "TOSSED",
                                               "RECORDING_DEVICE"),
                          multiple = "all") %>%
        dplyr::right_join(catch,., by = c("HAUL", "HAUL_ID", "SPECIES_CODE", "RECORDING_DEVICE"),
                          multiple = "all") %>%
        dplyr::left_join(specimen) %>%
        distinct() %>%
        dplyr::rename(SPN = HAUL) %>%
        dplyr::right_join(potlifts, ., by = c("VESSEL", "SPN"), multiple = "all") %>%
        dplyr::mutate(VESSEL = ifelse(VESSEL == 162, "Summer Bay", "Silver Spray")) %>%
        dplyr::filter(c(is.na(LAT_DD) & is.na(LON_DD) & is.na(SPN)) == FALSE) %>% # bad/gear testing hauls will have NA
        dplyr::select(CRUISE, VESSEL, SPN, POT_ID, BUOY, LAT_DD, LON_DD, DATE_HAUL, TIME_HAUL, SOAK_TIME, DEPTH_F,
                      SPECIES_CODE, SEX, LENGTH, WIDTH, SAMPLING_FACTOR, SHELL_CONDITION, EGG_COLOR, EGG_CONDITION, 
                      CLUTCH_SIZE, WEIGHT, DISEASE_CODE, DISEASE_DORSAL, DISEASE_VENTRAL, DISEASE_LEGS,  
                      CHELA_HEIGHT, MERUS_LENGTH, COMMENTS, NOTES.x) -> specimen_table
      
  # Changing mature barren females to immature for Summer Bay (determined all females
  # coded as mature barren were immature after seeing more immature females in the 
  # last couple of strings in the survey, comparing to mature females, Silver Spray data
  # (the other research vessel), and various other references)
      specimen_table %>%
        dplyr::mutate(CLUTCH_SIZE = ifelse((VESSEL == "Summer Bay" & SEX == 2 & EGG_CONDITION == 0 
                                     & EGG_COLOR == 0), 0, CLUTCH_SIZE)) -> specimen_table
      
  # Process specimen table for Oracle, save
      specimen_table %>%
        dplyr::select(!c(LAT_DD, LON_DD, DATE_HAUL, TIME_HAUL, SOAK_TIME, DEPTH_F, NOTES.x)) %>%
        dplyr::rename(HAUL = POT_ID, STATION = BUOY) %>% # MAY CHANGE BUOY TO ACTUAL STATION
        write.csv("./DataForOracle/Processed_Specimen_Data.csv")
      
  # Process specimen table with all haul data, save
      specimen_table %>% 
        rename(NOTES = NOTES.x) %>%
        write.csv("./Output/CPS1_2023_Processed_Specimen_Data.csv")
      
  # Update catch summary table with new crab #s from sampling factor
      specimen_table %>%
        dplyr::group_by(CRUISE, VESSEL, SPN, POT_ID, SPECIES_CODE) %>%
        dplyr::reframe(NUMBER_CRAB = sum(SAMPLING_FACTOR)) %>%
        dplyr::right_join(catch %>% dplyr::rename(SPN = HAUL, N_ENTRIES = NUMBER_CRAB) %>%
                            dplyr::mutate(VESSEL = ifelse(VESSEL == 162, "Summer Bay", "Silver Spray"))) %>%
        dplyr::select(CRUISE, VESSEL, SPN, POT_ID, SPECIES_CODE, NUMBER_CRAB, N_ENTRIES) %>%
        na.omit() -> catch_summary # bad/gear testing potlifts will have NA for SPN and # crab
      
  # Process catch_summary table for Oracle, save
      catch_summary %>%
        dplyr::select(!N_ENTRIES) %>%
        dplyr::rename(HAUL = SPN) %>%
        write.csv("./DataForOracle/Processed_Catch_Summary.csv")

# CALCULATE AND MAP BBRKC CPUE ---------------------------------------------------------------------------------------------     
  
  # Make non-overlapping maturity/sex and legal/sublegal categories, bind together
      specimen_table %>%
        dplyr::mutate(MAT_SEX = dplyr::case_when((SPECIES_CODE == 69322 & SEX == 1 & LENGTH >= 120) ~ "Mature male",
                                   (SPECIES_CODE == 69322 & SEX == 1 & LENGTH < 120) ~ "Immature male",
                                   (SPECIES_CODE == 69322 & SEX == 2 & CLUTCH_SIZE >= 1) ~ "Mature female",
                                   (SPECIES_CODE == 69322 & SEX == 2 & CLUTCH_SIZE == 0) ~ "Immature female")) -> maturity
      
      specimen_table%>%
        dplyr::mutate(MAT_SEX = dplyr::case_when((SPECIES_CODE == 69322 & SEX == 1 & LENGTH >= 135) ~ "Legal male",
                                   (SPECIES_CODE == 69322 & SEX == 1 & LENGTH < 135) ~ "Sublegal male")) %>%
        dplyr::filter(is.na(MAT_SEX) == "FALSE") -> legal
      
      
      rbind(maturity, legal) -> mat_spec #bind
      
      
  # Calculate COUNT and CATCH_PER_HOUR per pot, CHANGE VESSEL TO ACTUAL NAME
      mat_spec %>%
        dplyr::mutate(CATCH_PER_HOUR = SAMPLING_FACTOR/SOAK_TIME) %>% 
        dplyr::group_by(VESSEL, SPN, POT_ID, BUOY, LAT_DD, LON_DD, MAT_SEX) %>%
        dplyr::reframe(COUNT = sum(SAMPLING_FACTOR),
                         CATCH_PER_HOUR = sum(CATCH_PER_HOUR)) -> positive_pot_cpue
      
    # Change vessel #s to names in potlifts file
      potlifts %>%
        dplyr::mutate(VESSEL = ifelse(VESSEL == 162, "Summer Bay", "Silver Spray")) -> potlifts
    
  # Expand potlifts file to all mat-sex categories and potlifts, join to positive catch file to get zeros 
      mat_sex_combos <- c("Mature male", "Immature male", "Mature female", "Immature female", "Legal male", "Sublegal male")
      
      positive_pot_cpue %>%
        dplyr::right_join(expand_grid(MAT_SEX = mat_sex_combos,
                               potlifts)) %>%
        replace_na(list(COUNT = 0, CPUE = 0)) %>%
        dplyr::select(VESSEL, SPN, POT_ID, BUOY, LAT_DD, LON_DD, DATE_SET, TIME_SET, DATE_HAUL, TIME_HAUL, SOAK_TIME,
                      MAT_SEX, COUNT, CATCH_PER_HOUR)-> pot_cpue
      
  # Save csv
      write.csv(pot_cpue, "./Output/CPS1_2023_potcatch.csv")
      
  # Transform catch and tagging data to correct crs
      pot_cpue %>%
        sf::st_as_sf(coords = c(x = "LON_DD", y = "LAT_DD"), crs = sf::st_crs(4326)) %>%
        sf::st_transform(crs = map.crs)-> pot_cpue_mapdat 

  
  # Remove boat names from tagging data
      tagging %>%
        sf::st_as_sf(coords = c(x = "LON_DD", y = "LAT_DD"), crs = sf::st_crs(4326)) %>%
        sf::st_transform(crs = map.crs) %>%
        dplyr::mutate(VESSEL = ifelse(VESSEL == "162", "Summer Bay", "Silver Spray"),
                      MAT_SEX = "Mature male",
                      shp = case_when((N > 1) ~ "n>1",
                                      (N == 1) ~ "n=1")) -> tagging_mapdat 
      
  #set up plotting features
      map_layers <- readRDS("./Data/akgfmaps_layers.rds")
      
      plot.boundary.untrans <- data.frame(y = c(54.5, 58.5), 
                                          x = c(-164.8, -159)) # plot boundary unprojected
      
      plot.boundary <- plot.boundary.untrans %>%
        sf::st_as_sf(coords = c(x = "x", y = "y"), crs = sf::st_crs(4326)) %>%
        sf::st_transform(crs = map.crs) %>%
        sf::st_coordinates() %>%
        as.data.frame() %>%
        dplyr::rename(x = X, y = Y) # plot boundary projected
      
      
      max.date <- max(pot_cpue_mapdat$DATE_HAUL) # label for most recent pot haul date
      
      # Set up shape mapping
      shapes <- c(0,15) #set shape mapping
      names(shapes) <- c("n=1", "n>1")
      
      # Set up labels
      mat_labs <- c("Mature female", "Immature female", "Mature male (>= 120mm)", "Legal male (>= 135mm)", "Immature male (< 120mm)", "Sublegal male")
      names(mat_labs) <- c("Mature female", "Immature female", "Mature male", "Legal male", "Immature male", "Sublegal male")
      
      mat_labs <- data.frame(lab = c("Mature female", "Immature female", "Mature male (>= 120mm)", 
                                     "Legal male (>= 135mm)", "Immature male (< 120mm)", "Sublegal male"),
                             MAT_SEX = c("Mature female", "Immature female", "Mature male", "Legal male", "Immature male", "Sublegal male"))
      
      # Specify palette
      pal <- viridis::mako(10) # set palette
      
      # Plot without tagging
      mat_sex_combos %>%
        purrr::map(~ggplot() +
                     geom_sf(data = st_transform(map_layers$bathymetry, map.crs), color=alpha("grey70")) +
                     geom_sf(data = st_as_sf(BB_strata), fill = NA, mapping = aes(color = "black"), linewidth = 1) +
                     geom_sf(data = st_as_sf(RKCSA_sub), mapping = aes(color = "red"), fill = NA, alpha= 0.9, linewidth = 1) +
                     geom_sf(data = st_as_sf(RKCSA), fill = NA,  color = "red", alpha =0.5, linewidth = 1) +
                     geom_sf(data = st_transform(map_layers$akland, map.crs), fill = "grey80") +
                     geom_sf(data = filter(pot_cpue_mapdat, MAT_SEX == .x),
                             mapping = aes(size=COUNT, fill = COUNT, shape = COUNT == 0), alpha = 0.5, colour = "black")+
                     scale_shape_manual(values = c('TRUE' = 4, 'FALSE' = 21), guide = "none")+
                     scale_color_manual(values = c("black", "red"), 
                                        labels = c("EBS Summer Survey Boundary", "Red King Crab Savings Area"),
                                        name = "") +
                     scale_size_continuous(range = c(2, 10), limits = c(0, max(pot_cpue_mapdat$COUNT)), 
                                           breaks =seq(0, max(pot_cpue_mapdat$COUNT), by = 25))+ 
                     scale_fill_gradientn(breaks = seq(0, max(pot_cpue_mapdat$COUNT), by = 25),
                                          limits = c(0, max(pot_cpue_mapdat$COUNT)), 
                                          colors = c("gray", rev(pal[5:length(pal)])))+
                     scale_x_continuous(breaks = map_layers$lon.breaks)+
                     scale_y_continuous(breaks = map_layers$lat.breaks)+
                     labs(title = "2023 BBRKC Collaborative Pot Sampling", subtitle = paste(filter(mat_labs, MAT_SEX == .x)$lab))+
                     guides(size = guide_legend(title.position = "top", nrow = 2, override.aes = list(shape = c(4, rep(21, 5)))),
                            fill = guide_legend(),
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
                           plot.subtitle = element_text(size = 12))) -> BBRKC.maps
      
      # Plot with tagging
      mat_sex_combos %>%
        purrr::map(~ggplot() +
                     geom_sf(data = st_transform(map_layers$bathymetry, map.crs), color=alpha("grey70")) +
                     geom_sf(data = st_as_sf(BB_strata), fill = NA, mapping = aes(color = "black"), linewidth = 1) +
                     geom_sf(data = st_as_sf(RKCSA_sub), mapping = aes(color = "red"), fill = NA, alpha= 0.9, linewidth = 1) +
                     geom_sf(data = st_as_sf(RKCSA), fill = NA,  color = "red", alpha =0.5, linewidth = 1) +
                     geom_sf(data = st_transform(map_layers$akland, map.crs), fill = "grey80") +
                     geom_sf(data = filter(pot_cpue_mapdat, MAT_SEX == .x),
                             mapping = aes(size=COUNT, fill = COUNT, shape = COUNT == 0), alpha = 0.5, colour = "black")+
                     geom_sf(data = filter(tagging_mapdat, MAT_SEX == .x),
                             mapping = aes(shape = as.factor(shp)), size= 4, stat="identity", position="identity")+
                     scale_shape_manual(values = c('TRUE' = 4, 'FALSE' = 21, "n=1" = 0, "n>1" = 15),
                                        labels = c("n>1", "n>1", "n>1", "n=1", "n>1"),
                                        breaks = c("n>1", "n>1", "n>1", "n=1", "n>1"))+
                     scale_color_manual(values = c("black", "red"), 
                                        labels = c("EBS Summer Survey Boundary", "Red King Crab Savings Area"),
                                        name = "") +
                     scale_size_continuous(range = c(2, 10), limits = c(0, max(pot_cpue_mapdat$COUNT)))+ 
                     scale_fill_gradientn(limits = c(0, max(pot_cpue_mapdat$COUNT)), 
                                          colors = c("gray", rev(pal[5:length(pal)])))+
                     scale_x_continuous(breaks = map_layers$lon.breaks)+
                     scale_y_continuous(breaks = map_layers$lat.breaks)+
                     labs(title = "2023 BBRKC Collaborative Pot Sampling", subtitle = paste(filter(mat_labs, MAT_SEX == .x)$lab))+
                     guides(size = guide_legend(title.position = "top", nrow = 2, override.aes = list(shape = c(4, rep(21, 5)))), 
                            fill = guide_legend(), shape = guide_legend(title = "TAG RELEASE", title.position = "top", nrow =2),
                            #,
                            #override.aes = list(shape = c(0, 15))),
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
                           plot.subtitle = element_text(size = 12))) -> tagging_map
      
      
      # Save plots
      ggsave(plot = BBRKC.maps[[1]], "./Figures/BBRKC_matmale.png", height=7, width=10, units="in")
      ggsave(plot = tagging_map[[1]], "./Figures/BBRKC_matmaletagging.png", height=7, width=10, units="in")
      ggsave(plot = BBRKC.maps[[2]], "./Figures/BBRKC_immale.png", height=7, width=10, units="in")
      ggsave(plot = BBRKC.maps[[3]], "./Figures/BBRKC_matfem.png", height=7, width=10, units="in")
      ggsave(plot = BBRKC.maps[[4]], "./Figures/BBRKC_imfem.png", height=7, width=10, units="in")
      ggsave(plot = BBRKC.maps[[5]], "./Figures/BBRKC_legalmale.png", height=7, width=10, units="in")
      ggsave(plot = BBRKC.maps[[6]], "./Figures/BBRKC_sublegalmale.png", height=7, width=10, units="in")
 
# ERROR CHECKING -----------------------------------------------------------------------------------------------------
      
  # Write function
    error_chk <- function(specimen_table, catch_summary, potlifts, pot_cpue){
        
        print("CHECKING VESSEL AND CRUISE...")
        
        # 1) Does cruise number match 202201?"
        if(unique(specimen_table$CRUISE == 202201) == FALSE){
          print("ERROR: wrong cruise number")
        }
        
        
        # "2) Does the vessel # match the vessels utilized in the survey?"
        if(FALSE %in% (unique(specimen_table$VESSEL) %in% c("Summer Bay", "Silver Spray")) == TRUE){
          print("ERROR: vessel numbers entered do not match survey vessels")
        }
        
        
        print("Inventory of catch by cruise and vessel")
        specimen_table %>%
          dplyr::group_by(VESSEL, CRUISE) %>%
          dplyr::reframe(COUNT = sum(SAMPLING_FACTOR),
                         N = n()) %>%
          as.data.frame() -> invent1
        print(invent1)
        
        print("CHECKING SPECIES CODES AND SEX...")
        
        # 3) Do species codes match RKC code 69322?"
        if(unique(specimen_table$SPECIES_CODE == 69322) == FALSE){
          print("ERROR: wrong code entered for RKC")
        }
        
        
        # 4) Are sex codes assigned to either 1 or 2?"
        if(unique(specimen_table$SEX %in% c(1:2)) == FALSE){
          print("ERROR: sex code not 1 or 2")
        }
        
        print("Inventory of catch by cruise, vessel, species code, and sex")
        specimen_table %>%
          dplyr::group_by(VESSEL, CRUISE, SPECIES_CODE, SEX) %>%
          dplyr::reframe(COUNT = sum(SAMPLING_FACTOR), 
                         N = n()) %>%
          as.data.frame() -> invent2
        print(invent2)
        
        
        print("CHECKING SHELL CONDITION, EGG COLOR, EGG CONDITION, CLUTCH_SIZE...")
        
        #5) Are egg color codes valid for females?"
        if(FALSE %in% (unique(filter(specimen_table, SEX != 1)$EGG_COLOR %in% c(0, 2:6))) == TRUE){
          print("ERROR: invalid female egg color code (not 0 or 2:6)")
        }
        
        #6) Are egg condition codes valid for females?"
        if(FALSE %in% (unique(filter(specimen_table, SEX != 1)$EGG_CONDITION %in% c(0:4))) == TRUE){
          print("ERROR: invalid female egg condition code (not 0:4)")
        } 
        
        # 7) Are clutch size codes valid for females?
        if(FALSE %in% (unique(filter(specimen_table, SEX != 1)$CLUTCH_SIZE %in% c(0:6))) == TRUE){
          print("ERROR: invalid female clutch size code (not 0:6)")
        }
        
        # 8) Any egg, egg condition, or clutch size codes assigned to males?
        if(TRUE %in% (unique((specimen_table$SEX == 1 & is.na(specimen_table$EGG_COLOR | specimen_table$EGG_CONDITION | 
                                                   specimen_table$CLUTCH_SIZE) == FALSE))) == TRUE){
          print("ERROR: egg, egg condition, or clutch size code assigned to male")
        } 
        
        
        # 9) Any questionable egg condition x shell condition combinations for females?"
        
        # Checking shell condition = 0 and egg condition = 1
        if(TRUE %in% (unique(filter(specimen_table, SEX != 1)$SHELL_CONDITION == 0 & 
                  filter(specimen_table, SEX != 1)$EGG_CONDITION == 1)) == TRUE){
          print("ERROR: female with shell condition = 0 and egg condition = 1")
        } 
        
        # Checking shell condition = 1 and egg condition >1
        if(TRUE %in% (unique(filter(specimen_table, SEX != 1)$SHELL_CONDITION == 1 & 
                  filter(specimen_table, SEX != 1)$EGG_CONDITION > 1)) == TRUE){
          print("ERROR: female with shell condition = 1 and egg condition >1")
        }
        
        # Checking shell condition = 3, 4, or 5 and egg condition = 1
        if(TRUE %in% (unique(filter(specimen_table, SEX != 1)$SHELL_CONDITION %in% c(3:5) & 
                  filter(specimen_table, SEX != 1)$EGG_CONDITION == 1)) == TRUE){
          print("ERROR: female with shell condition 3:5 and egg condition = 1")
        }
        
        # Checking shell condition = 1 and egg condition >=2
        if(TRUE %in% (unique(filter(specimen_table, SEX != 1)$SHELL_CONDITION ==1 & 
                  filter(specimen_table, SEX != 1)$EGG_CONDITION >=2)) == TRUE){
          print("ERROR: female with shell condition = 1 and egg condition >=2")
        }
        
        #10) Any females without egg color, egg condition, or clutch codes?
        if(TRUE %in% (unique(specimen_table$SEX == 2 & (is.na(specimen_table$EGG_COLOR | specimen_table$EGG_CONDITION 
                                                   | specimen_table$CLUTCH_SIZE) == TRUE))) == TRUE){
          print("ERROR: female missing egg color, egg condition, or clutch code")
        }
        
        
        print("Inventory of shell condition")
        specimen_table %>%
          dplyr::group_by(SPECIES_CODE, SEX, SHELL_CONDITION) %>%
          dplyr::reframe(COUNT = sum(SAMPLING_FACTOR),
                         N = n()) %>%
          as.data.frame() -> invent3
        print(invent3)
        
        print("Inventory of female shell condition and egg codes")
        specimen_table %>%
          dplyr::filter(SEX == 2) %>%
          dplyr::group_by(SPECIES_CODE, SHELL_CONDITION, EGG_COLOR, EGG_CONDITION, CLUTCH_SIZE) %>%
          dplyr::reframe(COUNT = sum(SAMPLING_FACTOR),
                         N = n()) %>%
          as.data.frame()-> invent4
        print(invent4)
        
        print("CHECKING LENGTHS AND WIDTHS...")
        
        # 11) Any missing lengths for RKC?"
        if(unique(is.na(specimen_table$LENGTH)) == TRUE){
          print("ERROR: missing length for RKC")
        } 
        
        # 12) Any small female crab with a clutch size?"
        if(unique(filter(specimen_table, SEX != 1)$LENGTH < 65 & 
                  filter(specimen_table, SEX != 1)$CLUTCH_SIZE > 0) == TRUE){
          print("ERROR: female <65 with clutch size >0")
        } 
        
        # 13) Any small crab with old shell condition?"
        if(unique(specimen_table$LENGTH < 60 & specimen_table$SHELL_CONDITION > 2) == TRUE){
          print("ERROR: crab < 60 with shell condition >2")
        } 
        
        # 14) Any widths entered for RKC? 
        if(unique(is.na(specimen_table$WIDTH)) == FALSE) {
          print("ERROR: width entered for RKC when length is needed")
        }
        
        print("What are the minimum and maximum lengths reported by sex?")
        specimen_table %>%
          dplyr::group_by(SPECIES_CODE, SEX) %>%
          dplyr::reframe(MIN_LENGTH = min(LENGTH),
                         MAX_LENGTH = max(LENGTH)) %>%
          as.data.frame() -> A15
        print(A15)
        
        print("CHECKING DISEASE CODES...")
        
        # 16) Any black mat recorded but without % coverage entry?
        if(unique(is.na(specimen_table$DISEASE_CODE) == FALSE &
                  specimen_table$DISEASE_CODE == 1 & (is.na(specimen_table$DISEASE_DORSAL) == TRUE &
                                                      is.na(specimen_table$DISEASE_LEGS) == TRUE &
                                                      is.na(specimen_table$DISEASE_VENTRAL) == TRUE)) == TRUE){
          print("ERROR: black mat recorded without % coverage")
        } 
        
        # 17) Any bitter crab recorded for RKC and/or any bitter crab recorded with entries in % coverage?
        if(unique(is.na(specimen_table$DISEASE_CODE) == FALSE & specimen_table$DISEASE_CODE == 2 |
                  is.na(specimen_table$DISEASE_CODE) == FALSE & specimen_table$DISEASE_CODE == 2 & 
                  (is.na(specimen_table$DISEASE_DORSAL) == FALSE &
                   is.na(specimen_table$DISEASE_LEGS) == FALSE &
                   is.na(specimen_table$DISEASE_VENTRAL) == FALSE)) == TRUE){
          print("ERROR: bitter crab recorded for RKC and/or bitter crab recorded with % coverage")
        } 
        
        # 18) Any disease code not recorded but entries in % coverage?
        if(unique(is.na(specimen_table$DISEASE_CODE) == TRUE & (is.na(specimen_table$DISEASE_DORSAL) == FALSE &
                                                                is.na(specimen_table$DISEASE_LEGS) == FALSE &
                                                                is.na(specimen_table$DISEASE_VENTRAL) == FALSE)) == TRUE){
          print("ERROR: disease code not recorded but % cover entered")
        } 
        
        
        # 19) Any disease codes >9? 
        if(unique(is.na(specimen_table$DISEASE_CODE) == "FALSE" & specimen_table$DISEASE_CODE>9) == TRUE){
          print("ERROR: disease code >9")
        }
        
        print("CHECKING SAMPLING FACTOR...")
        
        # 20) What is the maximum sampling factor by by sex?
        print("What is the maximum sampling factor by sex?")
        specimen_table %>%
          dplyr::group_by(SPECIES_CODE, SEX) %>%
          dplyr::reframe(MIN_SAMPLING_FACTOR = min(SAMPLING_FACTOR),
                         MAX_SAMPLING_FACTOR= max(SAMPLING_FACTOR)) %>%
          as.data.frame() -> A20
        print(A20)
        
        # 21) Any sampling factors < 1?"
        if(unique(specimen_table$SAMPLING_FACTOR < 1) == TRUE){
          print("ERROR: minimum sampling factor < 1")
        } 
        
        
        print("CHECKING HAUL, STATION, AND BUOY IDs...")
        pot_cpue %>%
          dplyr::select(VESSEL, SPN, POT_ID, BUOY, LAT_DD, LON_DD) %>%
          distinct() %>%
          as.data.frame()-> xx
        
        potlifts %>%
          dplyr::select(VESSEL, SPN, POT_ID, BUOY, LAT_DD, LON_DD) %>%
          distinct() -> yy
        
        if(TRUE %in% is.na(suppressMessages(right_join(xx, yy, keep = TRUE))) == TRUE){
          print("ERROR: pot, station, and/or buoy IDs do not match between potlifts table and pot cpue table")
        }
        
        print("COMPARING SPECIMEN TABLE WITH CATCH SUMMARY...")
        specimen_table %>%
          dplyr::group_by(CRUISE, VESSEL, POT_ID, SPECIES_CODE) %>%
          dplyr::reframe(NUMBER_CRAB = sum(SAMPLING_FACTOR),
                         N_ENTRIES = n()) -> spec_sum
        
        catch_summary %>%
          dplyr::select(CRUISE, VESSEL, POT_ID, SPECIES_CODE, NUMBER_CRAB, N_ENTRIES) -> catch_sum
        
        # 22) Does the number of crab and number of entries match between the specimen table and catch summary?
        if(TRUE %in% is.na(suppressMessages(right_join(spec_sum, catch_sum))) == TRUE){
          print("ERROR: number of crab and number of entries do not 
                 match between summarized specimen table and catch summary")
        } 
        
        print("CHECKING COORDINATES...")
        
        # Transform potlifts data to correct crs
        potlifts %>%
          sf::st_as_sf(coords = c(x = "LON_DD", y = "LAT_DD"), crs = sf::st_crs(4326)) %>%
          sf::st_transform(crs = map.crs)-> potlifts_mapdat 
        
        #set up plotting features
        map_layers <- map_layers <- readRDS("./Data/akgfmaps_layers.rds")
        
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
        
        # Map coordinates
        ggplot() +
          geom_sf(data = st_as_sf(BB_strata), fill = NA, color = "black", linewidth = 1) +
          geom_sf(data = st_as_sf(RKCSA_sub), fill = NA, color = "red", alpha= 0.5, linewidth = 1) +
          geom_sf(data = st_as_sf(RKCSA), fill = NA, color = "red", alpha =0.5, linewidth = 1) +
          geom_sf(data = map_layers$akland, fill = "grey80") +
          geom_sf(data = potlifts_mapdat, shape = 19, size = 1.5, colour = "black", stat="identity", 
                  position="identity")+
          scale_x_continuous(breaks = breaks.x, labels = paste0(breaks.x, "°W")) + 
          scale_y_continuous(breaks = breaks.y, labels = paste0(breaks.y, "°N")) +
          ggtitle("Coordinate check")+
          coord_sf(xlim = plot.boundary$x,
                   ylim = plot.boundary$y) +
          theme_bw() -> coords
        
        print(coords)
        
        
    }    
      
  # Run function
    error_chk(specimen_table, catch_summary, potlifts, pot_cpue) 

# BYCATCH -----------------------------------------------------------------------
  # Process data
  bycatch <- rbind(read.csv("./Data/SS_BYCATCH.csv") %>%
                              select(!c(RKC.Male.in.NonSurveyPots, RKC.Female.in.NonSurvey.Pots,
                                        Pollock, Starry.Flounder)), 
                   read.csv("./Data/SB_BYCATCH .csv")) %>%
              mutate(SPN = as.numeric(as.character(SPN)), VESSEL = ifelse(VESSEL == 162, "Summer Bay", "Silver Spray")) %>%
              filter(is.na(SPN) == FALSE, !(VESSEL == "Summer Bay" & SPN %in% 300:311)) %>%
              right_join(potlifts %>% select(c(VESSEL, SPN, LON_DD, LAT_DD)), .) %>%
              replace(., is.na(.), 0) %>%
          sf::st_as_sf(coords = c(x = "LON_DD", y = "LAT_DD"), crs = sf::st_crs(4326)) %>%
          sf::st_transform(crs = map.crs)
    
    # Make map labels
    bc_labs <- data.frame(labs = c("Pacific Cod", "Yellowfin Sole"),
                          name = c("PacificCod", "YellowfinSole"))

    #Plot YFS
    ggplot() +
      geom_sf(data = st_transform(map_layers$bathymetry, map.crs), color=alpha("grey70")) +
      geom_sf(data = st_as_sf(BB_strata), fill = NA, mapping = aes(color = "black"), linewidth = 1) +
      geom_sf(data = st_as_sf(RKCSA_sub), mapping = aes(color = "red"), fill = NA, alpha= 0.9, linewidth = 1) +
      geom_sf(data = st_as_sf(RKCSA), fill = NA,  color = "red", alpha =0.5, linewidth = 1) +
      geom_sf(data = st_transform(map_layers$akland, map.crs), fill = "grey80") +
      geom_sf(data = bycatch,
              mapping = aes(size=YellowfinSole, fill = YellowfinSole, shape = YellowfinSole == 0), 
              alpha = 0.5, colour = "black")+
      scale_shape_manual(values = c('TRUE' = 4, 'FALSE' = 21), guide = "none") +
      scale_color_manual(values = c("black", "red"), 
                         labels = c("EBS Summer Survey Boundary", "Red King Crab Savings Area"),
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
      geom_sf(data = st_as_sf(BB_strata), fill = NA, mapping = aes(color = "black"), linewidth = 1) +
      geom_sf(data = st_as_sf(RKCSA_sub), mapping = aes(color = "red"), fill = NA, alpha= 0.9, linewidth = 1) +
      geom_sf(data = st_as_sf(RKCSA), fill = NA,  color = "red", alpha =0.5, linewidth = 1) +
      geom_sf(data = st_transform(map_layers$akland, map.crs), fill = "grey80") +
      geom_sf(data = bycatch,
              mapping = aes(size=PacificCod, fill = PacificCod, shape = PacificCod == 0), 
              alpha = 0.5, colour = "black")+
      scale_shape_manual(values = c('TRUE' = 4, 'FALSE' = 21), guide = "none") +
      scale_color_manual(values = c("black", "red"), 
                         labels = c("EBS Summer Survey Boundary", "Red King Crab Savings Area"),
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
    
    # Save figures
    ggsave(plot = yfs_map, "./Figures/yfs_map.png", height=7, width=10, units="in")
    ggsave(plot = pcod_map, "./Figures/pcod_map.png", height=7, width=10, units="in")
    
# CRAB COUNTS BY ZONE --------------------------------------------------------------------
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
    

    rbind(RKCSA_N, RKCSA_W, RKCSA_S, RKCSA_all, RKCSA_E) %>%
      filter(MAT_SEX != c("Legal male", "Sublegal male")) %>%
      group_by(MAT_SEX) %>%
      reframe(xx <- sum(SUM_COUNT)) -> tt
    

# CPT PLOTS -------------------------------------------------------------------------
  #1A: Empty map with bathy, CPS1 outline, and RKCSA -------
    ggplot() +
      geom_sf(data = st_transform(map_layers$bathymetry, map.crs), color=alpha("grey70")) +
      geom_sf(data = st_as_sf(CPS1_bound), fill = NA, aes(color = "black"), linewidth = 1)+
      geom_sf(data = st_as_sf(RKCSA_sub), mapping = aes(color = "red"), fill = NA, alpha= 0.9, linewidth = 1) +
      geom_sf(data = st_as_sf(RKCSA), fill = NA,  color = "red", alpha =0.5, linewidth = 1) +
      geom_sf(data = st_transform(map_layers$akland, map.crs), fill = "grey80") +
      scale_x_continuous(breaks = c(-165, -160), labels = paste0(c(165, 160), "°W"))+
      scale_y_continuous(breaks = c(56, 58), labels = paste0(c(56, 58), "°N"))+
      labs(title = "2023 BBRKC Collaborative Pot Sampling", subtitle = "Survey extent")+
      scale_color_manual(values = c("black", "red"), 
                         labels = c("CPS1 survey boundary", "Red King Crab Savings Area"),
                         name = "") +
      coord_sf(xlim = plot.boundary$x,
               ylim = plot.boundary$y) +
      #ggplot2::scale_fill_manual(name = "Temperature (?C)", values = viridis_pal(option = viridis_option)(n_temp_breaks),
      #labels = c(expression(""<=-1), "-0.9-0", "0.1-1", "1.1-2", "2.1-3", "3.1-4",
      #"4.1-5", "5.1-6", "6.1-7", "7.1-8", ">8.1"), drop = FALSE) +
      scale_fill_viridis(name = "Temperature (°C)", option = "plasma")+
      geom_sf_text(sf::st_as_sf(data.frame(lab= c("50m", "100m"), 
                                           x = c(-161.5, -165), y = c(58.3, 56.1)),
                                coords = c(x = "x", y = "y"), crs = sf::st_crs(4326)) %>%
                     sf::st_transform(crs = map.crs),
                   mapping = aes(label = lab))+
      guides(color = guide_legend(nrow = 2)) +
      theme_bw() +
      theme(axis.title = element_blank(),
            axis.text = element_text(size = 10),
            legend.text = element_text(size = 10),
            legend.title = element_text(size = 10),
            legend.position = "bottom",
            legend.direction = "horizontal",
            plot.title = element_text(face = "bold", size = 15),
            plot.subtitle = element_text(size = 12)) -> cpt_1A
    
    ggsave(plot = cpt_1A, "./Figures/emptymap.png", height=7, width=10, units="in")
    
  #1B: Temp map with bathy, CPS1 outline, and RKCSA -------
    map_layers <- readRDS("./Data/akgfmaps_layers.rds")
    map.crs <- "EPSG:3338"
    
    temp <- read.csv("./Data/2023_BBRKC_ALL_TEMPS_SAL.csv") %>%
      filter(is.na(AveTemp) == FALSE) %>%
      rename()
    
    # Load boundary
    st_read("./Data/CPS1_survey_boundary.shp") %>%
      st_transform(map.crs) %>%
      vect() -> CPS1_bound
    
    # Transform catch and tagging data to correct crs
    pot_cpue %>%
      sf::st_as_sf(coords = c(x = "LON_DD", y = "LAT_DD"), crs = sf::st_crs(4326)) %>%
      sf::st_transform(crs = map.crs)-> pot_cpue_mapdat 
    
    # Make raster for interpolation
    plot.boundary.untrans <- data.frame(y = c(54.5, 58.5), 
                                        x = c(-164.8, -159)) # plot boundary unprojected
    
    plot.boundary <- plot.boundary.untrans %>%
      sf::st_as_sf(coords = c(x = "x", y = "y"), crs = sf::st_crs(4326)) %>%
      sf::st_transform(crs = map.crs) %>%
      sf::st_coordinates() %>%
      as.data.frame() %>%
      dplyr::rename(x = X, y = Y) # plot boundary projected
    
    # Set up interpolation raster
    sp_interp.raster <- raster::raster(st_as_sf(CPS1_bound), res = 1000)
    
    st_read(survey_gdb,layer="BristolBaySurveyStrata") %>%
      st_transform(map.crs) -> BB_rast
    
    # Transform data for interpolation
    in.crs <- "+proj=longlat +datum=NAD83"
    sp_interp.df <- unique(temp)
    sp::coordinates(sp_interp.df) <- c(x = "Longitude_pot", y = "Latitude_pot")
    sp::proj4string(sp_interp.df) <- sp::CRS(in.crs)
    sp_interp.df <- sp::spTransform(sp_interp.df, sp::CRS(map.crs))
    
    # Set up a new IDW for ordinary kriging 
    idw_vgm_fit <- gstat::gstat(formula = AveTemp ~ 1, 
                                locations = sp_interp.df, 
                                nmax = Inf)
    
    # Ordinary Kriging: VGM
    ste.vgfit <- gstat::fit.variogram(variogram(idw_vgm_fit), 
                                      vgm(c( "Bes")))
    
    ste_fit <- gstat::gstat(formula = AveTemp ~ 1, 
                            locations = sp_interp.df, 
                            model = ste.vgfit, 
                            nmax = Inf)
    
    ste.predict <- predict(ste_fit, as(sp_interp.raster, "SpatialGrid"))
    
    st_read("./Data/CPS1_survey_boundary.shp") %>%
      st_transform(map.crs) -> CPS1_poly
    
    # write unmasked surfaces to raster, stacked by year
    ste.predict %>%
      raster::raster(.) %>%
      mask(CPS1_poly) ->  temp_rast
    
    # extract interpolated data from raster to data frame
    coords<-coordinates(temp_rast)
    
    temp_df<-na.omit(data.frame(coords, temperature = temp_rast@data@values))
    
    temp_breaks <- c(-Inf, seq(-1,4,1), Inf)
    viridis_option <- "H" # viridis turbo palette
    n_temp_breaks <- length(temp_breaks)-1
    
    # Filter potlifts data by pots that had adfg temp loggers
    potlifts %>%
      mutate(NOTES = ifelse((VESSEL == "Silver Spray" & SPN == 44), paste("adfg ", NOTES), NOTES)) %>% #missing adfg note
      filter(NOTES != "" & !grepl(" hotspot release", NOTES) & NOTES != "hotspot release") %>%
      temp %>%
      sf::st_as_sf(coords = c(x = "Longitude_pot", y = "Latitude_pot"), crs = sf::st_crs(4326)) %>%
      sf::st_transform(crs = map.crs) -> temploggers
    
    # Plot
    ggplot() +
      geom_tile(data = temp_df, aes(x = x, y = y, fill = temperature))+
      scale_fill_viridis(name = "TEMPERATURE (°C)", option = "plasma", 
                         guide = guide_colorbar(title.position = "top"))+
      geom_sf(data = st_transform(map_layers$bathymetry, map.crs), color=alpha("grey70")) +
      geom_sf(data = st_as_sf(CPS1_bound), fill = NA, aes(color = "black"), linewidth = 1)+
      #geom_sf(data = st_as_sf(BB_strata), fill = NA, mapping = aes(color = "black"), linewidth = 1) +
      geom_sf(data = st_as_sf(RKCSA_sub), mapping = aes(color = "red"), fill = NA, alpha= 0.9, linewidth = 1) +
      geom_sf(data = st_as_sf(RKCSA), fill = NA,  color = "red", alpha =0.5, linewidth = 1) +
      geom_sf(data = temploggers, color = "white", size = 2)+
      geom_sf(data = st_transform(map_layers$akland, map.crs), fill = "grey80") +
      scale_x_continuous(breaks = c(-165, -160), labels = paste0(c(165, 160), "°W"))+
      scale_y_continuous(breaks = c(56, 58), labels = paste0(c(56, 58), "°N"))+
      labs(title = "2023 BBRKC Collaborative Pot Sampling", subtitle = "Temperature logger locations")+
      scale_color_manual(values = c("black", "red"), 
                         labels = c("CPS1 survey boundary", "Red King Crab Savings Area"),
                         name = "") +
      coord_sf(xlim = plot.boundary$x,
               ylim = plot.boundary$y) +
      #ggplot2::scale_fill_manual(name = "Temperature (?C)", values = viridis_pal(option = viridis_option)(n_temp_breaks),
      #labels = c(expression(""<=-1), "-0.9-0", "0.1-1", "1.1-2", "2.1-3", "3.1-4",
      #"4.1-5", "5.1-6", "6.1-7", "7.1-8", ">8.1"), drop = FALSE) +
      geom_sf_text(sf::st_as_sf(data.frame(lab= c("50m", "100m"), 
                                           x = c(-161.5, -165), y = c(58.3, 56.1)),
                                coords = c(x = "x", y = "y"), crs = sf::st_crs(4326)) %>%
                     sf::st_transform(crs = map.crs),
                   mapping = aes(label = lab))+
      guides(color = guide_legend(nrow = 2)) +
      theme_bw() +
      theme(axis.title = element_blank(),
            axis.text = element_text(size = 10),
            legend.text = element_text(size = 10),
            legend.title = element_text(size = 10),
            legend.position = "bottom",
            legend.direction = "horizontal",
            plot.title = element_text(face = "bold", size = 15),
            plot.subtitle = element_text(size = 12)) -> temp_map_continuous
    
    ggsave(plot = temp_map_continuous, "./Figures/tempcont.png", height=7, width=10, units="in")
    
  
  #1C: All crab count and temp map --------------------------------------------------
    pot_cpue %>%
      group_by(SPN, LAT_DD, LON_DD) %>%
      reframe(COUNT = sum(COUNT)) %>%
      sf::st_as_sf(coords = c(x = "LON_DD", y = "LAT_DD"), crs = sf::st_crs(4326)) %>%
      sf::st_transform(crs = map.crs) -> pot_cpue_allcrabmapdat
    
    ggplot() +
      geom_tile(data = temp_df, aes(x = x, y = y, fill = temperature))+
      scale_fill_viridis(name = "TEMPERATURE (°C)", option = "plasma", 
                         guide = guide_colorbar(title.position = "top"))+
      ggnewscale::new_scale_fill()+
      geom_sf(data = st_transform(map_layers$bathymetry, map.crs), color=alpha("grey70")) +
      geom_sf(data = st_as_sf(CPS1_bound), fill = NA, aes(color = "black"), linewidth = 1)+
      #geom_sf(data = st_as_sf(BB_strata), fill = NA, mapping = aes(color = "black"), linewidth = 1) +
      geom_sf(data = st_as_sf(RKCSA_sub), mapping = aes(color = "red"), fill = NA, alpha= 0.9, linewidth = 1) +
      geom_sf(data = st_as_sf(RKCSA), fill = NA,  color = "red", alpha =0.5, linewidth = 1) +
      geom_sf(data = st_transform(map_layers$akland, map.crs), fill = "grey80") +
      geom_sf(data = pot_cpue_allcrabmapdat,
              mapping = aes(size=COUNT, fill = COUNT, shape = COUNT == 0), alpha = 0.5, colour = "black")+
      scale_shape_manual(values = c('TRUE' = 4, 'FALSE' = 21), guide = "none")+
      scale_color_manual(values = c("black", "red"), 
                         labels = c("CPS1 survey boundary", "Red King Crab Savings Area"),
                         name = "") +
      scale_size_continuous(range = c(2, 10), limits = c(0, max(pot_cpue_allcrabmapdat$COUNT)), 
                            breaks =seq(0, max(pot_cpue_allcrabmapdat$COUNT), by = 50))+ 
      scale_fill_gradientn(breaks = seq(0, max(pot_cpue_allcrabmapdat$COUNT), by = 50),
                           limits = c(0, max(pot_cpue_allcrabmapdat$COUNT)), 
                           colors = c("gray", rev(pal[5:length(pal)])))+
      scale_x_continuous(breaks = c(-165, -160), labels = paste0(c(165, 160), "°W"))+
      scale_y_continuous(breaks = c(56, 58), labels = paste0(c(56, 58), "°N"))+
      labs(title = "2023 BBRKC Collaborative Pot Sampling", subtitle = "All crab")+
      guides(size = guide_legend(title.position = "top", nrow = 2, order = 1, override.aes = list(shape = c(4, rep(21, 6)))),
             fill = guide_legend(order = 1),
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
            plot.subtitle = element_text(size = 12)) -> all_crab.temp_map
    
    ggsave(plot = all_crab.temp_map, "./Figures/all_crab.temp_map.png", height=7, width=10, units="in")
  
  #1D: Crab count by mat/sex and temp map -----------------------------------------------------------
    mat_sex_combos %>%
      purrr::map(~ggplot() +
                   geom_tile(data = temp_df, aes(x = x, y = y, fill = temperature))+
                   scale_fill_viridis(name = "TEMPERATURE (°C)", option = "plasma", 
                                      guide = guide_colorbar(title.position = "top"))+
                   ggnewscale::new_scale_fill()+
                   geom_sf(data = st_transform(map_layers$bathymetry, map.crs), color=alpha("grey70")) +
                   geom_sf(data = st_as_sf(CPS1_bound), fill = NA, aes(color = "black"), linewidth = 1)+
                   #geom_sf(data = st_as_sf(BB_strata), fill = NA, mapping = aes(color = "black"), linewidth = 1) +
                   geom_sf(data = st_as_sf(RKCSA_sub), mapping = aes(color = "red"), fill = NA, alpha= 0.9, linewidth = 1) +
                   geom_sf(data = st_as_sf(RKCSA), fill = NA,  color = "red", alpha =0.5, linewidth = 1) +
                   geom_sf(data = st_transform(map_layers$akland, map.crs), fill = "grey80") +
                   geom_sf(data = filter(pot_cpue_mapdat, MAT_SEX == .x),
                           mapping = aes(size=COUNT, fill = COUNT, shape = COUNT == 0), alpha = 0.5, colour = "black")+
                   scale_shape_manual(values = c('TRUE' = 4, 'FALSE' = 21), guide = "none")+
                   scale_color_manual(values = c("black", "red"), 
                                      labels = c("CPS1 survey boundary", "Red King Crab Savings Area"),
                                      name = "") +
                   scale_size_continuous(range = c(2, 10), limits = c(0, max(pot_cpue_mapdat$COUNT)), 
                                         breaks =seq(0, max(pot_cpue_mapdat$COUNT), by = 25))+ 
                   scale_fill_gradientn(breaks = seq(0, max(pot_cpue_mapdat$COUNT), by = 25),
                                        limits = c(0, max(pot_cpue_mapdat$COUNT)), 
                                        colors = c("gray", rev(pal[5:length(pal)])))+
                   scale_x_continuous(breaks = c(-165, -160), labels = paste0(c(165, 160), "°W"))+
                   scale_y_continuous(breaks = c(56, 58), labels = paste0(c(56, 58), "°N"))+
                   labs(title = "2023 BBRKC Collaborative Pot Sampling", subtitle = paste(filter(mat_labs, MAT_SEX == .x)$lab))+
                   guides(size = guide_legend(title.position = "top", nrow = 2, order = 1, override.aes = list(shape = c(4, rep(21, 5)))),
                          fill = guide_legend(order = 1),
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
                         plot.subtitle = element_text(size = 12))) -> BBRKC.temp.maps
    
    # Save plots
    ggsave(plot = BBRKC.temp.maps[[1]], "./Figures/BBRKC.temp_matmale.png", height=7, width=10, units="in")
    ggsave(plot = BBRKC.temp.maps[[2]], "./Figures/BBRKC.temp_immale.png", height=7, width=10, units="in")
    ggsave(plot = BBRKC.temp.maps[[3]], "./Figures/BBRKC.temp_matfem.png", height=7, width=10, units="in")
    ggsave(plot = BBRKC.temp.maps[[4]], "./Figures/BBRKC.temp_imfem.png", height=7, width=10, units="in")
    ggsave(plot = BBRKC.temp.maps[[5]], "./Figures/BBRKC.temp_legalmale.png", height=7, width=10, units="in")
    ggsave(plot = BBRKC.temp.maps[[6]], "./Figures/BBRKC.temp_sublegalmale.png", height=7, width=10, units="in")
    
  #1E: Bycatch and temp map ---------------------------------------------------------
    #Plot YFS
    ggplot() +
      geom_tile(data = temp_df, aes(x = x, y = y, fill = temperature))+
      scale_fill_viridis(name = "TEMPERATURE (°C)", option = "plasma", 
                         guide = guide_colorbar(title.position = "top"))+
      ggnewscale::new_scale_fill()+
      geom_sf(data = st_transform(map_layers$bathymetry, map.crs), color=alpha("grey70")) +
      geom_sf(data = st_as_sf(CPS1_bound), fill = NA, aes(color = "black"), linewidth = 1)+
      geom_sf(data = st_as_sf(RKCSA_sub), mapping = aes(color = "red"), fill = NA, alpha= 0.9, linewidth = 1) +
      geom_sf(data = st_as_sf(RKCSA), fill = NA,  color = "red", alpha =0.5, linewidth = 1) +
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
      guides(size = guide_legend(title.position = "top", title = "COUNT", order = 1, nrow = 2, override.aes = list(shape = c(4, rep(21, 4)))),
             fill = guide_legend(title = "COUNT", order = 1),
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
            plot.subtitle = element_text(size = 12)) -> yfs.temp_map
    
    # Plot Pcod
    ggplot() +
      geom_tile(data = temp_df, aes(x = x, y = y, fill = temperature))+
      scale_fill_viridis(name = "TEMPERATURE (°C)", option = "plasma", 
                         guide = guide_colorbar(title.position = "top"))+
      ggnewscale::new_scale_fill()+
      geom_sf(data = st_transform(map_layers$bathymetry, map.crs), color=alpha("grey70")) +
      geom_sf(data = st_as_sf(CPS1_bound), fill = NA, aes(color = "black"), linewidth = 1)+
      #geom_sf(data = st_as_sf(BB_strata), fill = NA, mapping = aes(color = "black"), linewidth = 1) +
      geom_sf(data = st_as_sf(RKCSA_sub), mapping = aes(color = "red"), fill = NA, alpha= 0.9, linewidth = 1) +
      geom_sf(data = st_as_sf(RKCSA), fill = NA,  color = "red", alpha =0.5, linewidth = 1) +
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
      guides(size = guide_legend(title.position = "top",title = "COUNT", order = 1, nrow = 2, override.aes = list(shape = c(4, rep(21, 3)))),
             fill = guide_legend(title = "COUNT", order = 1),
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
            plot.subtitle = element_text(size = 12)) -> pcod.temp_map
    
    # Save figures
    ggsave(plot = yfs.temp_map, "./Figures/yfs.temp_map.png", height=7, width=10, units="in")
    ggsave(plot = pcod.temp_map, "./Figures/pcod.temp_map.png", height=7, width=10, units="in")
  #1F: High resolution bathymetry and temp map ---------------------------
    ggplot() +
      geom_tile(data = temp_df, aes(x = x, y = y, fill = temperature))+
      scale_fill_viridis(name = "TEMPERATURE (°C)", option = "plasma", 
                         guide = guide_colorbar(title.position = "top"))+
      ggnewscale::new_scale_fill()+
      geom_sf(data = CPS1_bathy, color=alpha("white")) +
      geom_sf(data = st_as_sf(CPS1_bound), fill = NA, aes(color = "black"), linewidth = 1)+
      #geom_sf(data = st_as_sf(BB_strata), fill = NA, mapping = aes(color = "black"), linewidth = 1) +
      geom_sf(data = st_as_sf(RKCSA_sub), mapping = aes(color = "red"), fill = NA, alpha= 0.9, linewidth = 1) +
      geom_sf(data = st_as_sf(RKCSA), fill = NA,  color = "red", alpha =0.5, linewidth = 1) +
      geom_sf(data = st_transform(map_layers$akland, map.crs), fill = "grey80") +
      scale_shape_manual(values = c('TRUE' = 4, 'FALSE' = 21), guide = "none")+
      scale_color_manual(values = c("black", "red"), 
                         labels = c("CPS1 survey boundary", "Red King Crab Savings Area"),
                         name = "") +
      scale_size_continuous(range = c(2, 10), limits = c(0, max(pot_cpue_allcrabmapdat$COUNT)), 
                            breaks =seq(0, max(pot_cpue_allcrabmapdat$COUNT), by = 50))+ 
      scale_fill_gradientn(breaks = seq(0, max(pot_cpue_allcrabmapdat$COUNT), by = 50),
                           limits = c(0, max(pot_cpue_allcrabmapdat$COUNT)), 
                           colors = c("gray", rev(pal[5:length(pal)])))+
      scale_x_continuous(breaks = c(-165, -160), labels = paste0(c(165, 160), "°W"))+
      scale_y_continuous(breaks = c(56, 58), labels = paste0(c(56, 58), "°N"))+
      labs(title = "2023 BBRKC Collaborative Pot Sampling", subtitle = "High resolution bathymetry")+
      guides(size = guide_legend(title.position = "top", nrow = 2, order = 1, override.aes = list(shape = c(4, rep(21, 6)))),
             fill = guide_legend(order = 1),
             color = guide_legend(nrow = 2))+
      coord_sf(xlim = plot.boundary$x,
               ylim = plot.boundary$y) +
      geom_shadowtext(data = (sf::st_as_sf(data.frame(lab= c("35m", "45m", "55m", "65m", "75m", "85m"), 
                                           x = c(-160.6, -161.6, -161.6, -161.2, -163.5, -163.6), 
                                           y = c(56.1, 56, 57.5, 57.2, 56.76, 56.38)),
                                coords = c(x = "x", y = "y"), crs = sf::st_crs(4326)) %>%
                     sf::st_transform(crs = map.crs) %>%
                     cbind(st_coordinates(.))),
                   mapping = aes(label = lab, x = X, y = Y), color = "black", bg.color = "white", size = 3.5)+
      theme_bw() +
      theme(axis.title = element_blank(),
            axis.text = element_text(size = 10),
            legend.text = element_text(size = 10),
            legend.title = element_text(size = 10),
            legend.position = "bottom",
            legend.direction = "horizontal",
            plot.title = element_text(face = "bold", size = 15),
            plot.subtitle = element_text(size = 12)) -> temp.highresbathy_map
     
    ggsave(plot = temp.highresbathy_map, "./Figures/temp.highresbathy_map.png", height=7, width=10, units="in")
    
  # 1G: High resolution bathymetry, temp, and all crab ------------------------------------------------
    ggplot() +
      geom_tile(data = temp_df, aes(x = x, y = y, fill = temperature))+
      scale_fill_viridis(name = "TEMPERATURE (°C)", option = "plasma", 
                         guide = guide_colorbar(title.position = "top"))+
      ggnewscale::new_scale_fill()+
      geom_sf(data = CPS1_bathy, color=alpha("grey70")) +
      geom_sf(data = st_as_sf(CPS1_bound), fill = NA, aes(color = "black"), linewidth = 1)+
      #geom_sf(data = st_as_sf(BB_strata), fill = NA, mapping = aes(color = "black"), linewidth = 1) +
      geom_sf(data = st_as_sf(RKCSA_sub), mapping = aes(color = "red"), fill = NA, alpha= 0.9, linewidth = 1) +
      geom_sf(data = st_as_sf(RKCSA), fill = NA,  color = "red", alpha =0.5, linewidth = 1) +
      geom_sf(data = st_transform(map_layers$akland, map.crs), fill = "grey80") +
      geom_sf(data = pot_cpue_allcrabmapdat,
              mapping = aes(size=COUNT, fill = COUNT, shape = COUNT == 0), alpha = 0.5, colour = "black")+
      scale_shape_manual(values = c('TRUE' = 4, 'FALSE' = 21), guide = "none")+
      scale_color_manual(values = c("black", "red"), 
                         labels = c("CPS1 survey boundary", "Red King Crab Savings Area"),
                         name = "") +
      scale_size_continuous(range = c(2, 10), limits = c(0, max(pot_cpue_allcrabmapdat$COUNT)), 
                            breaks =seq(0, max(pot_cpue_allcrabmapdat$COUNT), by = 50))+ 
      scale_fill_gradientn(breaks = seq(0, max(pot_cpue_allcrabmapdat$COUNT), by = 50),
                           limits = c(0, max(pot_cpue_allcrabmapdat$COUNT)), 
                           colors = c("gray", rev(pal[5:length(pal)])))+
      scale_x_continuous(breaks = c(-165, -160), labels = paste0(c(165, 160), "°W"))+
      scale_y_continuous(breaks = c(56, 58), labels = paste0(c(56, 58), "°N"))+
      labs(title = "2023 BBRKC Collaborative Pot Sampling", subtitle = "All crab with high resolution bathymetry")+
      guides(size = guide_legend(title.position = "top", nrow = 2, order = 1, override.aes = list(shape = c(4, rep(21, 6)))),
             fill = guide_legend(order = 1),
             color = guide_legend(nrow = 2))+
      coord_sf(xlim = plot.boundary$x,
               ylim = plot.boundary$y) +
      #geom_shadowtext(data = (sf::st_as_sf(data.frame(lab= c("35m", "45m", "55m", "65m", "75m", "85m"), 
                                                      #x = c(-160.6, -161.6, -161.6, -161.2, -163.5, -163.6), 
                                                      #y = c(56.1, 56, 57.5, 57.2, 56.76, 56.38)),
                                           #coords = c(x = "x", y = "y"), crs = sf::st_crs(4326)) %>%
                                #sf::st_transform(crs = map.crs) %>%
                                #cbind(st_coordinates(.))),
                      #mapping = aes(label = lab, x = X, y = Y), color = "black", bg.color = "white", size = 3.5)+
      theme_bw() +
      theme(axis.title = element_blank(),
            axis.text = element_text(size = 10),
            legend.text = element_text(size = 10),
            legend.title = element_text(size = 10),
            legend.position = "bottom",
            legend.direction = "horizontal",
            plot.title = element_text(face = "bold", size = 15),
            plot.subtitle = element_text(size = 12)) -> all_crab.temp.highresbathy_map
    
    ggsave(plot = all_crab.temp.highresbathy_map, "./Figures/all_crab.temp.highresbathy_map.png", height=7, width=10, units="in")
    