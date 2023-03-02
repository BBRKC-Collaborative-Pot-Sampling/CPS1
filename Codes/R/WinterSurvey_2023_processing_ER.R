# PURPOSE -------------------------------------------------------------------------------------------------------------------
  # 1) To automate processing specimen tables and catch summaries from Winter 2023 pot survey for BBRKC
  # 2) To run error checks on processed specimen and catch summaries
  # 3) To calculate and map cpue by pot and mat/sex category for BBRKC

  # !!Need to add filter for gear testing codes and bad haul codes once known!! 

# INSTALL PACKAGES ----------------------------------------------------------------------------------------------------------
  #devtools::install_github("sean-rohan-NOAA/akgfmaps", build_vignettes = TRUE)
  #devtools::install_github("afsc-gap-products/coldpool")
  
  #install.packages(c("tidyverse", "gsubfn", "terra", "rgdal", "colorRamps", "coldpool"))

  #install_version("tidyverse", version = "1.3.2")
  
  #install_version("dplyr", version = "1.1.0")

# LOAD PACKAGES -------------------------------------------------------------------------------------------------------------
  library(tidyverse)
  library(akgfmaps)
  library(gsubfn)
  library(terra)
  library(rgdal)
  library(colorRamps)
  library(coldpool)

# LOAD DATA -----------------------------------------------------------------------------------------------------------------
  
  # Load summary catch and specimen tables
      catch <- list.files("./Data/FTP/Catch - FTP/") %>%
        purrr::map_df(~read.csv(paste0("./Data/FTP/Catch - FTP/", .x)))
      
      # Process catch data for Oracle and save
      catch %>%
        dplyr::select(HAUL_ID, SPECIES_CODE, NUMBER_CRAB, VESSEL, CRUISE, HAUL) %>%
        write.csv("./DataforOracle/Processed_Catch_Summary.csv")
      
      
      specimen <- list.files("./Data/FTP/Specimen - FTP/") %>%
        purrr::map_df(~read.csv(paste0("./Data/FTP/Specimen - FTP/", .x)))
      
  # Load raw data for processing below   
      raw_sample <- list.files("./Data/FTP/RawData - FTP/", pattern = "_SAMPLE_0") %>% # RECORDS of SAMPLE INFO
        purrr::map_df(~read.csv(paste0("./Data/FTP/RawData - FTP/", .x))) #E.G. SEX, SPECIES
      
      raw_sample_values <- list.files("./Data/FTP/RawData - FTP/", pattern = "_SAMPLE_VALUES") %>% #RECORDS OF # TOSSED
        purrr::map_df(~read.csv(paste0("./Data/FTP/RawData - FTP/", .x))) %>%
        mutate(TOSSED = ifelse(is.na(COUNT) == FALSE, COUNT,0)) %>%
        group_by(HAUL_ID, CATCH_SAMPLE_ID) %>%
        dplyr::reframe(TOSSED = sum(TOSSED)) 
      
      raw_specimen <- list.files("./Data/FTP/RawData - FTP/", pattern = "_SPECIMEN_0") %>% 
        purrr::map_df(~read.csv(paste0("./Data/FTP/RawData - FTP/", .x))) %>%
        dplyr::select(HAUL_ID, SPECIMEN_ID, CATCH_SAMPLE_ID, SPECIES_CODE)
      
      raw_specimen_bio<- list.files("./Data/FTP/RawData - FTP/", pattern = "_SPECIMEN_BIOMETRICS") %>% 
        purrr::map_df(~read.csv(paste0("./Data/FTP/RawData - FTP/", .x))) 
      
  # Read in potlifts and tagging data
      potlifts <- list.files("./Data/FTP/", pattern = "POTLIFTS", ignore.case = TRUE) %>% #MAY NEED TO CHANGE
        purrr::map_df(~read.csv(paste0("./Data/FTP/", .x))) 
      
      tagging <- list.files("./Data/FTP/", pattern = "TAGGING", ignore.case = TRUE) %>% #MAY NEED TO CHANGE
        purrr::map_df(~read.csv(paste0("./Data/FTP/", .x))) 
    
  # Read in spatial layers for mapping purposes 
      # Set crs
      map.crs <- coldpool:::ebs_proj_crs 
      
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
      
# PROCESS DATA ----------------------------------------------------------------------------------------------------------------
    
  # Calculate soak time and lat/lon in degrees decimal for all potlifts
    potlifts %>%
      dplyr::mutate(DATETIME_SET = as.POSIXct(paste(DATE_SET, TIME_SET), format = "%m/%d/%Y %H:%M"),
                    DATETIME_HAUL = as.POSIXct(paste(DATE_HAUL, TIME_HAUL), format = "%m/%d/%Y %H:%M"),
                    SOAK_TIME = as.numeric(difftime(DATETIME_HAUL, DATETIME_SET, units = "hours")),
                    LAT_DD = LAT_DEG + LAT_MIN/60,
                    LON_DD = (LON_DEG + LON_MIN/60)*-1) %>%
      dplyr::select(!c(DATETIME_SET, DATETIME_HAUL)) %>%
      dplyr::filter(is.na(VESSEL) == "FALSE") -> potlifts
    
  # Calculate lat/lon in degrees decimal for tagging release points and # of tags released per station
    tagging %>%
      dplyr::mutate(LAT_DD = LAT_DEG + LAT_MIN/60,
                    LON_DD = (LON_DEG + LON_MIN/60)*-1) %>%
      dplyr::filter(is.na(VESSEL) == "FALSE") %>%
      dplyr::group_by(VESSEL, LON_DD, LAT_DD) %>%
      dplyr::reframe(LAT_DD = LAT_DD, LON_DD = LON_DD, N = n()) -> tagging
    
  # Join raw_sample_values and raw_sample to get # tossed per haul, sex, and catch sample id
      samples <- right_join(raw_sample, raw_sample_values) %>%
        dplyr::select(HAUL_ID, CATCH_SAMPLE_ID, SPECIES_CODE, SPECIES_NAME, SEX, TOSSED)
    
  # Expand specimen biometric table, join to raw_specimen table to get catch sample ID, join with samples file to get 
  # number tossed
      raw_specimen_bio %>%
        dplyr::select(HAUL_ID, SPECIMEN_ID, BIOMETRIC_NAME, VALUE) %>%
        pivot_wider(., names_from = "BIOMETRIC_NAME", values_from = "VALUE") %>%
        dplyr::rename(SHELL_CONDITION = CRAB_SHELL_CONDITION, EGG_COLOR = CRAB_EGG_COLOR,
                      EGG_CONDITION = CRAB_EGG_CONDITION, CLUTCH_SIZE = CRAB_EGG_CLUTCH_SIZE,
                      LENGTH = CARAPACE_LENGTH) %>%
        right_join(., raw_specimen) %>%
        right_join(samples, ., by = c("HAUL_ID", "CATCH_SAMPLE_ID", "SPECIES_CODE", "SEX"), 
                   multiple = "all") -> specimen_sum
      
  # Calculate sampling factor from specimen summary table, join back with specimen_sum file to 
  # get specimen information, join with catch file to get vessel and pot #s, join with potlifts
  # file to get lat/lon, set/haul date and time for each pot (with positive catch)
      specimen_sum %>%
        dplyr::group_by(HAUL_ID, CATCH_SAMPLE_ID, SEX) %>%
        dplyr::reframe(KEPT = n(),
                  TOSSED = TOSSED,
                  SAMPLING_FACTOR = (KEPT + TOSSED)/KEPT) %>%
        distinct() %>%
        dplyr::right_join(specimen_sum, by = c("HAUL_ID", "CATCH_SAMPLE_ID", "SEX", "TOSSED"),
                          multiple = "all") %>%
        dplyr::right_join(catch,., by = c("SPECIES_CODE", "HAUL_ID"),
                          multiple = "all") %>%
        dplyr::left_join(specimen) %>%
        distinct() %>%
        dplyr::rename(POT = HAUL) %>% 
        dplyr::right_join(potlifts, ., by = c("VESSEL", "POT"), multiple = "all") %>%
        dplyr::select(CRUISE, VESSEL, POT, BUOY, LAT_DD, LON_DD, DATE_HAUL, TIME_HAUL, SOAK_TIME, DEPTH_F,
                      SPECIES_CODE, SEX, LENGTH, WIDTH, SAMPLING_FACTOR, SHELL_CONDITION, EGG_COLOR, EGG_CONDITION, 
                      CLUTCH_SIZE, WEIGHT, DISEASE_CODE, DISEASE_DORSAL, DISEASE_VENTRAL, DISEASE_LEGS,  
                      CHELA_HEIGHT, MERUS_LENGTH, COMMENTS) -> specimen_table
      
  # Process specimen table for Oracle, save
      specimen_table %>%
        dplyr::select(!c(LAT_DD, LON_DD, DATE_HAUL, TIME_HAUL, SOAK_TIME, DEPTH_F)) %>%
        dplyr::rename(HAUL = POT, STATION = BUOY) %>%
        write.csv("./DataForOracle/Processed_Specimen_Data.csv")
      
  # Update catch summary table with new crab #s from sampling factor
      specimen_table %>%
        dplyr::group_by(CRUISE, VESSEL, POT, SPECIES_CODE) %>%
        dplyr::reframe(NUMBER_CRAB = sum(SAMPLING_FACTOR)) %>%
        dplyr::right_join(catch %>% dplyr::rename(POT = HAUL, N_ENTRIES = NUMBER_CRAB)) %>%
        dplyr::select(HAUL_ID, SPECIES_CODE, CRUISE, VESSEL, POT, NUMBER_CRAB, N_ENTRIES) -> catch_summary
      
  # Process catch_summary table for Oracle, save
      catch_summary %>%
        dplyr::select(!N_ENTRIES) %>%
        dplyr::rename(HAUL = POT) %>%
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
      
      
  # Calculate COUNT and CPUE per pot 
      mat_spec %>%
        dplyr::mutate(CPUE = SAMPLING_FACTOR/SOAK_TIME) %>% 
        dplyr::group_by(VESSEL, POT, BUOY, LAT_DD, LON_DD, MAT_SEX) %>%
        dplyr::reframe(COUNT = sum(SAMPLING_FACTOR),
                         CPUE = sum(CPUE)) -> positive_pot_cpue
      
      
  # Expand potlifts file to all mat-sex categories and potlifts, join to positive catch file to get zeros 
      mat_sex_combos <- c("Mature male", "Immature male", "Mature female", "Immature female", "Legal male", "Sublegal male")
      
      positive_pot_cpue %>%
        dplyr::right_join(expand_grid(MAT_SEX = mat_sex_combos,
                               potlifts)) %>%
        replace_na(list(COUNT = 0, CPUE = 0)) %>%
        dplyr::select(VESSEL, POT, BUOY, LAT_DD, LON_DD, DATE_SET, TIME_SET, DATE_HAUL, TIME_HAUL, SOAK_TIME,
                      MAT_SEX, COUNT, CPUE)-> pot_cpue
      
  # Save csv
      write.csv(pot_cpue, "./Output/W2023_potcpue.csv")
      
  # Transform catch and tagging data to correct crs
      pot_cpue %>%
        sf::st_as_sf(coords = c(x = "LON_DD", y = "LAT_DD"), crs = sf::st_crs(4326)) %>%
        sf::st_transform(crs = map.crs)-> pot_cpue_mapdat 
      
      tagging %>%
        sf::st_as_sf(coords = c(x = "LON_DD", y = "LAT_DD"), crs = sf::st_crs(4326)) %>%
        sf::st_transform(crs = map.crs) %>%
        dplyr::mutate(VESSEL = ifelse(VESSEL == "162", "Summer Bay", "Silver Spray"),
                      MAT_SEX = "Mature male",
                      shp = case_when((VESSEL == "Silver Spray" & N > 1) ~ 15,
                                      (VESSEL == "Silver Spray" & N == 1) ~ 0,
                                      (VESSEL == "Summer Bay" & N > 1) ~ 25, 
                                      (VESSEL == "Summer Bay" & N == 1) ~ 6)) -> tagging_mapdat 
      
      
  #set up plotting features
      map_layers <- akgfmaps::get_base_layers(select.region = "bs.south", set.crs="auto") # get map layers
      
      
      plot.boundary.untrans <- data.frame(y = c(54.5, 58.5), 
                                          x = c(-164.8, -159)) 
      
      plot.boundary <- plot.boundary.untrans %>%
                          akgfmaps::transform_data_frame_crs(out.crs = map.crs) # specify plot boundary
      
      
      breaks.x <- map_layers$lon.breaks[(map_layers$lon.breaks >= plot.boundary.untrans$x[1] &  # set lon breaks
                                           map_layers$lon.breaks < plot.boundary.untrans$x[2]) == TRUE]
      
      breaks.y <- map_layers$lat.breaks[(map_layers$lat.breaks > plot.boundary.untrans$y[1] & # set lat breaks
                                           map_layers$lat.breaks < plot.boundary.untrans$y[2]) == TRUE]
      
      max.date <- max(pot_cpue_mapdat$DATE_HAUL) # label for most recent pot haul date
      
      pal <- viridis::mako(10) # set palette
      
  # Plot
      mat_sex_combos %>%
        purrr::map(~ggplot() +
                    #geom_sf(data = map_layers$bathymetry, color=alpha("grey70")) +
                    geom_sf(data = st_as_sf(BB_strata), fill = NA, color = "black", linewidth = 1) +
                    geom_sf(data = st_as_sf(RKCSA_sub), fill = NA, color = "red", alpha= 0.5, linewidth = 1) +
                    geom_sf(data = st_as_sf(RKCSA), fill = NA, color = "red", alpha =0.5, linewidth = 1) +
                    geom_sf(data = map_layers$akland, fill = "grey80") +
                    geom_sf(data = filter(pot_cpue_mapdat, MAT_SEX == .x),
                             mapping = aes(size=COUNT, fill = COUNT), shape = 21, colour = "black", stat="identity", position="identity")+
                    geom_sf(data = filter(tagging_mapdat, MAT_SEX == .x),
                             mapping = aes(shape = as.factor(shp)), size= 2.5, stat="identity", position="identity")+
                    scale_shape_manual(values = c(0, 2, 15, 17), 
                                       labels = c("Silver Spray (n=1)", "Summer Bay (n=1)", "Silver Spray (n>1)", "Summer Bay (n>1)"))+
                    scale_size_continuous(range = c(2, 10), limits = c(0, max(filter(pot_cpue_mapdat, MAT_SEX == .x)$COUNT)))+ 
                    scale_fill_gradientn(limits = c(0, max(filter(pot_cpue_mapdat, MAT_SEX == .x)$COUNT)), 
                                          colors = c("gray", rev(pal[5:length(pal)])))+
                    scale_x_continuous(breaks = breaks.x, labels = paste0(breaks.x, "°W")) + 
                    scale_y_continuous(breaks = breaks.y, labels = paste0(breaks.y, "°N")) +
                    ggtitle(paste("BBRKC", .x)) +
                    guides(size = guide_legend(title.position = "top"), 
                            fill = guide_legend(), shape = guide_legend(title = "TAG RELEASE", title.position = "top", nrow =2))+
                    coord_sf(xlim = plot.boundary$x,
                              ylim = plot.boundary$y) +
                    geom_text(data = akgfmaps::transform_data_frame_crs(
                                data.frame(label = paste("Last haul date:",max.date),x = c(-160.5), y = c(58.3)), 
                                out.crs = map.crs), 
                                mapping = aes(x = x, y = y, label = label, fontface = "bold"), 
                                size = 3.5, color = c("#40498EFF"))+
                    theme_bw() +
                    theme(axis.title = element_blank(),
                           axis.text = element_text(size = 10),
                           legend.text = element_text(size = 10),
                           legend.title = element_text(size = 10),
                           legend.position = "bottom",
                           legend.direction = "horizontal",
                           plot.title = element_text(face = "bold"))) -> BBRKC.maps
      
  # Save plots
      ggsave(plot = BBRKC.maps[[1]], "./Figures/BBRKC_matmale.png", height=7, width=10, units="in")
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
        if(unique(specimen_table$VESSEL) %in% c(162, 164) == FALSE){
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
        if(unique(filter(specimen_table, SEX != 1)$EGG_COLOR %in% c(0, 2:6)) == FALSE){
          print("ERROR: invalid female egg color code (not 0 or 2:6)")
        }
        
        #6) Are egg condition codes valid for females?"
        if(unique(filter(specimen_table, SEX != 1)$EGG_CONDITION %in% c(0:4) == FALSE)){
          print("ERROR: invalid female egg condition code (not 0:4)")
        } 
        
        # 7) Are clutch size codes valid for females?
        if(unique(filter(specimen_table, SEX != 1)$CLUTCH_SIZE %in% c(0:6)) == FALSE){
          print("ERROR: invalid female clutch size code (not 0:6)")
        }
        
        # 8) Any egg, egg condition, or clutch size codes assigned to males?
        if(unique((specimen_table$SEX == 1 & is.na(specimen_table$EGG_COLOR | specimen_table$EGG_CONDITION | 
                                                   specimen_table$CLUTCH_SIZE) == FALSE))== TRUE){
          print("ERROR: egg, egg condition, or clutch size code assigned to male")
        } 
        
        
        # 9) Any questionable egg condition x shell condition combinations for females?"
        
        # Checking shell condition = 0 and egg condition = 1
        if(unique(filter(specimen_table, SEX != 1)$SHELL_CONDITION == 0 & 
                  filter(specimen_table, SEX != 1)$EGG_CONDITION == 1) == TRUE){
          print("ERROR: female with shell condition = 0 and egg condition = 1")
        } 
        
        # Checking shell condition = 1 and egg condition >1
        if(unique(filter(specimen_table, SEX != 1)$SHELL_CONDITION == 1 & 
                  filter(specimen_table, SEX != 1)$EGG_CONDITION > 1) == TRUE){
          print("ERROR: female with shell condition = 1 and egg condition >1")
        }
        
        # Checking shell condition = 3, 4, or 5 and egg condition = 1
        if(unique(filter(specimen_table, SEX != 1)$SHELL_CONDITION %in% c(3:5) & 
                  filter(specimen_table, SEX != 1)$EGG_CONDITION == 1) == TRUE){
          print("ERROR: female with shell condition 3:5 and egg condition = 1")
        }
        
        # Checking shell condition = 1 and egg condition >=2
        if(unique(filter(specimen_table, SEX != 1)$SHELL_CONDITION ==1 & 
                  filter(specimen_table, SEX != 1)$EGG_CONDITION >=2) == TRUE){
          print("ERROR: female with shell condition = 1 and egg condition >=2")
        }
        
        #10) Any females without egg color, egg condition, or clutch codes?
        if(unique(specimen_table$SEX == 2 & (is.na(specimen_table$EGG_COLOR | specimen_table$EGG_CONDITION 
                                                   | specimen_table$CLUTCH_SIZE) == TRUE)) == TRUE){
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
          dplyr::select(VESSEL, POT, BUOY, LAT_DD, LON_DD) %>%
          distinct() %>%
          as.data.frame()-> xx
        
        potlifts %>%
          dplyr::select(VESSEL, POT, BUOY, LAT_DD, LON_DD) %>%
          distinct() -> yy
        
        if(TRUE %in% is.na(suppressMessages(right_join(xx, yy, keep = TRUE))) == TRUE){
          print("ERROR: pot, station, and/or buoy IDs do not match between potlifts table and pot cpue table")
        }
        
        print("COMPARING SPECIMEN TABLE WITH CATCH SUMMARY...")
        specimen_table %>%
          dplyr::group_by(CRUISE, VESSEL, POT, SPECIES_CODE) %>%
          dplyr::reframe(NUMBER_CRAB = sum(SAMPLING_FACTOR),
                         N_ENTRIES = n()) -> spec_sum
        
        catch_summary %>%
          dplyr::select(CRUISE, VESSEL, POT, SPECIES_CODE, NUMBER_CRAB, N_ENTRIES) -> catch_sum
        
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
        map_layers <- suppressWarnings(akgfmaps::get_base_layers(select.region = "bs.south", set.crs="auto")) # get map layers
        
        plot.boundary.untrans <- data.frame(y = c(54, 59.5), 
                                            x = c(-168, -158)) 
        
        plot.boundary <- plot.boundary.untrans %>%
          akgfmaps::transform_data_frame_crs(out.crs = map.crs) # specify plot boundary
        
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
      