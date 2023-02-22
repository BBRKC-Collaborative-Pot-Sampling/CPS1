# PURPOSE -----------------
# 1) To automate accessing the raw data folder and haul files (by pot) from the winter survey
# - Each pot is identified by the four digit number after "HAUL" and there can be potentially up to 20 pots a day
# - Each file name is appended by a time/date stamp. If multiple versions are sent of a file, the most recent one
#   needs to be read

# LOAD PACKAGES -------------------------------------------------------------------------------------------------------------
library(tidyverse)
library(dplyr)
library(akgfmaps)
library(gsubfn)
library(terra)
library(rgdal)

# PROCESSING AND PLOTTING FOR RKC -------------------------------------------------------------------------------------------
  # Read in catch data 
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
      dplyr::select(HAUL_ID, SPECIMEN_ID, CATCH_SAMPLE_ID, SPECIES_CODE)
    
    raw_specimen_bio<- list.files("./Data/FTP/RawData - FTP/", pattern = "_SPECIMEN_BIOMETRICS") %>% 
      purrr::map_df(~read.csv(paste0("./Data/FTP/RawData - FTP/", .x))) 
    
    potlifts <- list.files("./Data/FTP/Supplemental - FTP/", pattern = "POTLIFTS") %>% #MAY NEED TO CHANGE
      purrr::map_df(~read.csv(paste0("./Data/FTP/Supplemental - FTP/", .x)))
  
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
  
  # Join raw_sample_values and raw_sample to get # tossed by sample modifier (e.g., "All sizes", "Immature") 
    samples <- right_join(raw_sample, raw_sample_values) %>%
      dplyr::select(HAUL_ID, CATCH_SAMPLE_ID, SPECIES_CODE, SPECIES_NAME, SEX, TOSSED)
  
  # Expand specimen biometric table, join to raw_specimen table 
    raw_specimen_bio %>%
      dplyr::select(HAUL_ID, SPECIMEN_ID, BIOMETRIC_NAME, VALUE) %>%
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
  
  
  # Calculate COUNT and CPUE per pot 
    mat_spec %>%
      # Calculate soak time, divide by sampling factor by soak time
      dplyr::mutate(SET_TIME = as.POSIXct(SET_TIME, format = "%m/%d/%Y %H:%M"),
                    PULL_TIME = as.POSIXct(PULL_TIME, format = "%m/%d/%Y %H:%M"),
                    SOAK_TIME = as.numeric(difftime(PULL_TIME, SET_TIME, units = "mins")),
                    CPUE = SAMPLING_FACTOR/SOAK_TIME) %>% #TOTAL CATCH
      dplyr::group_by(VESSEL, POT, FLOAT, LATITUDE, LONGITUDE, MAT_SEX, COMMON_NAME) %>%
      dplyr::summarise(COUNT = sum(SAMPLING_FACTOR),
                       CPUE = sum(CPUE)) -> positive_pot_cpue
  
  
  # Expand potlifts file to all mat-sex categories and potlifts, join to positive catch file to get zeros 
    mat_sex_combos <- c("Mature male", "Immature male", "Mature female", "Immature female", "Legal male", "Sublegal male")
    
    positive_pot_cpue %>%
      right_join(expand_grid(MAT_SEX = mat_sex_combos,
                             potlifts)) %>%
      replace_na(list(COUNT = 0, CPUE = 0)) -> pot_cpue
  
  
  # Map catch for Bristol Bay 
    
    # Transform data to correct crs, add new column for zero catch data for mapping purposes
    pot_cpue %>%
      #filter(MAT_SEX %in% c("Legal male", "Sublegal male", "Mature female", "Immature female")) %>%
      #mutate(MAT_SEX = factor(MAT_SEX, levels = c("Legal male", "Sublegal male", "Mature female", "Immature female")),
             #nonzero = ifelse(COUNT == 0, "N", "Y")) %>%
      mutate(nonzero = ifelse(COUNT == 0, "N", "Y")) %>%
      sf::st_as_sf(coords = c(x = "LONGITUDE", y = "LATITUDE"), 
                   crs = sf::st_crs(4326)) %>% 
      sf::st_transform(crs = map.crs) -> pot_cpue_mapdat #make sure data is in correct crs
    
    #set up plotting features
    map_layers <- akgfmaps::get_base_layers(select.region = "bs.south", set.crs="auto") # get map layers
    
    plot.boundary <- data.frame(y = c(55.5, 58.5), 
                                x = c(-165, -159.5)) %>%
                                akgfmaps::transform_data_frame_crs(out.crs = map.crs) # specify plot boundary
    
    max.date <- strapplyc(max(pot_cpue$PULL_TIME), "\\d+/\\d+/\\d+", simplify = TRUE) # label for most recent pot pull date
    
    breaks <- pretty(pot_cpue_mapdat$COUNT) # calculate plotting breaks
    limits <- c(0, max(pot_cpue_mapdat$COUNT)) # set plotting limits
    
    # Plot
    mat_sex_combos %>%
      purrr::map(~ggplot(filter(pot_cpue_mapdat, MAT_SEX == .x)) +
          #geom_sf(data = map_layers$bathymetry, color=alpha("grey70")) +
          geom_sf(data = st_as_sf(BB_strata), fill = NA, color = "black", linewidth = 1) +
          geom_sf(data = st_as_sf(RKCSA_sub), fill = NA, color = "red", linewidth = 1) +
          geom_sf(data = st_as_sf(RKCSA), fill = NA, color = "red", linewidth = 1) +
          geom_sf(data = map_layers$akland, fill = "grey80") +
          geom_sf(mapping = aes(size=COUNT, shape = nonzero), stat="identity", position="identity")+
          scale_size(range = c(2, 10), breaks = breaks, limits = limits)+
          scale_shape_manual(values = c(10, 19))+
          scale_x_continuous(breaks = map_layers$lon.breaks) + 
          scale_y_continuous(breaks = map_layers$lat.breaks) +
          #facet_wrap(~MAT_SEX, nrow = 2, ncol = 2) +
          ggtitle(paste("BBRKC", .x)) +
          guides(size = guide_legend(override.aes = list(shape =c(10, rep(19, (length(breaks)-2)))), 
                                       title = "COUNT"), shape = "none") +
          coord_sf(xlim = plot.boundary$x,
                   ylim = plot.boundary$y) +
          geom_text(data = akgfmaps::transform_data_frame_crs(
                       data.frame(label = paste("Last pull date:",max.date),
                       x = c(-160.5),
                       y = c(58.3)), 
                       out.crs = coldpool:::ebs_proj_crs), 
            mapping = aes(x = x, y = y, label = label, fontface = "bold"), size = 3.5, color = c("blue"))+
          theme_bw() +
          theme(#panel.border = element_rect(color = "black", fill = NA),
                #panel.background = element_rect(fill = NA, color = "black"),
                #legend.key = element_rect(fill = NA, color = "grey70"),
                #legend.position = c(0.90, 0.18),
                axis.title = element_blank(),
                axis.text = element_text(size = 10),
                legend.text = element_text(size = 10),
                legend.title = element_text(size = 10),
                #plot.background = element_rect(fill = NA, color = NA),
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
    
    
# PROCESSING AND PLOTTING FOR BYCATCH -------------------------------------------------------------------------------------------