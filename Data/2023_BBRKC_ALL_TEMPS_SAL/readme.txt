READ ME

.readme: ADFG (CPS1_2023BB_RKC_LoggerData_ADFG.csv)
2023 Logger Data includes data from both the Silver Spray and Summer Bay. The Timestamp had to be adjusted by 1 hour to take
into account Daylight Savings Time. The first 3 columns are the recorded time from the loggers and the second 3 columns is the
time adjusted by 1 hour.
In order to bring the data together, a field for SPN2 was created so that the Silver Spray SPNs were numeric for non-core pots. 
The Summer Bay SPNs were all numeric, so the SPN2 field does not apply to this boat. Latitude and Longitude were also used to 
indicate a unique numeric identifier for a specific pot.
A field for Core Pots was also added to indicate whether the logger was in a survey pot or a non-survey pot.  This is indicated
by a Y for Core pot and N for Non-Core pot.
Core == Survey Pot
Non-Core == Tagging pot set at end of survey in a hotspot area, near locations of previously sampled high-catching survey pots  
(Silver Spray & Summer Bay) or camera deployment (Silver Spray, logger in camera pot or separate pot near camera)

.readme: ABSC (AMCC_FV_NRT_d113_5950_44d0.csv)
2023 ABSC Logger Data includes data from both the Silver Spray and Summer Bay.
Link to the ERDDAP database and you can query/download the data. (The data is labeled under Alaska Marine Conservation 
Council as they are the organization that purchased the sensors for our use) 
http://54.183.206.151/erddap/tabledap/AMCC_FV_NRT.html

.readme (Pot_IDs_and_Coordinates.csv)
PI L Zacher prepared grid survey points prepared prior to field research. J Weems added decimal degrees columns for mapping.

.readme (2023_BBRKC_ALL_TEMPS_SAL.csv)
Final export file for summerized bottom T/S data from 2023_BBRKC_AveTempSal_EDA.r script