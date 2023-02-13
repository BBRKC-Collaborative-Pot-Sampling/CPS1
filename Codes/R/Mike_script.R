
library(tidyverse)
sample<- list.files("./testData/EBS_Data/SAMPLE/")
sample_values<- list.files("./testData/EBS_Data/SAMPLE_VALUES")
setwd("C:/Users/jon.richar/Work/GitRepos/WinterSurvey2023/testData/EBS_Data/SAMPLE")

df <-list.files(path = "./testData/EBS_Data/SAMPLE/") %>% 
map_df(~read_csv(.))
df



library(data.table)
df <- 
  list.files(path = "./testData/EBS_Data/SAMPLE/", pattern = "*.csv") %>% 
  map_df(~fread(.))
df

















dates <- temps  <- matrix()
experiment.file <- data.frame()
sample<- list.files("./testData/EBS_Data/SAMPLE/")
sample_values<- list.files("./testData/EBS_Data/SAMPLE_VALUES")
length(sample)
sample

## begin with full North Pacific grid --------------------------

# create blank df to hold time series of sst and anomalies wrt 1950-1999
sample_series<- sample_value_series <- data.frame()

setwd()
for(i in 1:length(sample)){ # start i loop (each CMIP6 model)
  i <- 1
  path <- paste("./testData/EBS_Data/SAMPLE/", sample[i], sep="")
  x<-read_csv(sample[i])
}
path
read_csv(sample[1])
----------------------------------


