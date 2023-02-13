library(tidyverse)
test<-function(filename) {
  sample_data<-read_csv(filename)
  data<-sample_data %>% 
    select(CATCH_SAMPLE_ID,SPECIES_CODE,SEX,SAMPLE_MODIFIER_ID)
}

test_run<-test("C:/Users/jon.richar/Work/GitRepos/WinterSurvey2023/testData/SimulatedHaul/RawData/SAM_79_HAUL0999_SAMPLE_02081040.csv")

test_run
