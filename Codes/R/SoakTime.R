# Load packages
library(data.table)
library(lubridate)

# Create example data frame
dat <- fread("    Set.Time         Pull.Time
          2019-06-07 08:23    2019-06-07 07:53
          2019-06-07 08:53    2019-06-07 07:53
          2019-06-07 09:23    2019-06-07 07:53
          2019-06-07 09:53    2019-06-07 07:53
2019-06-07 10:23    2019-06-07 07:53
2019-06-07 10:53    2019-06-07 07:53
2019-06-07 11:53    2019-06-07 07:53
2019-06-07 13:53    2019-06-07 07:53
2019-06-07 15:53    2019-06-07 07:53
2019-06-07 17:53    2019-06-07 07:53")

# Convert to POSIXct class
dat <- dat[, lapply(.SD, function(x) ymd_hm(x))]

# See the class of each column
str(dat)


# Create a new column shows the time differences as hours
dat[, Time.After.Dose := as.double(Sample.Time - Dose.Time, units = "hours")]
print(dat)