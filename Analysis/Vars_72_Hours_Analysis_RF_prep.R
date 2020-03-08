## Analysis
## Meg Hutch
## 03.08.2020 - 72 Hours - RF Data Prep

library(dplyr)
library(tableone)
library(car)
library(ggplot2)
library(caret)

ds <- read.csv("/Users/User/Box Sync/Projects/Mimic_HSIP/Mimic_Data/Cardiac_Arrest/To_Analyze/vars_72_Hours.csv")

#### Create training and test sets with and without including temperature

# Prep data without including temp
ds_no_temp <- ds %>% select(-cooling, -min_temp, -mean_temp, -max_temp, -sd_temp, 
                    -var_temp, -entropy_temp)

# only keep full observations
ds_no_temp <- na.omit(ds_no_temp)
table(ds_no_temp$mortality) # nicely split!

# keep temp
ds_full <- na.omit(ds)

table(ds_full$mortality) # small sample size

##########################################################################
# check distributions
table(ds_no_temp$mortality)

set.seed(0308)
trainIndex <- createDataPartition(ds_no_temp$mortality, p = .70, 
                                  list = FALSE, 
                                  times = 1)
head(trainIndex)

mcTrain <- ds_no_temp[trainIndex,]
mcTest  <- ds_no_temp[-trainIndex,]

table(mcTrain$mortality)
table(mcTest$mortality)

round(table(mcTrain$mortality)/80*100)
round(table(mcTest$mortality)/33*100)

# Sample size may be too small for it to be perfectly balanced

# save training and testing sets
write.csv(ds_no_temp, file = "C:/Users/User/Box Sync/Projects/Mimic_HSIP/To_Analyze/ds_no_temp_72.csv", row.names = FALSE)
write.csv(ds_full, file = "C:/Users/User/Box Sync/Projects/Mimic_HSIP/To_Analyze/ds_full.csv_72", row.names = FALSE)

write.csv(mcTrain, file = "C:/Users/User/Box Sync/Projects/Mimic_HSIP/To_Analyze/mcTrain_no_temp_72.csv", row.names = FALSE)
write.csv(mcTest, file = "C:/Users/User/Box Sync/Projects/Mimic_HSIP/To_Analyze/mcTest_no_temp_72.csv", row.names = FALSE)
