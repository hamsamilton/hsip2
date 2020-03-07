library(dplyr)

# Aggregating measurements
bp <- read.csv("/Users/User/Box Sync/Projects/Mimic_HSIP/Mimic_Data/Cardiac_Arrest/bp_96.csv")
temp <- read.csv("/Users/User/Box Sync/Projects/Mimic_HSIP/Mimic_Data/Cardiac_Arrest/temp_96.csv")
spo2 <- read.csv("/Users/User/Box Sync/Projects/Mimic_HSIP/Mimic_Data/Cardiac_Arrest/spo2_96.csv")
rr <- read.csv("/Users/User/Box Sync/Projects/Mimic_HSIP/Mimic_Data/Cardiac_Arrest/rr_96.csv")
hr <- read.csv("/Users/User/Box Sync/Projects/Mimic_HSIP/Mimic_Data/Cardiac_Arrest/hr_96.csv")
demo_graphics <- read.csv("/Users/User/Box Sync/Projects/Mimic_HSIP/Mimic_Data/Cardiac_Arrest/ca_demographics.csv")

# join by id
ds <- full_join(demo_graphics, bp, by = "subject_id")
ds <- full_join(ds, temp, by = "subject_id")
ds <- full_join(ds, spo2, by = "subject_id")
ds <- full_join(ds, rr, by = "subject_id")
ds <- full_join(ds, hr, by = "subject_id")

# consider min temps below 93 as hypopthermic cooling
ds$cooling <- ifelse(ds$min_temp <= 93, 1, 0)
table(ds$cooling)

# full data
#ds_full <- na.omit(ds)
#table(ds_full$mortality)

# export final datasets
write.csv(ds, "/Users/User/Box Sync/Projects/Mimic_HSIP/Mimic_Data/Cardiac_Arrest/To_Analyze/vars_96_Hours.csv", row.names = FALSE)
