library(dplyr)

# Aggregating measurements
bp_80 <- read.csv("/Users/User/Box Sync/Projects/Mimic_HSIP/Mimic_Data/Cardiac_Arrest/bp_80.csv")
temp_80 <- read.csv("/Users/User/Box Sync/Projects/Mimic_HSIP/Mimic_Data/Cardiac_Arrest/temp_80.csv")
spo2_80 <- read.csv("/Users/User/Box Sync/Projects/Mimic_HSIP/Mimic_Data/Cardiac_Arrest/spo2_80.csv")
rr_80 <- read.csv("/Users/User/Box Sync/Projects/Mimic_HSIP/Mimic_Data/Cardiac_Arrest/rr_80.csv")
hr_80 <- read.csv("/Users/User/Box Sync/Projects/Mimic_HSIP/Mimic_Data/Cardiac_Arrest/hr_80.csv")
demo_graphics <- read.csv("/Users/User/Box Sync/Projects/Mimic_HSIP/Mimic_Data/Cardiac_Arrest/ca_demographics.csv")

# join by id
ds <- full_join(demo_graphics, bp_80, by = "subject_id")
ds <- full_join(ds, temp_80, by = "subject_id")
ds <- full_join(ds, spo2_80, by = "subject_id")
ds <- full_join(ds, rr_80, by = "subject_id")
ds <- full_join(ds, hr_80, by = "subject_id")

# consider min temps below 93 as hypopthermic cooling
ds$cooling <- ifelse(ds$min_temp <= 93, 1, 0)

# full data
ds_full <- na.omit(ds)
table(ds_full$mortality)

# export final datasets
write.csv(ds, "/Users/User/Box Sync/Projects/Mimic_HSIP/Mimic_Data/Cardiac_Arrest/To_Analyze/vars_80.csv", row.names = FALSE)
