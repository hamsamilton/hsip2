# Aggregating measurements
bp_200 <- read.csv("/Users/User/Box Sync/Projects/Mimic_HSIP/Mimic_Data/Cardiac_Arrest/bp_200.csv")
temp_200 <- read.csv("/Users/User/Box Sync/Projects/Mimic_HSIP/Mimic_Data/Cardiac_Arrest/temp_200.csv")
spo2_200 <- read.csv("/Users/User/Box Sync/Projects/Mimic_HSIP/Mimic_Data/Cardiac_Arrest/spo2_200.csv")
rr_200 <- read.csv("/Users/User/Box Sync/Projects/Mimic_HSIP/Mimic_Data/Cardiac_Arrest/rr_200.csv")
hr_200 <- read.csv("/Users/User/Box Sync/Projects/Mimic_HSIP/Mimic_Data/Cardiac_Arrest/hr_200.csv")
demo_graphics <- read.csv("/Users/User/Box Sync/Projects/Mimic_HSIP/Mimic_Data/Cardiac_Arrest/ca_demographics.csv")


# join by id
ds <- full_join(demo_graphics, bp_200, by = "subject_id")
ds <- full_join(ds, temp_200, by = "subject_id")
ds <- full_join(ds, spo2_200, by = "subject_id")
ds <- full_join(ds, rr_200, by = "subject_id")
ds <- full_join(ds, hr_200, by = "subject_id")

# export final datasets
write.csv(ds, "/Users/User/Box Sync/Projects/Mimic_HSIP/Mimic_Data/Cardiac_Arrest/To_Analyze/vars_200.csv", row.names = FALSE)
