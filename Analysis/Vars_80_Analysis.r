## Analysis
## Meg Hutch
## 03.06.2020 - 80 window size

ds <- read.csv("/Users/User/Box Sync/Projects/Mimic_HSIP/Mimic_Data/Cardiac_Arrest/To_Analyze/vars_80.csv")

# Univariate Entropy
lm1 <- glm(mortality ~ entropy_bp, family = "binomial", data = ds)
summary(lm1)

lm2 <- glm(mortality ~ entropy_rr, family = "binomial", data = ds)
summary(lm2)

lm3 <- glm(mortality ~ entropy_temp, family = "binomial", data = ds)
summary(lm3)

lm4 <- glm(mortality ~ entropy_hr, family = "binomial", data = ds)
summary(lm4)

# remove the outlier to avoid error message
#lm5 <- glm(mortality ~ entropy_spo2, family = "binomial", data = ds)
lm5 <- glm(mortality ~ entropy_spo2, family = "binomial", data = ds %>% filter(!subject_id == 3139))
summary(lm5)

# all entropy values
lm6 <- glm(mortality ~ entropy_bp + entropy_rr + entropy_temp + 
             entropy_hr + entropy_spo2, family = "binomial", 
           data = ds %>% filter(!subject_id == 3139))
summary(lm6)

# no temp - small sample
lm7 <- glm(mortality ~ entropy_bp + entropy_rr + 
             entropy_hr + entropy_spo2, family = "binomial", 
           data = ds %>% filter(!subject_id == 3139))
summary(lm7)

# control for cooling only with temp - no interactions
lm7 <- glm(mortality ~ entropy_bp + entropy_rr + entropy_temp:cooling + 
             entropy_hr + entropy_spo2, family = "binomial", 
           data = ds %>% filter(!subject_id == 3139))
summary(lm7)
