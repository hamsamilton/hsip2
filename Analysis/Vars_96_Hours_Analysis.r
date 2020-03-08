## Analysis
## Meg Hutch
## 03.07.2020 - 96 Hours

library(dplyr)
library(tableone)
library(car)
library(ggplot2)

ds <- read.csv("/Users/User/Box Sync/Projects/Mimic_HSIP/Mimic_Data/Cardiac_Arrest/To_Analyze/vars_96_Hours.csv")

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
lm5 <- glm(mortality ~ entropy_spo2, family = "binomial", data = ds)
summary(lm5)

# all entropy values
lm6 <- glm(mortality ~ entropy_bp + entropy_rr + entropy_temp + 
             entropy_hr + entropy_spo2, family = "binomial", data = ds)
summary(lm6)

# no temp 
lm7 <- glm(mortality ~ entropy_bp + entropy_rr + 
             entropy_hr + entropy_spo2, family = "binomial", data = ds)
summary(lm7)

# control for cooling only with temp - no interactions
lm7 <- glm(mortality ~ entropy_bp + entropy_rr + entropy_temp:cooling + 
             entropy_hr + entropy_spo2, family = "binomial", data = ds)
summary(lm7)

lm8 <- glm(mortality ~ entropy_bp + entropy_rr + entropy_temp:cooling + 
             entropy_hr + entropy_spo2, family = "binomial", data = ds)
summary(lm8)

lm9 <- glm(mortality ~ entropy_bp + entropy_rr + entropy_hr + entropy_spo2 + 
max_bp + max_hr + min_rr + min_spo2, family = "binomial", data = ds)


## Entropy forward selection
ds_full <- na.omit(ds)
table(ds_full$mortality)


fit0 <- lm(mortality ~ 1, data=ds_full)

f.lower <- ~ 1
f.upper <- ~ entropy_bp + entropy_rr + entropy_hr + entropy_spo2 + entropy_temp
fit.stepF <- step(fit0, scope = f.upper, direction="forward")

lm1 <- glm(mortality ~ entropy_spo2 + entropy_temp, family = "binomial", 
           data = ds_full)
summary(lm1)

# Forward selection without temp
ds_full <- na.omit(ds)
fit0 <- lm(mortality ~ 1, data=ds_full)

f.lower <- ~ 1
f.upper <- ~ entropy_bp + entropy_rr + entropy_hr + entropy_spo2
fit.stepF <- step(fit0, scope = f.upper, direction="forward")

lm1 <- glm(mortality ~ entropy_spo2 + entropy_bp, family = "binomial", 
           data = ds_full)
summary(lm1)

# Forward selection control for cooling
ds_full <- na.omit(ds)
fit0 <- lm(mortality ~ 1, data=ds_full)

f.lower <- ~ 1
f.upper <- ~ entropy_bp + entropy_rr + entropy_hr + entropy_spo2 + entropy_temp + 
  cooling + entropy_temp:cooling
fit.stepF <- step(fit0, scope = f.upper, direction="forward")

lm1 <- glm(mortality ~ entropy_spo2 + entropy_temp, family = "binomial", 
           data = ds_full)
summary(lm1)



#####################################################################
## Purposeful model selection; only patients with full data
varList <- names(ds_full[, 3:35])

## Entropy forward selection
ds_full <- na.omit(ds)
table(ds_full$mortality)

# univariate analysis
tab1 <- CreateTableOne(vars=varList, strata = "mortality", data=ds_full,)
print(tab1, showAllLevels = TRUE) 

# Keep only variabels < 0.250 from univariate analysis
#or if the variable is known to be clinically important
mod1 <- glm(mortality ~ gender + age + mean_temp + var_temp + entropy_temp +
              entropy_spo2 + max_rr + sd_rr + var_rr + entropy_rr,
            family = "binomial", data = ds_full)

#use p < .05 or p < 0.10 to decide on importance
summary(mod1)
Anova(mod1)

mod2 <- glm(mortality ~ entropy_spo2 + gender, family = "binomial", data = ds_full)
summary(mod2)

# Compare Models
anova(mod1, mod2, test="LRT")

# determine whether coefficients are < 20% different
# if Beta > 20%, it means that one of the excluded variables is
# important in a sense that it provides a needed adjustment of
# the effect of the variables in the model
betaComp <- compareCoefs(mod1, mod2, se=FALSE)
betaAdj <- betaComp[,1]
betaUnadj <- betaComp[,2]
deltaBeta <- (betaAdj-betaUnadj)/betaAdj*100; deltaBeta

# Continue this process (Steps 2 and 3), one or two variables at a
#time, until all important "adjusting" variables are back in the
#model. Make sure excluded variables are not clinically important

mod3 <- glm(mortality ~ entropy_spo2 + gender + age, family = "binomial", data = ds_full)
summary(mod3)
betaComp <- compareCoefs(mod1, mod3, se=FALSE)
betaAdj <- betaComp[,1]
betaUnadj <- betaComp[,2]
deltaBeta <- (betaAdj-betaUnadj)/betaAdj*100; deltaBeta

mod4 <- glm(mortality ~ entropy_spo2 + gender + age + mean_temp, family = "binomial", data = ds_full)
summary(mod4)
betaComp <- compareCoefs(mod1, mod4, se=FALSE)
betaAdj <- betaComp[,1]
betaUnadj <- betaComp[,2]
deltaBeta <- (betaAdj-betaUnadj)/betaAdj*100; deltaBeta

mod5 <- glm(mortality ~ entropy_spo2 + gender + age + var_temp, family = "binomial", data = ds_full)
summary(mod5)
betaComp <- compareCoefs(mod1, mod5, se=FALSE)
betaAdj <- betaComp[,1]
betaUnadj <- betaComp[,2]
deltaBeta <- (betaAdj-betaUnadj)/betaAdj*100; deltaBeta

mod6 <- glm(mortality ~ entropy_spo2 + gender + age + entropy_temp, family = "binomial", data = ds_full)
summary(mod6)
betaComp <- compareCoefs(mod1, mod6, se=FALSE)
betaAdj <- betaComp[,1]
betaUnadj <- betaComp[,2]
deltaBeta <- (betaAdj-betaUnadj)/betaAdj*100; deltaBeta

mod7 <- glm(mortality ~ entropy_spo2 + gender + age + max_rr, family = "binomial", data = ds_full)
summary(mod7)
betaComp <- compareCoefs(mod1, mod7, se=FALSE)
betaAdj <- betaComp[,1]
betaUnadj <- betaComp[,2]
deltaBeta <- (betaAdj-betaUnadj)/betaAdj*100; deltaBeta

mod8 <- glm(mortality ~ entropy_spo2 + gender + age + sd_rr, family = "binomial", data = ds_full)
summary(mod8)
betaComp <- compareCoefs(mod1, mod8, se=FALSE)
betaAdj <- betaComp[,1]
betaUnadj <- betaComp[,2]
deltaBeta <- (betaAdj-betaUnadj)/betaAdj*100; deltaBeta

mod9 <- glm(mortality ~ entropy_spo2 + gender + age + var_rr, family = "binomial", data = ds_full)
summary(mod9)
betaComp <- compareCoefs(mod1, mod9, se=FALSE)
betaAdj <- betaComp[,1]
betaUnadj <- betaComp[,2]
deltaBeta <- (betaAdj-betaUnadj)/betaAdj*100; deltaBeta

mod10 <- glm(mortality ~ entropy_spo2 + gender + age + entropy_rr, family = "binomial", data = ds_full)
summary(mod10)
betaComp <- compareCoefs(mod1, mod10, se=FALSE)
betaAdj <- betaComp[,1]
betaUnadj <- betaComp[,2]
deltaBeta <- (betaAdj-betaUnadj)/betaAdj*100; deltaBeta

## It seems like age has a large affect, but this should be clinically important
# will use the smallest coefficient % change for age, and continue adding variables

mod4.1 <- glm(mortality ~ entropy_spo2 + gender + age + mean_temp + var_temp,
              family = "binomial", data = ds_full)
summary(mod4.1)
betaComp <- compareCoefs(mod1, mod4.1, se=FALSE)
betaAdj <- betaComp[,1]
betaUnadj <- betaComp[,2]
deltaBeta <- (betaAdj-betaUnadj)/betaAdj*100; deltaBeta


mod4.2 <- glm(mortality ~ entropy_spo2 + gender + age + mean_temp + 
                entropy_temp,
              family = "binomial", data = ds_full)
summary(mod4.2)
betaComp <- compareCoefs(mod1, mod4.2, se=FALSE)
betaAdj <- betaComp[,1]
betaUnadj <- betaComp[,2]
deltaBeta <- (betaAdj-betaUnadj)/betaAdj*100; deltaBeta

# Finally we get all the imporant variables back into the model with 
# coefficients < 20%

# mod 4.2 now considered mod 3
mod3 <- mod4.2

# add back each non-significant variable one at a time
# This checks for variables that are not by themselves associated
# with outcome, but may be important in the presence of other
# variables
summary(update(mod3, ~. + min_bp))$coef[7, ]
summary(update(mod3, ~. + max_bp))$coef[7, ]
summary(update(mod3, ~. + mean_bp))$coef[7, ]
summary(update(mod3, ~. + sd_bp))$coef[7, ]
summary(update(mod3, ~. + var_bp))$coef[7, ]
summary(update(mod3, ~. + entropy_bp))$coef[7, ]
summary(update(mod3, ~. + min_temp))$coef[7, ] # 0.09
summary(update(mod3, ~. + max_temp))$coef[7, ] # 0.055
summary(update(mod3, ~. + sd_temp))$coef[7, ]
summary(update(mod3, ~. + min_spo2))$coef[7, ]
summary(update(mod3, ~. + mean_spo2))$coef[7, ]
summary(update(mod3, ~. + sd_spo2))$coef[7, ]
summary(update(mod3, ~. + var_spo2))$coef[7, ]
summary(update(mod3, ~. + min_rr))$coef[7, ]
summary(update(mod3, ~. + mean_rr))$coef[7, ]
summary(update(mod3, ~. + min_hr))$coef[7, ]
summary(update(mod3, ~. + max_hr))$coef[7, ]
summary(update(mod3, ~. + mean_hr))$coef[7, ]
summary(update(mod3, ~. + sd_hr))$coef[7, ]
summary(update(mod3, ~. + var_hr))$coef[7, ]
summary(update(mod3, ~. + entropy_hr))$coef[7, ]
summary(update(mod3, ~. + cooling))$coef[7, ]

# add in max temp - almost significant
mod3.1 <- update(mod3, ~. + max_temp)
summary(mod3.1)

# check betas - change quite a bit! Don't include max_temp
betaComp <- compareCoefs(mod3, mod3.1, se=FALSE)
betaAdj <- betaComp[,1]
betaUnadj <- betaComp[,2]
deltaBeta <- (betaAdj-betaUnadj)/betaAdj*100; deltaBeta

# Save model 3 as model 4 - this is our preliminary main effects model
mod4 <- mod3

# Step 5: Examine more closely the variables in Model 4
# For each continuous variable, check whether logit increases or
# decreases linearly, i.e. the linearity assumption 
# check whether direction of the effect makes sense


logit <- function(pi) log(pi/(1-pi))

yhat_entropy_spo2 <- loess(mortality ~ entropy_spo2, data=ds_full) %>% predict
#yhat_gender <- loess(mortality ~ gender, data=ds_full) %>% predict
yhat_age <- loess(mortality ~ age, data=ds_full) %>% predict
yhat_mean_temp <- loess(mortality ~ mean_temp, data=ds_full) %>% predict
yhat_entropy_temp <- loess(mortality ~ entropy_temp, data=ds_full) %>% predict

# Check for linearity - none of these are linear
qplot(ds_full$entropy_spo2, logit(yhat_entropy_spo2)) + geom_line() + theme_bw(base_size = 20)
qplot(ds_full$age, logit(yhat_age)) + geom_line() + theme_bw(base_size = 20)
qplot(ds_full$mean_temp, logit(yhat_mean_temp)) + geom_line() + theme_bw(base_size = 20)
qplot(ds_full$entropy_temp, logit(yhat_entropy_temp)) + geom_line() + theme_bw(base_size = 20)

# We have violated linearity assumptions. Convert to categorical variable 
ds_full$ientropy_spo2 <- cut(ds_full$entropy_spo2,  
                             breaks = c(0, 0.33, 0.45, 0.9, 1.4))
ds_full$iage <- cut(ds_full$age,  
                             breaks = c(0, 58, 73, 90))
ds_full$imean_temp <- ifelse(ds_full$mean_temp <= 100, "low_temp", "high_temp")
ds_full$ientropy_temp <- ifelse(ds_full$entropy_temp <= 0.50, "low_temp_entropy", "high_temp_entropy")

# remodel with the new categorical variables
mod4.2 <- glm(mortality ~ ientropy_spo2 + iage + imean_temp + ientropy_temp +
                gender, family = "binomial", data = ds_full)
summary(mod4.2)

# This is now the main effects model

# We need to remove temperature, otherwise we can't have this many predictors
