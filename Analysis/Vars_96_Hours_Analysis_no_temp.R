## Analysis
## Meg Hutch
## 03.07.2020 - 96 Hours - no temp

library(dplyr)
library(tableone)
library(car)
library(ggplot2)

ds <- read.csv("/Users/User/Box Sync/Projects/Mimic_HSIP/Mimic_Data/Cardiac_Arrest/To_Analyze/vars_96_Hours.csv")

# remove temps
ds <- ds %>% select(-cooling, -min_temp, -mean_temp, -max_temp, -sd_temp, 
                    -var_temp, -entropy_temp)

# only keep full observations
ds_full <- na.omit(ds)
table(ds_full$mortality) # we can have up to 5 predictors 

##########################################################################3
# Univariate on each entropy predictor separately
summary(glm(mortality ~ entropy_bp, family = "binomial", data = ds_full))
summary(glm(mortality ~ entropy_hr, family = "binomial", data = ds_full))
summary(glm(mortality ~ entropy_spo2, family = "binomial", data = ds_full))
summary(glm(mortality ~ entropy_rr, family = "binomial", data = ds_full))

#####################################################################
## Purposeful model selection; only patients with full data
varList <- names(ds_full[, 3:28])

## Step 1 - Univariate analysis
tab1 <- CreateTableOne(vars=varList, strata = "mortality", data=ds_full,)
print(tab1, showAllLevels = TRUE) 

# Keep only variabels < 0.250 from univariate analysis
#or if the variable is known to be clinically important
mod1 <- glm(mortality ~ gender + age + entropy_spo2 + max_rr + sd_rr + 
              var_rr, family = "binomial", data = ds_full)

#use p < .05 or p < 0.10 to decide on importance
summary(mod1)
Anova(mod1)

mod2 <- glm(mortality ~  sd_rr + var_rr, family = "binomial", data = ds_full)
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

# All coefficeints < 20% change
# mod 2 now considered mod 3
mod3 <- mod2

# add back each non-significant variable one at a time
# This checks for variables that are not by themselves associated
# with outcome, but may be important in the presence of other
# variables
summary(update(mod3, ~. + min_bp))$coef[4, ]
summary(update(mod3, ~. + max_bp))$coef[4, ]
summary(update(mod3, ~. + mean_bp))$coef[4, ]
summary(update(mod3, ~. + sd_bp))$coef[4, ]
summary(update(mod3, ~. + var_bp))$coef[4, ]
summary(update(mod3, ~. + entropy_bp))$coef[4, ]
summary(update(mod3, ~. + min_spo2))$coef[4, ]
summary(update(mod3, ~. + max_spo2))$coef[4, ]
summary(update(mod3, ~. + mean_spo2))$coef[4, ]
summary(update(mod3, ~. + sd_spo2))$coef[4, ]
summary(update(mod3, ~. + var_spo2))$coef[4, ]
summary(update(mod3, ~. + min_rr))$coef[4, ]
summary(update(mod3, ~. + mean_rr))$coef[4, ]
summary(update(mod3, ~. + entropy_rr))$coef[4, ]
summary(update(mod3, ~. + min_hr))$coef[4, ]
summary(update(mod3, ~. + max_hr))$coef[4, ] # 0.099
summary(update(mod3, ~. + mean_hr))$coef[4, ]
summary(update(mod3, ~. + sd_hr))$coef[4, ]
summary(update(mod3, ~. + var_hr))$coef[4, ]
summary(update(mod3, ~. + entropy_hr))$coef[4, ]


# Save model 3 as model 4 - this is our preliminary main effects model
mod4 <- mod3

# Step 5: Examine more closely the variables in Model 4
# For each continuous variable, check whether logit increases or
# decreases linearly, i.e. the linearity assumption 
# check whether direction of the effect makes sense


logit <- function(pi) log(pi/(1-pi))

yhat_sd_rr <- loess(mortality ~ sd_rr, data=ds_full) %>% predict
yhat_var_rr <- loess(mortality ~ var_rr, data=ds_full) %>% predict

# Check for linearity - none of these are linear
qplot(ds_full$sd_rr, logit(yhat_sd_rr)) + geom_line() + theme_bw(base_size = 20)
qplot(ds_full$var_rr, logit(yhat_var_rr)) + geom_line() + theme_bw(base_size = 20)

# We have violated linearity assumptions. Convert to categorical variable 
ds_full$isd_rr <- cut(ds_full$sd_rr,  
                             breaks = c(0, 3.7, 10))
ds_full$ivar_rr <- cut(ds_full$var_rr,  
                      breaks = c(0, 15, 85))

# remodel with the new categorical variables
mod4.2 <- glm(mortality ~ isd_rr + ivar_rr, family = "binomial", 
              data = ds_full)
summary(mod4.2)

mod5 <- mod4.2

# This is now the main effects model

## Step 6: Consider Interactions
# Any interaction considered should make sense from a clinical perspective
# Add one interaction at a time (a single interaction may
# involve multiple terms if one of the variables is categorical)
# I Use LRT to test for interactions, and only include interactions
# with p<.05 (or even <.01)

f.upper <- ~ gender + age + min_bp + max_bp + mean_bp + 
  sd_bp + var_bp + entropy_bp + min_spo2 + max_spo2 + mean_spo2 + sd_spo2 + 
  var_spo2 + entropy_spo2 + mean_spo2 + sd_spo2 + var_spo2 + entropy_spo2 + 
  min_rr + max_rr + mean_rr + sd_rr + var_rr + entropy_rr + min_hr + max_hr +
  mean_hr + sd_hr + var_hr + entropy_hr
f.upper2 <- update(f.upper, ~ (.)^2)
f.upper2

summary(update(mod5, ~ . + isd_rr:age))
summary(update(mod5, ~ . + ivar_rr:age)) # significant
summary(update(mod5, ~ . + min_bp:age))
summary(update(mod5, ~ . + max_bp:age))
summary(update(mod5, ~ . + mean_bp:age))
summary(update(mod5, ~ . + sd_bp:age))
summary(update(mod5, ~ . + var_bp:age))
summary(update(mod5, ~ . + min_spo2:age))
summary(update(mod5, ~ . + max_spo2:age))
summary(update(mod5, ~ . + mean_spo2:age))
summary(update(mod5, ~ . + sd_spo2:age))
summary(update(mod5, ~ . + var_spo2:age))
summary(update(mod5, ~ . + mean_spo2:age))
summary(update(mod5, ~ . + sd_spo2:age))
summary(update(mod5, ~ . + var_spo2:age))
summary(update(mod5, ~ . + min_rr:age))
summary(update(mod5, ~ . + max_rr:age))
summary(update(mod5, ~ . + mean_rr:age)) # < 0.10
summary(update(mod5, ~ . + min_hr:age))
summary(update(mod5, ~ . + max_hr:age))
summary(update(mod5, ~ . + mean_hr:age))
summary(update(mod5, ~ . + sd_hr:age))
summary(update(mod5, ~ . + var_hr:age))
summary(update(mod5, ~ . + entropy_bp:age))
summary(update(mod5, ~ . + entropy_rr:age))
summary(update(mod5, ~ . + entropy_hr:age))
summary(update(mod5, ~ . + entropy_spo2:age))

# gender
summary(update(mod5, ~ . + isd_rr:gender))
summary(update(mod5, ~ . + ivar_rr:gender)) 
summary(update(mod5, ~ . + min_bp:gender))
summary(update(mod5, ~ . + max_bp:gender))
summary(update(mod5, ~ . + mean_bp:gender))
summary(update(mod5, ~ . + sd_bp:gender))
summary(update(mod5, ~ . + var_bp:gender))
summary(update(mod5, ~ . + min_spo2:gender))
summary(update(mod5, ~ . + max_spo2:gender))
summary(update(mod5, ~ . + mean_spo2:gender))
summary(update(mod5, ~ . + sd_spo2:gender))
summary(update(mod5, ~ . + var_spo2:gender))
summary(update(mod5, ~ . + mean_spo2:gender))
summary(update(mod5, ~ . + sd_spo2:gender))
summary(update(mod5, ~ . + var_spo2:gender))
summary(update(mod5, ~ . + min_rr:gender))
summary(update(mod5, ~ . + max_rr:gender))
summary(update(mod5, ~ . + mean_rr:gender)) 
summary(update(mod5, ~ . + min_hr:gender))
summary(update(mod5, ~ . + max_hr:gender))
summary(update(mod5, ~ . + mean_hr:gender))
summary(update(mod5, ~ . + sd_hr:gender))
summary(update(mod5, ~ . + var_hr:gender))
summary(update(mod5, ~ . + entropy_bp:gender))
summary(update(mod5, ~ . + entropy_rr:gender))
summary(update(mod5, ~ . + entropy_hr:gender))
summary(update(mod5, ~ . + entropy_spo2:gender)) # < 0.10

# Entropy interactions
summary(update(mod5, ~ . + entropy_bp:entropy_rr))
summary(update(mod5, ~ . + entropy_spo2:entropy_rr))
summary(update(mod5, ~ . + entropy_hr:entropy_rr))

summary(update(mod5, ~ . + entropy_spo2:entropy_bp))
summary(update(mod5, ~ . + entropy_hr:entropy_bp))

summary(update(mod5, ~ . + entropy_hr:entropy_spo2))


# Add all interactions that are found significant when added one at a time 
# to Model 5. Perform Step 2 (remove non-significant interactions, and refit
# the model). 
# Do not remove main effects, only interactions
# This is the preliminary final model (Model 6)

mod5.1 <- glm(mortality ~ isd_rr + ivar_rr + ivar_rr:age + mean_rr:age +
                entropy_spo2:gender, family = "binomial", data = ds_full)

anova(mod5, mod5.1, test = "LRT") # mod5.1 is almost significant better than mod5


# We should check final model fit 

# check linearity

# check residuals 